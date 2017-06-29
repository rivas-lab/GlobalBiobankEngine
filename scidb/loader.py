import glob
import itertools
import logging
import numpy
import os
import scidbpy
import subprocess
import tempfile

try:
    from weakref import finalize
except ImportError:
    from backports.weakref import finalize

from future.moves.itertools import zip_longest


import config


logging.basicConfig(level=logging.INFO, format='%(funcName)s:%(message)s')
logger = logging.getLogger(__file__)


class Loader:
    def __init__(self):
        self.db = scidbpy.connect()
        self.fifo_names = Loader.make_fifos()
        self.instances = [str(i) for i in range(config.SCIDB_INSTANCE_NUM)]

        finalize(self,
                 self.remove_fifos,
                 tuple(self.fifo_names))

    def load_qc(self):
        qc = []
        for rec in config.QC_FILES:
            with open(rec['file']) as fp:
                lines = fp.readlines()
                skip = rec.get('header', 0)
                qc.extend(
                    l.strip().split()[0]
                    for l in lines[:-skip if skip else None])
                logger.info('QC:file:%s lines:%d skip:%d',
                            rec['file'], len(lines), skip)

        # TODO fix when upload_data=list() implemented
        self.db.input(upload_data=numpy.array(qc)).store(config.QC_ARRAY)
        logger.info('Array:%s records:%d', config.QC_ARRAY, len(qc))

    def load_icd_map(self):
        icd_lst = list(set(Loader.get_icd(*Loader.get_icd_parts(fn))
                for fn in itertools.chain(
                        glob.iglob(config.ICD_GLOB),
                        glob.iglob(config.QT_GLOB))))
        self.icd_id_map = dict(zip(icd_lst, range(len(icd_lst))))
        self.db.input(config.ICD_INDEX_SCHEMA,
                      upload_data=numpy.array(icd_lst)).store(
                          config.ICD_INDEX_ARRAY)

    def load_icd(self):
        self.db.create_array(config.ICD_ARRAY, config.ICD_SCHEMA)
        file_iter = [glob.iglob(config.ICD_GLOB)] * config.SCIDB_INSTANCE_NUM

        for file_names in zip_longest(*file_iter):

            # Last chunk might have None, strip them out
            if None in file_names:
                file_names = [file_name
                              for file_name in file_names if file_name]

            # Get conditional expressions
            (icd_id_cond, icdind_cond) = self.get_icd_cond(file_names)

            # Make pipes
            logger.info('Pipes:starting %d...', len(file_names))
            pipes = [Loader.make_pipe(file_name, fifo_name)
                     for (file_name, fifo_name) in zip(file_names,
                                                       self.fifo_names)]

            query = config.ICD_LOAD_QUERY.format(
                icd=config.ICD_ARRAY,
                qc=config.QC_ARRAY,
                paths=';'.join(self.fifo_names[:len(file_names)]),
                instances=';'.join(self.instances[:len(file_names)]),
                icd_id_cond=icd_id_cond,
                icdind_cond=icdind_cond)

            logger.info('Query:starting...')
            self.db.iquery(query)
            self.remove_versions(config.ICD_ARRAY)
            logger.info('Query:done')

            logger.info('Pipes:return code:%s',
                        ','.join(str(pipe.poll()) for pipe in pipes))

    def load_qt(self):
        file_iter = [glob.iglob(config.QT_GLOB)] * config.SCIDB_INSTANCE_NUM

        for file_names in zip_longest(*file_iter):

            # Last chunk might have None, strip them out
            if None in file_names:
                file_names = [file_name
                              for file_name in file_names if file_name]

            # Get conditional expressions
            (icd_id_cond, icdind_cond) = self.get_icd_cond(file_names)

            # Make pipes
            logger.info('Pipes:starting %d...', len(file_names))
            pipes = [Loader.make_pipe(file_name, fifo_name)
                     for (file_name, fifo_name) in zip(file_names,
                                                       self.fifo_names)]

            query = config.QT_LOAD_QUERY.format(
                icd=config.ICD_ARRAY,
                qc=config.QC_ARRAY,
                paths=';'.join(self.fifo_names[:len(file_names)]),
                instances=';'.join(self.instances[:len(file_names)]),
                icd_id_cond=icd_id_cond,
                icdind_cond=icdind_cond)

            logger.info('Query:starting...')
            self.db.iquery(query)
            self.remove_versions(config.ICD_ARRAY)
            logger.info('Query:done')

            logger.info('Pipes:return code:%s',
                        ','.join(str(pipe.poll()) for pipe in pipes))

    def load_icd_info(self):
        self.db.create_array(config.ICD_INFO_ARRAY, config.ICD_INFO_SCHEMA)
        self.db.iquery(config.ICD_INFO_LOAD_QUERY.format(
            path=config.ICD_INFO_FILE,
            icd_info=config.ICD_INFO_ARRAY,
            icd_index=config.ICD_INDEX_ARRAY))
        logger.info('Array:%s', config.ICD_INFO_ARRAY)

    def get_icd_cond(self, file_names):
        """Build SciDB conditional expressions ("iif") to map file names to
        "icd_id", "icdind" using "src_instance_id"

        """
        icd_id_cond = 'null'
        icdind_cond = 'null'
        for (inst, file_name) in zip(range(len(file_names)), file_names):
            (prefix, suffix, intadd) = Loader.get_icd_parts(file_name)

            icd_id = self.icd_id_map[Loader.get_icd(prefix, suffix, intadd)]
            icd_id_cond = ('iif(src_instance_id = {inst}, ' +
                           '{icd_id}, {prev})').format(
                               inst=inst, icd_id=icd_id, prev=icd_id_cond)

            icdind_cond = ('iif(src_instance_id = {inst}, ' +
                           '\'{val}\', {prev})').format(
                               inst=inst,
                               val=str(intadd) + suffix,
                               prev=icdind_cond)
        return (icd_id_cond, icdind_cond)

    def remove_arrays(self):
        if not Loader.confirm('Remove and recreate arrays'):
            return
        for array in (config.QC_ARRAY,
                      config.ICD_INDEX_ARRAY,
                      config.ICD_ARRAY,
                      config.ICD_INFO_ARRAY):
            try:
                self.db.remove(array)
            except:
                pass

    def remove_versions(self, name):
        self.db.remove_versions(name,
                                self.db.versions(name)[:]['version_id'].max())

    @staticmethod
    def get_icd_parts(file_name):
        name = os.path.basename(file_name)
        parts = name.split('.')

        # Figure out prefix and intadd
        if 'brainmri' in parts:
            prefix = 'BRMRI'
            intadd = 10
        elif 'additionalimaging' in parts:
            prefix = 'ADD'
            intadd = 20
        elif 'initialdata' in parts:
            prefix = 'INI'
            intadd = 30
        elif '_FH2' in name:
            prefix = 'FH'
            intadd = 40
        elif 'RH' in name:
            prefix = 'RH'
            intadd = 50
        elif 'cancer' in name:
            prefix = 'cancer'
            intadd = 60
        elif 'HC' in name:
            prefix = 'HC'
            intadd = 70

        # Figure out suffix
        suffix = parts[1].split(
                '_')[0].strip()[1:].strip(
                    'RH').strip(
                        'FH').strip(
                            'cancer').strip(
                                'HC').split(
                                    '_FH2')[0]

        return (prefix, suffix, intadd)

    @staticmethod
    def get_icd(prefix, suffix, intadd):
        return prefix + suffix

    @staticmethod
    def make_pipe(file_name, fifo_name):
        cmd = 'zcat {} > {}'.format(file_name, fifo_name)
        proc = subprocess.Popen(cmd, shell=True)

        logger.info('Spawn:%s pid:%d', cmd, proc.pid)
        return proc

    @staticmethod
    def make_fifo():
        fifo_dir = tempfile.mkdtemp()
        fifo_name = os.path.join(fifo_dir, 'fifo')
        os.mkfifo(fifo_name)

        logger.info('FIFO:%s', fifo_name)
        return fifo_name

    @staticmethod
    def remove_fifo(fifo_name):
        logger.info('Remove:%s', fifo_name)
        os.remove(fifo_name)
        os.rmdir(os.path.dirname(fifo_name))

    @staticmethod
    def make_fifos():
        return [Loader.make_fifo() for i in range(config.SCIDB_INSTANCE_NUM)]

    @staticmethod
    def remove_fifos(fifo_names):
        for fn in fifo_names:
            Loader.remove_fifo(fn)

    @staticmethod
    def confirm(prompt=None):
        ans = raw_input('{}{} confirm with "y":'.format(
            prompt, ',' if prompt else ''))
        if ans.strip().lower() == 'y':
            return True
        return False


if __name__ == '__main__':
    loader = Loader()
    loader.remove_arrays()
    loader.load_qc()
    loader.load_icd_map()
    loader.load_icd()
    loader.load_qt()
    loader.load_icd_info()
