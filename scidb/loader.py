import config
import glob
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


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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
            with open(rec['filename']) as fp:
                lines = fp.readlines()
                skip = rec.get('header', 0)
                qc.extend(
                    l.strip().split()[0]
                    for l in lines[:-skip if skip else None])
                logger.info('QC:file:%s lines:%d skip:%d',
                            rec['filename'], len(lines), skip)

        # TODO fix when upload_data=list() implemented
        self.db.input(upload_data=numpy.array(qc)).store(config.QC_ARRAY)
        logger.info('Array:%s records:%d', config.QC_ARRAY, len(qc))

    def load_icd_map(self):
        icd_lst = list(set(Loader.get_icd(fn)
                           for fn in glob.iglob(config.ICD_GLOB)))
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

            # Map filenames to icd_id based on src_instance_id
            icd_id_cond = 'null'
            for (inst, file_name) in zip(range(len(file_names)), file_names):
                icd_id = self.icd_id_map[Loader.get_icd(file_name)]
                icd_id_cond = ('iif(src_instance_id = {inst}, ' +
                               '{icd_id}, {prev})').format(
                                   inst=inst, icd_id=icd_id, prev=icd_id_cond)

            # Make pipes
            logger.info('Pipes:starting %d...', len(file_names))
            pipes = [Loader.make_pipe(file_name, fifo_name)
                     for (file_name, fifo_name) in zip(file_names,
                                                       self.fifo_names)]

            query = config.ICD_QUERY.format(
                icd=config.ICD_ARRAY,
                qc=config.QC_ARRAY,
                paths=';'.join(self.fifo_names[:len(file_names)]),
                instances=';'.join(self.instances[:len(file_names)]),
                icd_id_cond=icd_id_cond)

            logger.info('Query:starting...')
            self.db.iquery(query)
            self.remove_versions(config.ICD_ARRAY)
            logger.info('Query:done')

            logger.info('Pipes:return code:%s',
                        ','.join(str(pipe.poll()) for pipe in pipes))

    def remove_arrays(self):
        for array in (config.QC_ARRAY,
                      config.ICD_INDEX_ARRAY,
                      config.ICD_ARRAY):
            try:
                self.db.remove(array)
            except:
                pass

    def remove_versions(self, name):
        self.db.remove_versions(name,
                                self.db.versions(name)[:]['version_id'].max())

    @staticmethod
    def get_icd(filename):
        return os.path.basename(filename).split('.')[1]

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


if __name__ == '__main__':
    loader = Loader()
    loader.remove_arrays()
    loader.load_qc()
    loader.load_icd_map()
    loader.load_icd()
