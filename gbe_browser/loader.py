import StringIO
import argparse
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

    def store(self, file_name, query, array_name):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(file_name, fifo_name)

        logger.info('Query:running...')
        self.db.iquery(query.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', array_name)

    # -- -
    # -- - QC - --
    # -- -
    def store_qc(self):
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

        # TODO fix when upload_data=list() implemented in SciDB-Py
        self.db.input(upload_data=numpy.array(qc)).store(config.QC_ARRAY)
        logger.info('Array:%s records:%d', config.QC_ARRAY, len(qc))

    # -- -
    # -- - ICD - --
    # -- -
    def set_icd_qt_lists(self):
        if config.ICD_INFO_ARRAY in dir(self.db.arrays):
            self.is_first_time = False

            logger.info('Query:running...')
            icd_info = self.db.iquery(config.ICD_INFO_ICD_QUERY,
                                      fetch=True,
                                      atts_only=True,
                                      schema=config.ICD_INFO_ICD_SCHEMA)
            logger.info('Found:%d records in %s',
                        len(icd_info), config.ICD_INFO_ARRAY)

            print(icd_info)
            icd_orphan = [str(i['icd_idx']['val'])
                          for i in icd_info
                          if i['icd']['null'] != 255]
            logger.info('Found:%d orphan ICDs, consider removing them',
                        len(icd_orphan))
            logger.info('Found:orphan:icd_idx:%s', ','.join(icd_orphan))

            self.icd_idx_map = dict(
                (i['icd']['val'], i['icd_idx']['val'])
                for i in icd_info
                if i['icd']['null'] == 255)
            self.icd_exist = set(self.icd_idx_map.keys())
            icd_chrom = set(i['icd']['val']
                            for i in icd_info
                            if (i['icd']['null'] == 255 and
                                (i['num_chrom']['null'] != 255 or
                                 i['num_chrom']['val'] < config.CHROM_MAX)))
        else:
            self.is_first_time = True

            self.db.create_array(config.ICD_INFO_ARRAY, config.ICD_INFO_SCHEMA)
            self.db.create_array(config.ICD_ARRAY, config.ICD_SCHEMA)

            self.icd_idx_map = {}
            self.icd_exist = set()

        logger.info('Found:%d files in %s',
                    sum(1 for _ in glob.iglob(
                        os.path.join(config.ICD_PATH, '*'))),
                    config.ICD_PATH)

        icd_files = [(Loader.get_icd(*Loader.get_icd_parts(fn)), fn)
                      for fn in glob.iglob(config.ICD_GLOB)]
        logger.info('Found:%d ICD files matching %s',
                    len(icd_files),
                    config.ICD_GLOB)
        self.icd_lst = [(icd, fn)
                        for (icd, fn) in icd_files
                        if self.is_first_time or icd not in self.icd_exist]
        logger.info('Select:%d ICD files of new ICDs', len(self.icd_lst))

        if not self.is_first_time:
            icd_lst_chrom = [(icd, fn)
                             for (icd, fn) in icd_files
                             if icd in icd_chrom]
            logger.info('Select:%d ICD files ' +
                        'for ICDs with less than %s chromosomes',
                        len(icd_lst_chrom),
                        config.CHROM_MAX)
            self.icd_lst += icd_lst_chrom
            logger.info('ICDs to be loaded: {}'.format(
                ','.join(set(icd for (icd, fn) in self.icd_lst))))

        qt_files = [(Loader.get_icd(*Loader.get_icd_parts(fn)), fn)
                    for fn in glob.iglob(config.QT_GLOB)]
        logger.info('Found:%d QT files matching %s',
                    len(qt_files),
                    config.QT_GLOB)

        self.qt_lst = [(icd, fn)
                       for (icd, fn) in qt_files
                       if self.is_first_time or icd not in self.icd_exist]
        logger.info('Select:%d QT files of new ICDs', len(self.qt_lst))

        if not self.is_first_time:
            qt_lst_chrom = [(icd, fn)
                             for (icd, fn) in qt_files
                             if icd in icd_chrom]
            logger.info('Select:%d QT files ' +
                        'for ICDs with less than %s chromosomes',
                        len(qt_lst_chrom),
                        config.CHROM_MAX)
            self.qt_lst += qt_lst_chrom
            logger.info('QTs to be loaded: {}'.format(
                ','.join(set(icd for (icd, fn) in self.qt_lst))))

        if not self.icd_lst and not self.qt_lst:
            raise Exception(('No {}ICD or QT files found. ' +
                             'Patterns used:\n{}\n{}').format(
                                '' if self.is_first_time else 'new ',
                                config.ICD_GLOB,
                                config.QT_GLOB))

    def append_icd_info(self):
        icd_qt_lst = set(
            icd
            for (icd, fn) in itertools.chain(self.icd_lst, self.qt_lst)
            if icd not in self.icd_exist)

        if self.is_first_time:
            icd_st = 0
        else:
            icd_st = self.db.iquery(config.ICD_INFO_IDX_MAX_QUERY,
                                    schema=config.ICD_INFO_IDX_MAX_SCHEMA,
                                    fetch=True,
                                    atts_only=True)[0]['max']['val'] + 1

        self.icd_idx_map.update(
            dict(zip(icd_qt_lst,
                     range(icd_st, icd_st + len(icd_qt_lst)))))
        if icd_qt_lst:
            self.db.iquery(
                config.ICD_INFO_APPEND_QUERY.format(
                    fn='{fn}',
                    start=icd_st,
                    stop=icd_st + len(icd_qt_lst)),
                upload_data=StringIO.StringIO('\n'.join(icd_qt_lst)))
        logger.info('Array:%s', config.ICD_INFO_ARRAY)

    def insert_icd_info(self):
        self.db.iquery(config.ICD_INFO_INSERT_QUERY)
        logger.info('Array:%s', config.ICD_INFO_ARRAY)

    def insert_icd(self):
        self.icd_lst.sort(key=lambda k: k[1])
        file_iter = ([(fn for (icd, fn) in self.icd_lst)] *
                     config.SCIDB_INSTANCE_NUM)

        for file_names in zip_longest(*file_iter):

            # Last chunk might have None, strip them out
            if None in file_names:
                file_names = [file_name
                              for file_name in file_names if file_name]

            # Get conditional expressions
            (icd_idx_cond, icdind_cond) = self.get_icd_cond(file_names)

            # Make pipes
            logger.info('Pipes:starting %d...', len(file_names))
            pipes = [Loader.make_pipe(file_name, fifo_name)
                     for (file_name, fifo_name) in zip(file_names,
                                                       self.fifo_names)]

            query = config.ICD_INSERT_QUERY.format(
                paths=';'.join(self.fifo_names[:len(file_names)]),
                instances=';'.join(self.instances[:len(file_names)]),
                icd_idx_cond=icd_idx_cond,
                icdind_cond=icdind_cond)

            logger.info('Query:running...')
            self.db.iquery(query)
            self.remove_versions(config.ICD_ARRAY)
            logger.info('Query:done')

            logger.info('Pipes:return code:%s',
                        ','.join(str(pipe.poll()) for pipe in pipes))
        logger.info('Array:%s', config.ICD_ARRAY)

    def insert_qt(self):
        self.qt_lst.sort(key=lambda k: k[1])
        file_iter = ([(fn for (icd, fn) in self.qt_lst)] *
                     config.SCIDB_INSTANCE_NUM)

        for file_names in zip_longest(*file_iter):

            # Last chunk might have None, strip them out
            if None in file_names:
                file_names = [file_name
                              for file_name in file_names if file_name]

            # Get conditional expressions
            (icd_idx_cond, icdind_cond) = self.get_icd_cond(file_names)

            # Make pipes
            logger.info('Pipes:starting %d...', len(file_names))
            pipes = [Loader.make_pipe(file_name, fifo_name)
                     for (file_name, fifo_name) in zip(file_names,
                                                       self.fifo_names)]

            query = config.QT_INSERT_QUERY.format(
                paths=';'.join(self.fifo_names[:len(file_names)]),
                instances=';'.join(self.instances[:len(file_names)]),
                icd_idx_cond=icd_idx_cond,
                icdind_cond=icdind_cond)

            logger.info('Query:running...')
            self.db.iquery(query)
            self.remove_versions(config.ICD_ARRAY)
            logger.info('Query:done')

            logger.info('Pipes:return code:%s',
                        ','.join(str(pipe.poll()) for pipe in pipes))
        logger.info('Array:%s', config.ICD_ARRAY)

    def get_icd_cond(self, file_names):
        """Build SciDB conditional expressions ("iif") to map file names to
        "icd_idx", "icdind" using "src_instance_id"

        """
        icd_idx_cond = 'null'
        icdind_cond = 'null'
        for (inst, file_name) in zip(range(len(file_names)), file_names):
            (prefix, suffix, intadd, ind) = Loader.get_icd_parts(file_name)

            icd_idx = self.icd_idx_map[
                Loader.get_icd(prefix, suffix, intadd, ind)]
            icd_idx_cond = ('iif(src_instance_id = {inst}, ' +
                            '{icd_idx}, {prev})').format(
                                inst=inst, icd_idx=icd_idx, prev=icd_idx_cond)

            icdind_cond = ('iif(src_instance_id = {inst}, ' +
                           '\'{val}\', {prev})').format(
                               inst=inst,
                               val=str(intadd) + ind,
                               prev=icdind_cond)
        return (icd_idx_cond, icdind_cond)

    def store_affyid_index(self):
        logger.info('Query:running...')
        self.db.iquery(config.AFFYID_INDEX_STORE_QUERY)
        logger.info('Array:%s', config.AFFYID_INDEX_ARRAY)

    def store_icd_affyid(self):
        logger.info('Query:running...')
        self.db.iquery(config.ICD_AFFYID_STORE_QUERY)
        logger.info('Array:%s', config.ICD_AFFYID_ARRAY)

    # -- -
    # -- - GENE - --
    # -- -
    def store_gene_index(self):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(
            config.GENE_FILE,
            fifo_name,
            "zcat {file_name} | grep '\tgene\t' | " +
            "cut --fields=9 | cut --delimiter=' ' --fields=2 | "
            "sort | uniq > {fifo_name}")

        logger.info('Query:running...')
        self.db.iquery(
            config.GENE_INDEX_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', config.GENE_INDEX_ARRAY)

    def store_transcript_index(self):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(
            config.GENE_FILE,
            fifo_name,
            "zcat {file_name} | grep '\ttranscript\t' | " +
            "cut --fields=9 | cut --delimiter=' ' --fields=4 | " +
            "sort | uniq > {fifo_name}")

        logger.info('Query:running...')
        self.db.iquery(
            config.TRANSCRIPT_INDEX_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', config.TRANSCRIPT_INDEX_ARRAY)

    def store_dbnsfp(self):
        self.store(config.DBNSFP_FILE,
                   config.DBNSFP_STORE_QUERY,
                   config.DBNSFP_ARRAY)

    def store_canonical(self):
        self.store(config.CANONICAL_FILE,
                   config.CANONICAL_STORE_QUERY,
                   config.CANONICAL_ARRAY)

    def store_omim(self):
        self.store(config.OMIM_FILE,
                   config.OMIM_STORE_QUERY,
                   config.OMIM_ARRAY)

    def store_gene(self):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(
            config.GENE_FILE,
            fifo_name,
            "zcat {file_name} | grep '\tgene\t' > {fifo_name}")

        logger.info('Query:running...')
        self.db.iquery(config.GENE_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', config.GENE_ARRAY)

    def store_transcript(self):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(
            config.GENE_FILE,
            fifo_name,
            "zcat {file_name} | grep '\ttranscript\t' > {fifo_name}")

        logger.info('Query:running...')
        self.db.iquery(config.TRANSCRIPT_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', config.TRANSCRIPT_ARRAY)

    def store_exon(self):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(
            config.GENE_FILE,
            fifo_name,
            "zcat {file_name} | " +
            "grep --perl-regexp '\\t(exon|CDS|UTR)\\t' > {fifo_name}")

        logger.info('Query:running...')
        self.db.iquery(config.EXON_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', config.EXON_ARRAY)

    def store_gene_info(self):
        logger.info('Query:running...')
        self.db.iquery(config.GENE_C_TRANSCRIPT_INFO_STORE_QUERY)
        logger.info('Query:done')
        logger.info('Query:running...')
        self.db.iquery(config.GENE_TRANSCRIPT_INFO_STORE_QUERY)
        logger.info('Query:done')
        logger.info('Query:running...')
        self.db.iquery(config.GENE_EXON_INFO_STORE_QUERY)
        logger.info('Query:done')
        logger.info('Array:%s', config.GENE_ARRAY)

    def store_transcript_info(self):
        logger.info('Query:running...')
        self.db.iquery(config.TRANSCRIPT_EXON_INFO_STORE_QUERY)
        logger.info('Query:done')
        logger.info('Array:%s', config.TRANSCRIPT_ARRAY)

    # -- -
    # -- - VARIANT - --
    # -- -
    def store_variant(self):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(config.VARIANT_FILE, fifo_name)

        logger.info('Query:running...')
        self.db.iquery(config.VARIANT_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', config.VARIANT_ARRAY)

    def store_variant_index(self, key, query, array):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(
            config.VARIANT_FILE,
            fifo_name,
            'zcat {{file_name}} | python {cmd} {key} > {{fifo_name}}'.format(
                cmd=os.path.join(
                      os.path.dirname(os.path.realpath(__file__)),
                      'unnest.py'),
                key=key))

        logger.info('Query:running...')
        self.db.iquery(query.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', array)

    def store_variant_gene(self):
        self.store_variant_index('Gene',
                                 config.VARIANT_GENE_STORE_QUERY,
                                 config.VARIANT_GENE_ARRAY)

    def store_variant_transcript(self):
        self.store_variant_index('Feature',
                                 config.VARIANT_TRANSCRIPT_STORE_QUERY,
                                 config.VARIANT_TRANSCRIPT_ARRAY)

    # -- -
    # -- - COVERAGE - --
    # -- -
    def store_coverage(self):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(config.COVERAGE_FILE, fifo_name)

        logger.info('Query:running...')
        self.db.iquery(config.COVERAGE_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', config.COVERAGE_ARRAY)

    # -- -
    # -- - DBSNP - --
    # -- -
    def store_dbsnp(self):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(config.DBSNP_FILE, fifo_name)

        logger.info('Query:running...')
        self.db.iquery(
            config.DBSNP_BY_RSID_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', config.DBSNP_BY_RSID_ARRAY)

        logger.info('Query:running...')
        self.db.iquery(
            config.DBSNP_BY_CHROM_POS_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Array:%s', config.DBSNP_BY_CHROM_POS_ARRAY)

    # -- -
    # -- - BIM - --
    # -- -
    def store_bim(self):
        fifo_name = self.fifo_names[0]
        pipe = Loader.make_pipe(config.BIM_FILE, fifo_name)

        logger.info('Query:running...')
        self.db.iquery(config.BIM_STORE_QUERY.format(path=fifo_name))
        logger.info('Query:done')
        logger.info('Pipe:return code:%s', pipe.poll())
        logger.info('Array:%s', config.BIM_ARRAY)

    # -- - - --
    def remove_arrays(self):
        if not Loader.confirm('Remove and recreate arrays'):
            return
        for array in (config.QC_ARRAY,
                      config.ICD_INFO_ARRAY,
                      config.ICD_ARRAY,
                      config.AFFYID_INDEX_ARRAY,
                      config.ICD_AFFYID_ARRAY,
                      config.GENE_INDEX_ARRAY,
                      config.DBNSFP_ARRAY,
                      config.CANONICAL_ARRAY,
                      config.OMIM_ARRAY,
                      config.GENE_INDEX_ARRAY,
                      config.TRANSCRIPT_INDEX_ARRAY,
                      config.GENE_ARRAY,
                      config.TRANSCRIPT_ARRAY,
                      config.EXON_ARRAY,
                      config.VARIANT_ARRAY,
                      config.VARIANT_GENE_ARRAY,
                      config.VARIANT_TRANSCRIPT_ARRAY,
                      config.COVERAGE_ARRAY,
                      config.DBSNP_BY_RSID_ARRAY,
                      config.DBSNP_BY_CHROM_POS_ARRAY,
                      config.BIM_ARRAY):
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
        elif 'INI' in parts:
            prefix = 'INI'
            intadd = 30
        elif 'initialdata' in name:
            prefix = 'INI'
            intadd = 30
        elif 'INI' in name:
            prefix = 'INI'
            intadd = 30
        elif 'qt' in name:
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
        elif 'MED' in name:
            prefix = 'MED'
            intadd = 80

        # Figure out suffix
        ind = parts[1].split(
            '_')[0].strip()[1:].strip(
                'RH').strip(
                    'FH').strip(
                        'qt').strip(
                        'INI').strip(
                        'initialdata').strip(
                            'cancer').strip(
                                'HC').split(
                                    '_FH2')[0]

        suffix = parts[1].strip(
            'RH').strip(
                'FH').strip(
                    'qt').strip(
                    'INI').strip(
                    'initialdata').strip(
                        'MED').strip(
                            'cancer').strip(
                                'HC').split(
                                    '_FH2')[0]

        return (prefix, suffix, intadd, ind)

    @staticmethod
    def get_icd(prefix, suffix, intadd, ind):
        return prefix + suffix

    @staticmethod
    def make_pipe(file_name, fifo_name, pipe='zcat {file_name} > {fifo_name}'):
        cmd = pipe.format(file_name=file_name, fifo_name=fifo_name)
        proc = subprocess.Popen(cmd, shell=True)

        logger.info('Spawn:%s pid:%d', cmd, proc.pid)
        return proc

    @staticmethod
    def make_fifo():
        fifo_dir = tempfile.mkdtemp()
        os.chmod(fifo_dir, 0755)
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
    parser = argparse.ArgumentParser(description='SciDB Loader')
    parser.add_argument('action',
                        choices=['all', 'icd-incremental'],
                        help='Loading choice')

    args = parser.parse_args()
    loader = Loader()

    if args.action == 'all':
        loader.remove_arrays()

        loader.store_qc()
        loader.set_icd_qt_lists()
        loader.append_icd_info()
        loader.insert_icd_info()
        loader.insert_icd()
        loader.insert_qt()
        loader.store_affyid_index()
        loader.store_icd_affyid()

        loader.store_gene_index()
        loader.store_transcript_index()
        loader.store_dbnsfp()
        loader.store_canonical()
        loader.store_omim()
        loader.store_gene()
        loader.store_transcript()
        loader.store_exon()
        loader.store_gene_info()
        loader.store_transcript_info()

        loader.store_variant()
        loader.store_variant_gene()
        loader.store_variant_transcript()

        loader.store_coverage()
        loader.store_dbsnp()
        loader.store_bim()

    elif args.action == 'icd-incremental':
        loader.set_icd_qt_lists()
        loader.append_icd_info()
        loader.insert_icd_info()
        loader.insert_icd()
        loader.insert_qt()
        loader.store_affyid_index()
        loader.store_icd_affyid()
