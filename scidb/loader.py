import config
import glob
import logging
import numpy
import os
import scidbpy
import subprocess
import tempfile

from future.moves.itertools import zip_longest


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

db = scidbpy.connect()


def load_qc():
    qc = []
    for rec in config.QC_FILES:
        with open(rec['filename']) as fp:
            lines = fp.readlines()
            skip = rec.get('header', 0)
            qc.extend(
                l.strip().split()[0] for l in lines[:-skip if skip else None])
            logger.info('QC:file:%s lines:%d skip:%d',
                        rec['filename'], len(lines), skip)

    # TODO fix when upload_data=list() implemented
    db.input(upload_data=numpy.array(qc)).store(config.QC_ARRAY)
    logger.info('Array:%s records:%d', config.QC_ARRAY, len(qc))


def load_icd():
    icd_lst = list(set(get_icd(fn) for fn in glob.iglob(config.ICD_GLOB)))
    icd_id_map = dict(zip(icd_lst, range(len(icd_lst))))
    db.input(config.ICD_INDEX_SCHEMA,
             upload_data=numpy.array(icd_lst)).store(config.ICD_INDEX_ARRAY)

    db.create_array(config.ICD_ARRAY, config.ICD_SCHEMA)
    try:
        fifo_names = [make_fifo() for i in range(config.SCIDB_INSTANCE_NUM)]
        file_iter = [glob.iglob(config.ICD_GLOB)] * config.SCIDB_INSTANCE_NUM
        instances = [str(i) for i in range(config.SCIDB_INSTANCE_NUM)]

        for file_names in zip_longest(*file_iter):

            # Last chunk might have None, strip them out
            if None in file_names:
                file_names = [file_name
                              for file_name in file_names if file_name]

            # Map filenames to icd_id based on src_instance_id
            icd_id_cond = 'null'
            for (inst, file_name) in zip(range(len(file_names)), file_names):
                icd_id = icd_id_map[get_icd(file_name)]
                icd_id_cond = ('iif(src_instance_id = {inst}, ' +
                               '{icd_id}, {prev})').format(
                                   inst=inst, icd_id=icd_id, prev=icd_id_cond)

            # Make pipes
            logger.info('Pipes:starting %d...', len(file_names))
            pipes = [make_pipe(file_name, fifo_name)
                     for (file_name, fifo_name) in zip(file_names, fifo_names)]

            query = config.ICD_QUERY.format(
                icd=config.ICD_ARRAY,
                qc=config.QC_ARRAY,
                paths=';'.join(fifo_names[:len(file_names)]),
                instances=';'.join(instances[:len(file_names)]),
                icd_id_cond=icd_id_cond)

            logger.info('Query:starting...')
            db.iquery(query)
            remove_versions(config.ICD_ARRAY)
            logger.info('Query:done')

            logger.info('Pipes:return code:%s',
                        ','.join(str(pipe.poll()) for pipe in pipes))
    finally:
        for fifo_name in fifo_names:
            remove_fifo(fifo_name)


def get_icd(filename):
    return os.path.basename(filename).split('.')[1]


def make_pipe(file_name, fifo_name):
    cmd = 'zcat {} > {}'.format(file_name, fifo_name)
    proc = subprocess.Popen(cmd, shell=True)

    logger.info('Spawn:%s pid:%d', cmd, proc.pid)
    return proc


def make_fifo():
    fifo_dir = tempfile.mkdtemp()
    fifo_name = os.path.join(fifo_dir, 'fifo')
    os.mkfifo(fifo_name)

    logger.info('FIFO:%s', fifo_name)
    return fifo_name


def remove_arrays():
    for array in (config.QC_ARRAY,
                  config.ICD_INDEX_ARRAY,
                  config.ICD_ARRAY):
        try:
            db.remove(array)
        except:
            pass


def remove_fifo(fifo_name):
    os.remove(fifo_name)
    os.rmdir(os.path.dirname(fifo_name))


def remove_versions(name):
    db.remove_versions(name,
                       db.versions(name)[:]['version_id'].max())


if __name__ == '__main__':
    remove_arrays()
    load_qc()
    load_icd()
