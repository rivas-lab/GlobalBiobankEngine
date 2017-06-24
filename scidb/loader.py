import config
import logging
import numpy
import scidbpy


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

db = scidbpy.connect()


def load_qc_filters():
    qc_filters = []
    for rec in config.QC_FILTER_FILES:
        with open(rec['filename']) as fp:
            lines = fp.readlines()
            skip = rec.get('header', 0)
            qc_filters.extend(
                l.strip().split()[0] for l in lines[:-skip if skip else None])
            logger.info('QC:file:%s lines:%d skip:%d',
                        rec['filename'], len(lines), skip)

    # TODO fix when upload_data=list() implemented
    db.input(upload_data=numpy.array(qc_filters)).store('qc_filter')
    logger.info('Array:qc_filter records:%d', len(qc_filters))


if __name__ == '__main__':
    load_qc_filters()
