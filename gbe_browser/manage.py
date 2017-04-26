from flask.ext.script import Manager
from gbe import app
import gbe

manager = Manager(app)


@manager.command
def hello():
    print "hello"


@manager.command
def load_db():
    gbe.load_db()
    
@manager.command
def load_icd_info_stats():
    gbe.load_icd_info_stats()

@manager.command
def load_base_coverage():
    gbe.load_base_coverage()

@manager.command
def load_base_icdstats():
    gbe.load_base_icdstats()

@manager.command
def load_icd_stats():
    gbe.load_icd_stats()

@manager.command
def load_variants_file():
    gbe.load_variants_file()


@manager.command
def load_gene_models():
    gbe.load_gene_models()


@manager.command
def load_dbsnp_file():
    gbe.load_dbsnp_file()


@manager.command
def create_cache():
    gbe.create_cache()


@manager.command
def precalculate_metrics():
    gbe.precalculate_metrics()

if __name__ == "__main__":
    manager.run()

