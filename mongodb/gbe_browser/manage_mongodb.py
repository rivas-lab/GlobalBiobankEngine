from flask_script import Manager
from gbe_mongodb import app
import gbe_mongodb

manager = Manager(app)


@manager.command
def hello():
    print "hello"


@manager.command
def load_db():
    gbe_mongodb.load_db()
    
@manager.command
def load_icd_info_stats():
    gbe_mongodb.load_icd_info_stats()

@manager.command
def load_base_coverage():
    gbe_mongodb.load_base_coverage()

@manager.command
def load_base_icdstats():
    gbe_mongodb.load_base_icdstats()

@manager.command
def load_icd_stats():
    gbe_mongodb.load_icd_stats()

@manager.command
def load_variants_file():
    gbe_mongodb.load_variants_file()


@manager.command
def load_gene_models():
    gbe_mongodb.load_gene_models()


@manager.command
def load_dbsnp_file():
    gbe_mongodb.load_dbsnp_file()


@manager.command
def create_cache():
    gbe_mongodb.create_cache()


@manager.command
def precalculate_metrics():
    gbe_mongodb.precalculate_metrics()

if __name__ == "__main__":
    manager.run()

