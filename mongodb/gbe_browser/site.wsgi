#!//local/miniconda/bin/python

import sys
sys.path.insert(0, '/local/ibd2016meta/exac/exac_browser')
sys.stderr.write("site.wsgi: sys.path: %s \n" % ", ".join(sys.path))

sys.stderr.write("site.wsgi: from exac import app \n")
from exac import app  # as application
sys.stderr.write("site.wsgi: from werkzeug.debug import DebuggedApplication \n")
from werkzeug.debug import DebuggedApplication
sys.stderr.write("site.wsgi: application = DebuggedApplication(app, True) \n")
application = DebuggedApplication(app, True)

sys.stderr.write("site.wsgi: done \n")


import logging
logging.basicConfig(stream=sys.stderr)

