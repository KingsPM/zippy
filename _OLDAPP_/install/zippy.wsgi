activate_this = '/usr/local/zippy/venv/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))
import sys
import logging
logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, '~')
sys.path.append('/usr/local/zippy')
from zippy import app as application
