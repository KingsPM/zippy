#!venv/bin/python

'''runs flask development server'''

from zippy import app

app.debug = True
app.run(host='0.0.0.0')
