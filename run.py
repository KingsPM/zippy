#!flask/bin/python

from app import app

#app.run(debug=True)  # not accessible outside guest
app.run(host='0.0.0.0')
