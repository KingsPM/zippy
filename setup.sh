# this flask app is developed following the Flask tutorial at
# http://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-i-hello-world
# and the O'reilly book by Miguel Grinberg 'Flask Web Development'


sudo pip install virtualenv
virtualenv flask

flask/bin/pip install flask
flask/bin/pip install flask-login
flask/bin/pip install flask-openid
flask/bin/pip install flask-mail
flask/bin/pip install flask-sqlalchemy
flask/bin/pip install sqlalchemy-migrate
flask/bin/pip install flask-whooshalchemy
flask/bin/pip install flask-wtf
flask/bin/pip install flask-babel
flask/bin/pip install guess_language
flask/bin/pip install flipflop
flask/bin/pip install coverage

mkdir app
mkdir app/static
mkdir app/templates
mkdir tmp
