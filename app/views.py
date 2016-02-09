import os
from flask import Flask, render_template, request, redirect
from celery import Celery
# from flask.ext.wtf import Form
from werkzeug.utils import secure_filename
import subprocess 
from app import app

ALLOWED_EXTENSIONS = set(['txt', 'batch', 'vcf', 'bed'])
UPLOAD_FOLDER = '/home/vagrant/dev/zippy/uploads'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'

celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

@celery.task
def query_zippy(file):
    # some long running task here
	# filename = secure_filename(file.filename)
	# print filename
	# file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
	# print "file saved to ./uploads/%s" % filename
	# os.chdir('./app/')
	# print subprocess.call(['./zippy.py', 'get', '-o', outfile, '../uploads/%s'% filename], shell=False)
	return "running zippy..."


def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

@app.route('/')
@app.route('/index')
def index():
	return render_template('index.html')

@app.route('/no_file')
def no_file():
	return render_template('no_file.html')

@app.route('/file_uploaded')
def file_uploaded():
	return render_template('file_uploaded.html')

# @app.route('/upload/', methods=['POST'])
# def upload():
# 	file = request.files['filePath']
# 	print("request of file successful")
# 	if file and allowed_file(file.filename):
# 		filename = secure_filename(file.filename)
# 		print filename
# 		file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
# 		print "file saved to ./uploads/%s" % filename
# 		os.chdir('./app/')
# 		print subprocess.call(['./zippy.py', 'get', '--outfile', 'outfile', '../uploads/%s'% filename], shell=False)
# 		# print os.system('./zippy.py batch ././uploads/%s'% filename)'./zippy.py get -o', outfile, '../uploads/%s'% filename
# 		# print list_content.communicate()
# 		print "running zippy..."
# 		return redirect('/file_uploaded')
# 	else:
# 		print("file for upload not supplied or file-type not allowed")
# 		return redirect('/no_file')

@app.route('/upload/', methods=['POST'])
def upload():
	ourFile = request.files['filePath']
	print("request of file successful")
	if ourFile and allowed_file(ourFile.filename):
		query_zippy.apply_async(args=[ourFile])
		return redirect('/file_uploaded')
	else:
		print("file for upload not supplied or file-type not allowed")
		return redirect('/no_file')





