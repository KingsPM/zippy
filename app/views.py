import os
import json
from flask import Flask, render_template, request, redirect, send_from_directory
from celery import Celery
from werkzeug.utils import secure_filename
import subprocess 
from app import app
from zippy import zippyBatchQuery
from zippylib import ascii_encode_dict
from zippylib.database import PrimerDB
import hashlib

ALLOWED_EXTENSIONS = set(['txt', 'batch', 'vcf', 'bed'])
UPLOAD_FOLDER = '/home/vagrant/dev/zippy/uploads'
DOWNLOAD_FOLDER = '/home/vagrant/dev/zippy/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['DOWNLOAD_FOLDER'] = DOWNLOAD_FOLDER

app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'

celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

@celery.task(bind=True)
def query_zippy(uploadFile):
    filename = secure_filename(uploadFile.filename)
    print filename
    uploadFile.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
    print "file saved to ./uploads/%s" % filename
    os.chdir('./app/')
    print subprocess.call(['./zippy.py', 'get', '--outfile', outfile, '../uploads/%s'% filename], shell=False)
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

@app.route('/file_uploaded/<path:filename>')
def download_file(filename):
    return send_from_directory(app.config['DOWNLOAD_FOLDER'],
                               filename, as_attachment=True)

# @app.route('/upload/', methods=['POST'])
# def upload():
#   uploadFile = request.files['filePath']
#   print("request of file successful")
#   if uploadFile and allowed_file(uploadFile.filename):
#       filename = secure_filename(uploadFile.filename)
#       print filename
#       print app.config['UPLOAD_FOLDER']
#       uploadFile.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
#       print "file saved to ./uploads/%s" % filename
#       os.chdir('./app/')
#       print subprocess.call(['./zippy.py', 'get', '--outfile', 'outfile', '../uploads/%s'% filename], shell=False)
#       print "running zippy..."
#       return redirect('/file_uploaded')
#   else:
#       print("file for upload not supplied or file-type not allowed")
#       return redirect('/no_file')



@app.route('/upload/', methods=['POST'])
def upload():
    uploadFile = request.files['filePath']
    print("request of file successful")
    if uploadFile and allowed_file(uploadFile.filename):
        filename = secure_filename(uploadFile.filename)
        print filename
        uploadFile.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        fileHash = hashlib.sha1(open('./uploads/%s' % filename).read()).hexidigest()
        subprocess.check.call(['mkdir', '-p', './results/%s' % fileHash])
        print "file saved to ./uploads/%s" % filename
        configFile = './app/zippy.json'
        targetFile = './uploads/%s' % filename
        print configFile
        print targetFile
        with open(configFile) as conf:
            config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
        print 'About to run Zippy'
        #zippyBatchQuery(config, targets, design=True, outfile=None, db=None)
        arrayOfFiles = zippyBatchQuery(config, './uploads/%s' % filename, True, filename, db)
        for eachFile in arrayOfFiles:
            subprocess.check.call(['mv', eachFile, './results/%s/' % fileHash])
        # print arrayOfFiles
        

        # print subprocess.call(['/home/vagrant/dev/zippy/app/zippy.py', '-c', '/home/vagrant/dev/zippy/app/zippy.json', 'batch', '--outfile', 'outfile', '/home/vagrant/dev/zippy/uploads/%s'% filename], shell=False)
        # os.chdir('./app/')
        # print subprocess.call(['./zippy.py', 'batch', '--outfile', 'outfile', '../uploads/%s'% filename], shell=False)
        # return redirect('/file_uploaded')
        return render_template('file_uploaded.html', outputFiles=arrayOfFiles)
    else:
        print("file for upload not supplied or file-type not allowed")
        return redirect('/no_file')


# @app.route('/upload/', methods=['POST'])
# def upload():
#     ourFile = request.files['filePath']
#     print("request of file successful")
#     if ourFile and allowed_file(ourFile.filename):

#         query = query_zippy.apply_async(args=[ourFile])
#         print type(query), dir(query), query.id
#         print("has run through zippy")
#         return redirect('/file_uploaded')
#     else:
#         print("file for upload not supplied or file-type not allowed")
#         return redirect('/no_file')





