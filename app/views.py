import os
import re
import json
from flask import Flask, render_template, request, redirect, send_from_directory
from celery import Celery
from werkzeug.utils import secure_filename
import subprocess
from app import app
from zippy import zippyBatchQuery, zippyPrimerQuery
from zippylib import ascii_encode_dict
from zippylib.database import PrimerDB
import hashlib

ALLOWED_EXTENSIONS = set(['txt', 'batch', 'vcf', 'bed', 'csv'])
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['DOWNLOAD_FOLDER'] = 'results'
app.config['CONFIG_FILE'] = 'app/zippy.json'

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
    return send_from_directory(os.getcwd(), filename, as_attachment=True)

@app.route('/adhoc_result')
def adhoc_result(primerTable, resultList, missedIntervals):
    return render_template('file_uploaded.html', primerTable, resultList, missedIntervals)

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
        print "UPLOAD FILENAME", uploadFile.filename, filename

        # save file
        uploadedFile = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        uploadFile.save(uploadedFile)
        print "file saved to %s" % uploadedFile

        # open config file and database
        with open(app.config['CONFIG_FILE']) as conf:
            config = json.load(conf, object_hook=ascii_encode_dict)
            db = PrimerDB(config['database'])

        # create output folder
        downloadFolder = os.path.join(app.config['DOWNLOAD_FOLDER'], hashlib.sha1(open(uploadedFile).read()).hexdigest())
        subprocess.check_call(['mkdir', '-p', downloadFolder], shell=False)

        # run zippy
        print 'About to run Zippy'
        downloadFile = os.path.join(downloadFolder, filename)
        arrayOfFiles = zippyBatchQuery(config, uploadedFile, True, downloadFile, db)
        return render_template('file_uploaded.html', outputFiles=arrayOfFiles)
    else:
        print("file for upload not supplied or file-type not allowed")
        return redirect('/no_file')

@app.route('/adhoc_design/', methods=['POST'])
def adhocdesign():
    locus = request.form.get('locus')
    print locus
    store = request.form.get('store')
    print store
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
    # zippyPrimerQuery(config, targets, design=True, outfile=None, db=None, store=False)
    primerTable, resultList, missedIntervals = zippyPrimerQuery(config, locus, True, None, db, store)
    missedIntervalNames = []
    for interval in missedIntervals:
        missedIntervalNames.append(interval.name)
    # os.chdir('./app/')
    # print subprocess.call(['./zippy.py', 'get', locus, '--design', '--nostore'], shell=False)
    return render_template('/adhoc_result.html', primerTable=primerTable, resultList=resultList, missedIntervals=missedIntervalNames)


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
