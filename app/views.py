import os
import re
import json
from flask import Flask, render_template, request, redirect, send_from_directory
from celery import Celery
from werkzeug.utils import secure_filename
import subprocess
from app import app
from zippy import zippyBatchQuery, zippyPrimerQuery, updateLocation, searchByName
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

@app.route('/location_updated')
def location_updated(status):
    return render_template('location_updated.html', status)


@app.route('/upload/', methods=['POST'])
def upload():
    uploadFile = request.files['filePath']
    print("request of file successful")
    if uploadFile and allowed_file(uploadFile.filename):
        filename = secure_filename(uploadFile.filename)
        print "Uploaded: ", filename

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

        # run Zippy to design primers
        print 'About to run Zippy'
        shortName = os.path.splitext(filename)[0]
        downloadFile = os.path.join(downloadFolder, shortName)
        arrayOfFiles, missedIntervals = zippyBatchQuery(config, uploadedFile, True, downloadFile, db)
        missedIntervalNames = []
        for interval in missedIntervals:
            for i in interval:
                missedIntervalNames.append(i.name)
        print 'Missed intervals: ', missedIntervalNames
        return render_template('file_uploaded.html', outputFiles=arrayOfFiles, missedIntervals=missedIntervalNames)
    else:
        print("file for upload not supplied or file-type not allowed")
        return redirect('/no_file')

@app.route('/adhoc_design/', methods=['POST'])
def adhocdesign():
    locus = request.form.get('locus')
    # if locus:
    if re.match('\w{1,2}:\d+-\d+',locus):
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
    else:
        print "locus not given"
        return render_template('/adhoc_result.html', primerTable=[], resultList=[], missedIntervals=[])

@app.route('/update_location/', methods=['POST'])
def update_Location():
    location = request.form.get('location')
    if re.match('\w+\s\w+\s\w',location):
        with open(app.config['CONFIG_FILE']) as conf:
            config = json.load(conf, object_hook=ascii_encode_dict)
            db = PrimerDB(config['database'])
        updateStatus = updateLocation(location, db)
        print updateStatus
        return render_template('location_updated.html', status=updateStatus)
    else:
        print 'update location not given in correct format'
        return render_template('location_updated.html', status=None)

@app.route('/search_by_name/', methods=['POST'])
def searchName():
    searchName = request.form.get('searchName')
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
        searchResult = searchByName(searchName, db)
        for pairs in searchResult:
            for result in pairs:
                print result.name
    return render_template('searchname_result.html', searchResult=searchResult, searchName=searchName)

    



