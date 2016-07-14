#!/usr/local/env python

import sys
import os
import re
import json
import hashlib
import subprocess
from flask import Flask, render_template, request, redirect, send_from_directory
from celery import Celery
from werkzeug.utils import secure_filename
from . import app
from .zippy import zippyBatchQuery, zippyPrimerQuery, updateLocation, searchByName
from .zippylib import ascii_encode_dict
from .zippylib.primer import Location
from .zippylib.database import PrimerDB

ALLOWED_EXTENSIONS = set(['txt', 'batch', 'vcf', 'bed', 'csv'])

app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['DOWNLOAD_FOLDER'] = 'results'
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
app.config['CONFIG_FILE'] = os.path.join(APP_ROOT, 'zippy.json')

# app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
# app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'
#
# celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
# celery.conf.update(app.config)
#
# @celery.task(bind=True)
# def query_zippy(uploadFile):
#     filename = secure_filename(uploadFile.filename)
#     print filename
#     uploadFile.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
#     print "file saved to ./uploads/%s" % filename
#     os.chdir('./app/')
#     print subprocess.call(['./zippy.py', 'get', '--outfile', outfile, '../uploads/%s'% filename], shell=False)
#     return "running zippy..."


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
    # read form
    uploadFile = request.files['filePath']
    deep = request.form.get('deep')
    predesign = request.form.get('predesign')
    design = request.form.get('design')
    outfile = request.form.get('outfile')
    if uploadFile and allowed_file(uploadFile.filename):
        filename = secure_filename(uploadFile.filename)
        print >> sys.stderr, "Uploaded: ", filename

        # save file
        uploadedFile = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        uploadFile.save(uploadedFile)
        print >> sys.stderr, "file saved to %s" % uploadedFile

        # open config file and database
        with open(app.config['CONFIG_FILE']) as conf:
            config = json.load(conf, object_hook=ascii_encode_dict)
            db = PrimerDB(config['database'])

        # create output folder
        downloadFolder = os.path.join(app.config['DOWNLOAD_FOLDER'], hashlib.sha1(open(uploadedFile).read()).hexdigest())
        subprocess.check_call(['mkdir', '-p', downloadFolder], shell=False)

        # run Zippy to design primers
        shortName = os.path.splitext(filename)[0]
        downloadFile = os.path.join(downloadFolder, outfile) if outfile else os.path.join(downloadFolder, shortName)
        arrayOfFiles, missedIntervalNames = zippyBatchQuery(config, uploadedFile, design, downloadFile, db, predesign, deep)
        return render_template('file_uploaded.html', outputFiles=arrayOfFiles, missedIntervals=missedIntervalNames)
    else:
        print("file for upload not supplied or file-type not allowed")
        return redirect('/no_file')

@app.route('/adhoc_design/', methods=['POST'])
def adhocdesign():
    # read form data
    uploadFile = request.files['filePath']
    locus = request.form.get('locus')
    design = request.form.get('design')
    deep = request.form.get('deep')
    store = request.form.get('store')
    # if locus:
    if re.match('\w{1,2}:\d+-\d+',locus) or (uploadFile and allowed_file(uploadFile.filename)):
        # get target
        if uploadFile:
            filename = secure_filename(uploadFile.filename)
            print >> sys.stderr, "Uploaded: ", filename
            target = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            uploadFile.save(target)
            print >> sys.stderr, "file saved to %s" % target
        else:
            target = locus
        # read config
        with open(app.config['CONFIG_FILE']) as conf:
            config = json.load(conf, object_hook=ascii_encode_dict)
            db = PrimerDB(config['database'])
        # run Zippy
        primerTable, resultList, missedIntervals = zippyPrimerQuery(config, target, design, None, db, store, deep)
        # get missed and render template
        missedIntervalNames = []
        for interval in missedIntervals:
            missedIntervalNames.append(interval.name)
        return render_template('/adhoc_result.html', primerTable=primerTable, resultList=resultList, missedIntervals=missedIntervalNames)
    else:
        print "no locus or file given"
        return render_template('/adhoc_result.html', primerTable=[], resultList=[], missedIntervals=[])

@app.route('/update_location/', methods=['POST'])
def updatePrimerLocation():
    primername = request.form.get('primername')
    vessel = request.form.get('vessel')
    well = request.form.get('well')
    force = request.form.get('force')
    try:
        assert primername
        loc = Location(vessel, well)
    except:
        print 'Please fill in all fields (PrimerName VesselNumber Well)'
        return render_template('location_updated.html', status=None)
    # read config
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
    # run zippy and render
    updateStatus = updateLocation(primername, loc, db, force)
    return render_template('location_updated.html', status=updateStatus)

# @app.route('/update_location_from_table/', methods=['POST'])
@app.route('/update_location_from_table/<primerName>', methods=['POST'])
def updateLocationFromTable(primerName):
    combinedlocation = request.form.get('combinedlocation')
    force = request.form.get('force')
    splitloc = combinedlocation.split('-')
    vessel = splitloc[0]
    well = splitloc[1]
    try:
        assert primerName
        loc = Location(vessel, well)
    except:
        print 'Please fill in all fields (PrimerName VesselNumber Well)'
        return render_template('location_updated.html', status=None)
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
    print primerName, loc, db, force
    # run zippy and render
    updateStatus = updateLocation(primerName, loc, db, force)
    return render_template('location_updated.html', status=updateStatus)


@app.route('/search_by_name/', methods=['POST'])
def searchName():
    searchName = request.form.get('searchName')
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
        searchResult = searchByName(searchName, db)
    return render_template('searchname_result.html', searchResult=searchResult, searchName=searchName)
