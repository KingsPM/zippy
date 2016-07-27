#!/usr/local/env python

import sys
import os
import re
import json
import hashlib
import subprocess
from flask import Flask, render_template, request, redirect, send_from_directory, session, flash
from celery import Celery
from werkzeug.utils import secure_filename
from . import app
from .zippy import zippyBatchQuery, zippyPrimerQuery, updateLocation, searchByName, updatePrimerName, updatePrimerPairName, blacklistPair, readprimerlocations
from .zippylib import ascii_encode_dict
from .zippylib.primer import Location
from .zippylib.database import PrimerDB

ALLOWED_EXTENSIONS = set(['txt', 'batch', 'vcf', 'bed', 'csv'])

app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['DOWNLOAD_FOLDER'] = 'results'
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
app.config['CONFIG_FILE'] = os.path.join(APP_ROOT, 'zippy.json')
app.secret_key = 'someKey'

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


@app.route('/upload/', methods=['POST', 'GET'])
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

@app.route('/select_pair_to_update/<pairName>')
def pair_to_update(pairName):
    return render_template('update_pair.html', pairName=pairName)

@app.route('/update_pair_name/<pairName>', methods=['POST'])
def update_pair_name(pairName):
    newName = request.form.get('name')
    print pairName
    print newName
    if newName == pairName:
        flash('Pair renaming failed - new name is the same as current', 'warning')
        return render_template('update_pair.html', pairName=pairName)
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
        nameChange = updatePrimerPairName(pairName, newName, db)
        print nameChange
        if nameChange:
            flash('Pair "%s" renamed "%s"' % (pairName, newName), 'success')
        else:
            flash('Pair renaming failed', 'warning')
    return render_template('update_pair.html', pairName=newName)

@app.route('/select_primer_to_update/<primerName>/<primerLoc>')
def primer_to_update(primerName, primerLoc):
    print primerName
    print primerLoc
    primerInfo = primerName + '|' + primerLoc
    print primerInfo
    return redirect('/update_location_from_table/%s' % (primerInfo))

@app.route('/update_location_from_table/<primerInfo>', methods=['GET','POST'])
def updateLocationFromTable(primerInfo):
    if request.method == 'POST':
        primerName = primerInfo
        vessel = request.form.get('vessel')
        well = request.form.get('well')
        force = request.form.get('force')
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
        if updateStatus[0] == 'occupied':
            flash('Location already occupied by %s' % (' and '.join(updateStatus[1])), 'warning')
        elif updateStatus[0] == 'success':
            flash('%s location sucessfully set to %s' % (primerName, str(loc)), 'success')
        else:
            flash('%s location update to %s failed' % (primerName, str(loc)), 'warning')
        return render_template('update_location_from_table.html', primerName=primerName, primerLoc=vessel + '-' + well)
    else:
        print primerInfo
        splitInfo = primerInfo.split('|')
        primerName = splitInfo[0]
        primerLoc = splitInfo[1]
        print primerName
        print primerLoc
        return render_template('update_location_from_table.html', primerName=primerName, primerLoc=primerLoc)

@app.route('/select_primer_to_rename/<primerName>/<primerLoc>', methods=['POST'])
def primer_to_rename(primerName, primerLoc):
    newName = request.form.get('name')
    print primerName
    print primerLoc
    print newName
    primerInfo = primerName + '|' + primerLoc + '|' + newName
    print primerInfo
    return redirect('/update_primer_name/%s' % (primerInfo))

@app.route('/update_primer_name/<primerInfo>')
def update_name_of_primer(primerInfo):
    splitInfo = primerInfo.split('|')
    currentName = splitInfo[0]
    primerLoc = splitInfo[1]
    newName = splitInfo[2]
    if newName == currentName:
        flash('Primer renaming failed - new name is the same as current', 'warning')
        return render_template('update_location_from_table.html', primerName=newName, primerLoc=primerLoc)
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
        nameChange = updatePrimerName(currentName, newName, db)
        print nameChange
        if nameChange:
            flash('Primer "%s" renamed "%s"' % (currentName, newName), 'success')
        else:
            flash('Primer renaming failed', 'warning')
    return render_template('update_location_from_table.html', primerName=newName, primerLoc=primerLoc)

@app.route('/specify_searchname/', methods=['POST'])
def searchName():
    searchName = request.form.get('searchName')
    session['searchName'] = searchName
    return redirect('/search_by_name/')

@app.route('/search_by_name/')
def search_by_name():
    searchName = session['searchName']
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
        searchResult = searchByName(searchName, db)
    return render_template('searchname_result.html', searchResult=searchResult, searchName=searchName)

@app.route('/blacklist_pair/<pairname>', methods=['POST'])
def blacklist_pair(pairname):
    print 'This is the pairname: ' + pairname
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
        blacklisted = blacklistPair(pairname, db)
        for b in blacklisted:
            flash('%s added to blacklist' % (b,), 'success')
    return redirect('/search_by_name/')

@app.route('/upload_batch_locations/', methods=['POST'])
def upload_samplesheet():
    if request.method == 'POST':
        locationsheet = request.files['locationsheet']
        if not locationsheet:
            flash('No csvfile submitted. Please try again','warning')
        else:
            saveloc = 'uploads/'+locationsheet.filename
            locationsheet.save(saveloc)
            updateList = readprimerlocations(saveloc)
            for item in updateList:
                print 'Primer:locations to update: ', item
                with open(app.config['CONFIG_FILE']) as conf:
                    config = json.load(conf, object_hook=ascii_encode_dict)
                    db = PrimerDB(config['database'])
                    updateStatus = updateLocation(item[0], item[1], db, True) # Force is set to True, will force primers into any occupied locations
                    if updateStatus[0] == 'occupied':
                        flash('Location already occupied by %s' % (' and '.join(updateStatus[1])), 'warning')
                    elif updateStatus[0] == 'success':
                        flash('%s location sucessfully set to %s' % (item[0], str(item[1])), 'success')
                    else:
                        flash('%s location update to %s failed' % (item[0], str(item[1])), 'warning')
            print 'Updated locations using :', updateList
    return redirect('/index')
