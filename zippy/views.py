#!/usr/local/env python

import sys
import os
import re
import json
import bcrypt
import hashlib
import subprocess
from functools import wraps
from flask import Flask, render_template, request, redirect, send_from_directory, session, flash, url_for
from celery import Celery
from werkzeug.utils import secure_filename
from . import app
from .zippy import zippyBatchQuery, zippyPrimerQuery, updateLocation, searchByName, updatePrimerName, updatePrimerPairName, blacklistPair, deletePair, readprimerlocations
from .zippylib import ascii_encode_dict
from .zippylib.primer import Location
from .zippylib.database import PrimerDB

app.config['ALLOWED_EXTENSIONS'] = set(['txt', 'batch', 'vcf', 'bed', 'csv', 'tsv'])
app.secret_key = 'Zippy is the best handpuppet out there'
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['DOWNLOAD_FOLDER'] = 'results'
app.config['CONFIG_FILE'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zippy.json')
# read password (SHA1 hash, not the safest)
with open(app.config['CONFIG_FILE']) as conf:
    config = json.load(conf, object_hook=ascii_encode_dict)
    app.config['PASSWORD'] = config['password']

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

def login_required(func):
    @wraps(func)
    def wrap(*args, **kwargs):
        if 'logged_in' in session:
            return func(*args, **kwargs)
        else:
            flash('Authentication required!', 'warning')
            return redirect(url_for('login'))
    return wrap

@app.route('/')
@app.route('/index')
@login_required
def index():
    return render_template('index.html',designtiers=config['design']['tiers'])

# simple access control (login)
@app.route('/login', methods=['GET', 'POST'])
def login():
    error = None
    if request.method == 'POST':
        if bcrypt.hashpw(request.form['password'].rstrip().encode('utf-8'), app.config['PASSWORD']) == app.config['PASSWORD']:
            session['logged_in'] = True
            return redirect(url_for('index'))
        else:
            error = 'Wrong password. Please try again.'
    return render_template('login.html', error=error)

@app.route('/logout')
def logout():
    session.pop('logged_in', None)
    return redirect(url_for('login'))


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
    uploadFile = request.files['variantTable']
    uploadFile2 = request.files['missedRegions']
    uploadFile3 = request.files['singleGenes']
    tiers = map(int,request.form.getlist('tiers'))
    predesign = request.form.get('predesign')
    design = request.form.get('design')
    outfile = request.form.get('outfile')

    # save files
    uploadedFiles = []
    for uf in [uploadFile, uploadFile2, uploadFile3]:
        if uf and allowed_file(uf.filename):
            uploadedFiles.append(os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(uf.filename)))
            uf.save(uploadedFiles[-1])
            print >> sys.stderr, "file saved to %s" % uploadedFiles[-1]

    if uploadedFiles:
        # open config file and database
        with open(app.config['CONFIG_FILE']) as conf:
            config = json.load(conf, object_hook=ascii_encode_dict)
            db = PrimerDB(config['database'],dump=config['ampliconbed'])

        # create output folder
        filehash = hashlib.sha1(''.join([ open(uf).read() for uf in uploadedFiles ])).hexdigest()
        downloadFolder = os.path.join(app.config['DOWNLOAD_FOLDER'], filehash)
        subprocess.check_call(['mkdir', '-p', downloadFolder], shell=False)

        # run Zippy to design primers
        shortName = os.path.splitext(os.path.basename(uploadedFiles[0]))[0]
        downloadFile = os.path.join(downloadFolder, outfile) if outfile else os.path.join(downloadFolder, shortName)
        arrayOfFiles, missedIntervalNames = zippyBatchQuery(config, uploadedFiles, design, downloadFile, db, predesign, tiers)
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
    tiers = map(int,request.form.getlist('tiers'))
    gap = request.form.get('gap')
    store = request.form.get('store')

    print >> sys.stderr, 'tiers', tiers
    print >> sys.stderr, 'locus', locus
    print >> sys.stderr, 'gap', gap

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
            db = PrimerDB(config['database'],dump=config['ampliconbed'])
        # run Zippy
        primerTable, resultList, missedIntervals = zippyPrimerQuery(config, target, design, None, db, store, tiers, gap)

        print >> sys.stderr, primerTable

        # get missed and render template
        missedIntervalNames = []
        for interval in missedIntervals:
            missedIntervalNames.append(interval.name)
        return render_template('/adhoc_result.html', primerTable=primerTable, resultList=resultList, missedIntervals=missedIntervalNames)
    else:
        print >> sys.stderr, "no locus or file given"
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
        print >> sys.stderr, 'Please fill in all fields (PrimerName VesselNumber Well)'
        return render_template('location_updated.html', status=None)
    # read config
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'],dump=config['ampliconbed'])
    # run zippy and render
    updateStatus = updateLocation(primername, loc, db, force)
    return render_template('location_updated.html', status=updateStatus)

@app.route('/select_pair_to_update/<pairName>')
def pair_to_update(pairName):
    return render_template('update_pair.html', pairName=pairName)

@app.route('/update_pair_name/<pairName>', methods=['POST'])
def update_pair_name(pairName):
    newName = request.form.get('name')
    if newName == pairName:
        flash('New name is the same as current', 'warning')
        return render_template('update_pair.html', pairName=pairName)
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'],dump=config['ampliconbed'])
        if updatePrimerPairName(pairName, newName, db):
            flash('Pair "%s" renamed "%s"' % (pairName, newName), 'success')
        else:
            flash('Pair renaming failed', 'warning')
    return render_template('update_pair.html', pairName=newName)

@app.route('/select_primer_to_update/<primerName>/<primerLoc>')
def primer_to_update(primerName, primerLoc):
    primerInfo = primerName + '|' + primerLoc
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
            print >> sys.stderr, 'Please fill in all fields (PrimerName VesselNumber Well)'
            return render_template('location_updated.html', status=None)
        with open(app.config['CONFIG_FILE']) as conf:
            config = json.load(conf, object_hook=ascii_encode_dict)
            db = PrimerDB(config['database'],dump=config['ampliconbed'])
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
        splitInfo = primerInfo.split('|')
        primerName = splitInfo[0]
        primerLoc = splitInfo[1]
        return render_template('update_location_from_table.html', primerName=primerName, primerLoc=primerLoc)

@app.route('/select_primer_to_rename/<primerName>/<primerLoc>', methods=['POST'])
def primer_to_rename(primerName, primerLoc):
    newName = request.form.get('name')
    print >> sys.stderr, primerName
    print >> sys.stderr, primerLoc
    print >> sys.stderr, newName
    primerInfo = primerName + '|' + primerLoc + '|' + newName
    print >> sys.stderr, primerInfo
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
        db = PrimerDB(config['database'],dump=config['ampliconbed'])
        if updatePrimerName(currentName, newName, db):
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
        db = PrimerDB(config['database'],dump=config['ampliconbed'])
        searchResult = searchByName(searchName, db)
    return render_template('searchname_result.html', searchResult=searchResult, searchName=searchName)

@app.route('/blacklist_pair/<pairname>', methods=['POST'])
def blacklist_pair(pairname):
    print >> sys.stderr, 'This is the pairname: ' + pairname
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'],dump=config['ampliconbed'])
        blacklisted = blacklistPair(pairname, db)
        for b in blacklisted:
            flash('%s added to blacklist' % (b,), 'success')
    return redirect('/search_by_name/')

@app.route('/delete_pair/<pairname>', methods=['POST'])
def delete_pair(pairname):
    print >> sys.stderr, 'This is the pairname: ' + pairname
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'],dump=config['ampliconbed'])
        deleted = deletePair(pairname, db)
        for d in deleted:
            flash('%s deleted' % (d,), 'success')
    return redirect('/search_by_name/')

@app.route('/upload_batch_locations/', methods=['POST'])
def upload_samplesheet():
    if request.method == 'POST':
        locationsheet = request.files['locationsheet']
        if not locationsheet or not locationsheet.filename.endswith('.csv'):
            flash('Not a CSV file. Please try again','warning')
        else:
            filename = secure_filename(locationsheet.filename)
            saveloc = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            locationsheet.save(saveloc)
            updateList = readprimerlocations(saveloc)
            with open(app.config['CONFIG_FILE']) as conf:
                config = json.load(conf, object_hook=ascii_encode_dict)
                db = PrimerDB(config['database'],dump=config['ampliconbed'])
                for item in updateList:
                    updateStatus = updateLocation(item[0], item[1], db, True) # Force is set to True, will force primers into any occupied locations
                    if updateStatus[0] == 'occupied':
                        flash('Location already occupied by %s' % (' and '.join(updateStatus[1])), 'warning')
                    elif updateStatus[0] == 'success':
                        flash('%s location sucessfully set to %s' % (item[0], str(item[1])), 'success')
                    else:
                        flash('%s location update to %s failed' % (item[0], str(item[1])), 'warning')
            print >> sys.stderr, 'Updated locations using :', updateList
    return redirect('/index')
