#!/usr/bin/env python

import sys
import os
import re
import json
import hashlib
import subprocess
from flask import Flask, render_template, request, redirect, send_from_directory
from zippy import app
from .. import ascii_encode_dict
from ..primer import Location
from ..database import PrimerDB

# RESTful API for zippy (LOVD/SNPpy connections)

@app.route('/update_location/<primername>/<int:vessel>/<well>', methods=['PUT','GET','POST'])
def update_Location2(primername,vessel,well):
    # read config
    with open(app.config['CONFIG_FILE']) as conf:
        config = json.load(conf, object_hook=ascii_encode_dict)
        db = PrimerDB(config['database'])
    # run zippy and render
    updateStatus = updateLocation(primername, Location(vessel, well), db, False)
    return render_template('location_updated.html', status=updateStatus)
