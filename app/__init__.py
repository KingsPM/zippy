#!/usr/bin/env python

__doc__=="""Zippy"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "1.1"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import os
from flask import Flask, render_template, request, redirect
from celery import Celery
from werkzeug.utils import secure_filename

app = Flask(__name__)

# app.config.from_object('config')

from app import views
