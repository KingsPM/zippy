#/usr/bin/env python

import os
from flask import Flask, render_template, request, redirect
from celery import Celery
from werkzeug.utils import secure_filename

app = Flask(__name__)

# app.config.from_object('config')

from app import views
