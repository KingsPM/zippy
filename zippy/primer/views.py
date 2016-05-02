# -*- coding: utf-8 -*-
"""Public section, including homepage and signup."""
from flask import Blueprint, flash, redirect, render_template, request, url_for
from flask_login import login_required

from zippy.extensions import login_manager
from zippy.primer.forms import LocusForm
from zippy.user.models import User
from zippy.utils import flash_errors

blueprint = Blueprint('primer', __name__, static_folder='../static')

@blueprint.route('/design', methods=['GET'])
@login_required
def design():
    """Design Primer"""
    form = LocusForm(request.form)
    return render_template('primer/design.html', form=form)

# @blueprint.route('/worklist', methods=['GET'])
# def batch():
#     """Create worksheet for batch validation"""
#     return render_template('primer/batch.html', form=form)
#
# # admin methods
# @blueprint.route('/import', methods=['GET'])
# def design():
#     """Import from list"""
#     return render_template('public/home.html', form=form)
#
#
# @blueprint.route('/jobs', method=['GET','POST']
# def jobs():
#     """Jobs that have been requested by user"""
