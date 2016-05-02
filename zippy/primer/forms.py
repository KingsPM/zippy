# -*- coding: utf-8 -*-
"""Public forms."""
from flask_wtf import Form
from wtforms import PasswordField, StringField, DecimalField, BooleanField
from wtforms.validators import DataRequired, Length

from zippy.primer.models import Locus

class LocusForm(Form):
    """DesignPrimers"""
    chrom = StringField('Chromosome', validators=[DataRequired(), Length(min=1, max=25)])
    chromStart = DecimalField('Start', validators=[DataRequired()])
    chromEnd = DecimalField('End', validators=[DataRequired()])
    reverse = BooleanField('Reverse', validators=[])

    def __init__(self, *args, **kwargs):
        """Create instance."""
        super(LocusForm, self).__init__(*args, **kwargs)
        self.locus = None

    def validate(self):
        """Validate the form."""
        initial_validation = super(RegisterForm, self).validate()
        if not initial_validation:
            return False
        user = User.query.filter_by(username=self.username.data).first()
        if user:
            self.username.errors.append('Username already registered')
            return False
        user = User.query.filter_by(email=self.email.data).first()
        if user:
            self.email.errors.append('Email already registered')
            return False
        return True
