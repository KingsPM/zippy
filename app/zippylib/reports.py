#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__=="""Report Generator"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "2.0.0"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import os
import sys
import time, datetime
from math import ceil
from copy import deepcopy
from random import shuffle
import itertools
from hashlib import sha1
from collections import Counter
from . import PlateError, char_range, imageDir, githash
from .primer import parsePrimerName, PrimerPair, Primer

from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm, inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.rl_config import defaultPageSize

# Page Settings
PAGE_HEIGHT=defaultPageSize[1]; PAGE_WIDTH=defaultPageSize[0]
leftMargin = rightMargin = 2.6*cm
topMargin = 2*cm
bottomMargin = 2*cm

# MolPath version controlled document template
class MolPathTemplate(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        canvas.Canvas.__init__(self, *args, **kwargs)
        self.pages = []

    def showPage(self):
        self.pages.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        page_count = len(self.pages)
        for page in self.pages:
            self.__dict__.update(page)
            self.makePage(page_count)
            canvas.Canvas.showPage(self)
        canvas.Canvas.save(self)

    def makePage(self, page_count):
        self.setFont("Helvetica", 9)
        # header
        headerstring = "KINGS COLLEGE HOSPITAL, MOLECULAR PATHOLOGY"
        self.drawRightString(PAGE_WIDTH-rightMargin, PAGE_HEIGHT-topMargin/2.0, headerstring)
        # authorised by
        authorised = "James Bond"
        self.drawString(leftMargin, bottomMargin/2, "Authorised by %s" % authorised)
        # version
        self.drawCentredString(PAGE_WIDTH/2, bottomMargin/2, githash('zippy'))
        # page number
        page = "Page %s of %s" % (self._pageNumber, page_count)
        self.drawRightString(PAGE_WIDTH-rightMargin, bottomMargin/2, page)

# Report
class Report(object):
    def __init__(self,fi,title='This is the title',logo=None):
        # get document
        self.doc = SimpleDocTemplate(fi,pagesize=A4,
                        rightMargin=rightMargin,leftMargin=leftMargin,
                        topMargin=2*cm,bottomMargin=cm)
        # style sheets
        self.styles = getSampleStyleSheet()
        self.styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
        self.styles.add(ParagraphStyle(name='Warning', backColor = '#FFFF00', borderColor = "#000000"))
        self.styles.add(ParagraphStyle(name='tiny', fontSize=4, leading=4, alignment=TA_CENTER))
        self.styles.add(ParagraphStyle(name='small', fontSize=6, leading=6))
        # flowable elements
        self.elements = []
        # Header
        logo = Image(os.path.join(imageDir,logo),width=1.32*inch,height=0.7*inch) if logo else Paragraph('LOGO', self.styles["Heading1"])
        titl = Paragraph('%s' % title, self.styles["Heading1"])
        date = Paragraph('<font size=12>Generated on %s</font>' % time.ctime(), self.styles["Normal"])
        self.elements.append(Table([[logo,titl],['',date]],
            colWidths=[1.5*inch,4.5*inch],
            style=[ ('SPAN',(0,0),(0,1)), ('VALIGN',(0,0),(-1,-1),'BOTTOM'),
                ('VALIGN',(1,1),(-1,-1),'TOP'), ('ALIGN',(1,0),(-1,-1),'LEFT') ]))
        self.elements.append(Spacer(1, 12))
        # header fields
        TABLE_STYLE = TableStyle([
            ('ALIGN',(0,0),(-1,-1),'RIGHT'),
            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
            ('FONTSIZE',(0,1),(-1,-1),8),
            ('INNERGRID', (0,0), (1,0), 0.25, colors.black),
            ('BOX', (0,0), (1,0), 1, colors.black),
            ('BACKGROUND', (0,0), (0,0), colors.cyan),
            ('INNERGRID', (3,0), (4,0), 0.25, colors.black),
            ('BOX', (3,0), (4,0), 1, colors.black),
            ('BACKGROUND', (3,0), (3,0), colors.cyan),
            ('INNERGRID', (6,0), (7,0), 0.25, colors.black),
            ('BOX', (6,0), (7,0), 1, colors.black),
            ('BACKGROUND', (6,0), (6,0), colors.cyan)
            ])
        self.elements.append(Spacer(1, 2))
        data = [[ 'Date','','','Operator','','','Worklist','']]
        t = Table(data, \
            colWidths=[2.3*cm, 2.3*cm, 0.85*cm, 2.3*cm, 2.3*cm, 0.85*cm, 2.3*cm, 2.3*cm], rowHeights=0.6*cm)
        t.setStyle(TABLE_STYLE)
        self.elements.append(t)
        self.elements.append(Spacer(1, 12))

    def plateLayouts(self,plates):
        PLATE_STYLE = TableStyle(
            [
            ('FONTSIZE',(0,0),(-1,-1),5),
            ('VALIGN',(0,0),(-1,0),'BOTTOM'),
            ('VALIGN',(0,1),(-1,-1),'MIDDLE'),
            ('ALIGN',(0,1),(-1,-1),'LEFT'),
            ('ALIGN',(0,0),(-1,0),'CENTER'),
            ('INNERGRID', (1,1), (-1,-1), 0.25, colors.grey),
            ('BOX', (1,1), (-1,-1), 0.25, colors.grey),
            ('BACKGROUND',(0,1),(0,-1),colors.snow),
            ('BACKGROUND',(1,0),(-1,0),colors.snow),
            ('LEFTPADDING',(1,1),(-1,-1),1),
            ('RIGHTPADDING',(1,1),(-1,-1),1)
            ])

        # self.elements.append(Paragraph('Plate layouts', self.styles["Heading2"]))
        # self.elements.append(Spacer(1, 2))
        for plateNumber, data in enumerate(plates):
            self.elements.append(Paragraph('Plate%s' % str(plateNumber+1), self.styles["Heading4"]))
            data = [ [ chr(ord('A')+i) ] + r for i,r in enumerate(data) ]  # add row names
            data = [ [ str(c) if c else '' for c in range(len(data[0])) ] ] + data  # add column names
            # Make paragraphs
            for i,r in enumerate(data):
                for j,c in enumerate(r):
                    data[i][j] = Paragraph(c,self.styles['tiny']) if i and j else Paragraph(c,self.styles['small'])
            t = Table(data, colWidths=1.2*cm, rowHeights=0.5*cm)
            t.setStyle(PLATE_STYLE)
            self.elements.append(t)
            self.elements.append(Spacer(1, 12))

    def samplePrimerLists(self,s,p):
        TABLE_STYLE = TableStyle([
            ('FONTSIZE',(0,1),(-1,-1),8),  # body
            ('FONTSIZE',(0,0),(-1,0),10),  # title line
            ('VALIGN',(0,0),(-1,-1),'TOP'),
            ('ALIGN',(0,0),(-1,-1),'LEFT'),
            ('INNERGRID', (0,1), (1,len(s)), 0.25, colors.black),
            ('LINEABOVE', (0,1),(1,1),1,colors.black),
            ('BOX', (0,0), (1,len(s)), 1, colors.black),
            ('INNERGRID', (3,1), (-1,len(p)), 0.25, colors.black),
            ('LINEABOVE', (3,1),(-1,1),1,colors.black),
            ('BOX', (3,0), (-1,len(p)), 1, colors.black),
            ])
        doubleLine = ParagraphStyle('suffixes', fontSize=5, leading=5)  # suffix column
        centered = ParagraphStyle('locations', fontSize=8, leading=5, alignment=1)  # Location column
        data = [[str(len(s)),'Samples','',str(len(p)),'Primer Pair', 'Suffixes', 'Locations']]
        for i in range(max(len(s),len(p))):
            data.append([ '', s[i] if i<len(s) else '', '', ''] + \
            ([ p[i][0], Paragraph('<br/>'.join(p[i][1]),doubleLine), Paragraph(' '.join(map(str,p[i][2])),centered) ] \
            if i<len(p) else ['','','']))
        self.elements.append(Spacer(1, 2))
        t = Table(data, colWidths=[0.6*cm,5*cm,0.3*cm,0.6*cm,5.5*cm,1.6*cm,1.9*cm], rowHeights=0.6*cm)
        t.setStyle(TABLE_STYLE)
        self.elements.append(t)
        self.elements.append(Spacer(1, 12))

    def volumeLists(self,reactions,mastermix,qsolution,excess):
        # batch mix
        data = [['MasterMix', str((1.+excess)*reactions*mastermix)+' µl', '', 'Reactions', str(reactions), '' ],
            ['Q-Solution', str((1.+excess)*reactions*qsolution)+' µl', '', 'Excess', str((excess)*100)+' %', '' ],
            ['TOTAL', str((1.+excess)*reactions*(mastermix+qsolution))+' µl','','','']]
        t = Table(data, colWidths=[3*cm,2*cm,0.3*cm,3*cm,2*cm,5.2*cm], rowHeights=0.6*cm)
        t.setStyle(TableStyle([
            ('FONTSIZE',(0,0),(0,-1),10),
            ('FONTSIZE',(3,0),(3,-1),10),
            ('FONTSIZE',(1,0),(1,-1),8),
            ('FONTSIZE',(4,0),(4,-1),8),
            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
            ('ALIGN',(0,0),(-1,-1),'RIGHT'),
            ('INNERGRID', (0,0), (1,-1), 0.25, colors.black),
            ('INNERGRID', (3,0), (4,1), 0.25, colors.black),
            ('LINEABOVE', (0,-1),(1,-1), 1, colors.black),
            ('BOX', (0,0), (1,-1), 1, colors.black),
            ('BOX', (3,0), (4,1), 1, colors.black)
            ]))
        self.elements.append(Spacer(1, 12))
        self.elements.append(t)
        self.elements.append(Spacer(1, 12))

    def checkBoxes(self,titles):
        # right justified checkboxes with appropriate names
        TABLE_STYLE = TableStyle([
            ('ALIGN',(0,0),(-1,-1),'RIGHT'),
            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
            ('FONTSIZE',(0,1),(-1,-1),8),
            ('BOX', (0,0), (-1,-1), 1, colors.black),
            ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
            ('LINEABOVE', (0,1), (-1,1), 1, colors.black),
            ('BACKGROUND', (0,0), (-1,0), colors.bisque)
            ])

        self.elements.append(Paragraph('Checks', self.styles["Heading4"]))
        self.elements.append(Spacer(1, 2))
        data = [[ 'Task', 'Date', 'Operator']]
        for i in range(len(titles)):
            data.append([ titles[i], '', '' ])
        t = Table(data, colWidths=[7.5*cm,4*cm,4*cm], rowHeights=0.6*cm)
        t.setStyle(TABLE_STYLE)
        self.elements.append(t)
        self.elements.append(Spacer(1, 12))

    def build(self):
        self.doc.build(self.elements, canvasmaker=MolPathTemplate)
        #self.doc.build(self.elements, onFirstPage=myFirstPage, onLaterPages=myFirstPage)
        # self.doc.build(self.elements)

'''PCR test (PrimerPair, samplename]'''
class Test(object):
    def __init__(self,primerpair,sample,**kwargs):
        self.sample = sample
        self.primerpair = primerpair
        self.control = False

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def __str__(self):
        try:
            testString = [
                str(self.__class__),
                self.sample,
                self.primerpair,
                str(self.control),
                ','.join(map(str,self.primerpairobject.locations()))
            ]
        except:
            raise
        else:
            return ' '.join(testString)

    def __hash__(self):
        #hash((self.sample,self.primerpair,self.control))
        return hash(self.sample) ^ hash(self.primerpair) ^ hash(self.control)

    def __eq__(self,other):
        return hash(self) == hash(other)

    @property
    def primerpair(self):
        return self.primerpairobject

    @primerpair.getter
    def primerpair(self):
        return self.primerpairobject.name

    @primerpair.setter
    def primerpair(self,x):
        self.primerpairobject = x

    # returns human readable primer names (pairname and primer names if different)
    def operatorPrimerName(self):
        raise NotImplementedError

'''Worksheet data (list of tests)'''
class Worksheet(list):
    def __init__(self,elements,name='Worksheet'):
        list.__init__(self, elements)
        self.plates = []
        self.name = name
        self.date = datetime.datetime.now().isoformat()

    def __str__(self):
        return "<Worksheet (%d tests, %d variants)> " % (len(self),sum([len(e.variants) for e in self]))

    def count(self,attr,x):
        return len([ e for e in self if getattr(e,attr)==x ])

    '''generates IDT orderlist for primers without assigned storage location'''
    def orderCsv(self,fi,config):
        # get all primers without storage location
        primerpairs = set()
        for plate in self.plates:
            for r in range(len(plate.M)):
                for c in range(len(plate.M[r])):
                    if plate.M[r][c] is not None:
                        primerpairs.add(plate.M[r][c].primerpairobject)
        # print CSV with extra columns
        with open(fi,'w') as fh:
            for pp in sorted(list(primerpairs),key=lambda x: x.name):
                # add tag according to primername (check 1 fwd and 1 rev)
                tagorder = [1,0] if pp.reversed else [0,1]
                for i, p in enumerate(pp):
                    try:
                        tagseq = config['sequencetags'][p.tag]['tags'][tagorder[i]]
                    except:
                        tagseq = '' if not p.tag else p.tag+'-'
                    print >> fh, '\t'.join([ p.name, tagseq+p.seq ] + config['extracolumns'])

    '''print robot csv, (DestinationPlate,DestinationWell,SampleID,PrimerID)'''
    def robotCsv(self,fi,sep=','):
        with open(fi,'w') as fh:
            print >> fh, sep.join(['DestinationPlate','DestinationWell','PrimerID','SampleID','PrimerName'])
            for n, p in enumerate(self.plates):
                P = 'Plate'+str(n+1)
                for i,row in enumerate(p.M):
                    R = chr(ord('A')+i)
                    for j,cell in enumerate(row):
                        C = str(j+1)
                        if cell:
                            print >> fh, sep.join([P,R+C,cell.primerpairobject.uniqueid()[:10],cell.sample,cell.primerpair])

    '''assign tests to plate (smart work better only for big test sets)'''
    def fillPlates(self,size=[8,12],randomize=True,roworder='primerpair',includeSamples=True,includeControls=True):
        if randomize:
            shuffle(self)  # shuffle filling
        else:
            self = sorted(self, key=lambda x: x.sample)  # sort by sample
        # fill plates
        self.plates.append(Plate(*size))  # start plate
        for t in self:
            if t.control == includeControls or t.control != includeSamples:
                if self.plates[-1].isfull():
                    self.plates.append(Plate(*size))
                # add sample (automatically randomises controls)
                self.plates[-1].add(t,roworder)
        # compress plates
        for p in self.plates:
            p.compress(roworder)
        return

    '''add control samples'''
    def addControls(self,control="NTC"):
        controlsamples = {}
        for e in self:
            if not e.control and e.primerpair not in controlsamples.keys():
                x = deepcopy(e)
                x.control = True
                x.sample = control
                controlsamples[e.primerpair] = x
        self += controlsamples.values()
        return

    '''PDF worksheet'''
    def createWorkSheet(self,fi,primertest=False,**kwargs):
        logo = kwargs['logo'] if 'logo' in kwargs.keys() and kwargs['logo'] else None
        r = Report(fi,self.name,logo)
        # add plates
        samples, primers, plates = [], [], []
        for p in self.plates:
            s, p, m = p.platemap()  # gets samples, (pairname, (primersuffixes), (locations)), platemap
            samples += s
            primers += p  # PrimerPair Objects
            plates.append(m)
        # sample list (similar to plate order)
        sampleOrder = { s: self.plates[0]._bestRows(Test(PrimerPair([None,None],name='dummyprimer'),s),'sample')[0] \
            for s in set(samples) }
        orderedSamples = [ x[0] for x in sorted(sampleOrder.items(), key=lambda x: x[1]) ]
        # primer list (similar to plate order)
        primerOrder = { p: self.plates[0]._bestRows(Test(p,'dummy'),'primerpair')[0] \
            for p in set(primers) }
        orderedPrimers = [ (x[0].name, x[0].primerSuffixes(), tuple(x[0].locations())) for x in sorted(primerOrder.items(), key=lambda x: x[1]) ]
        # store ordered list of sample (str) and primers (primername, primersuffixes, locations)
        r.samplePrimerLists(orderedSamples,orderedPrimers)
        # reaction volume list
        r.volumeLists(sum([len(p) for p in self.plates]),kwargs['volumes']['mastermix'],kwargs['volumes']['qsolution'],kwargs['volumes']['excess'])
        # add checkboxes
        allTested = sum([ p.locations().count(None) for p in primers])==0
        checkTasks = ['New primers ordered', 'Plate orientation checked', 'Primer checked and storage assigned'] if primertest \
            else ['Plate orientation checked'] if allTested \
            else ['New primers tested', 'Plate orientation checked']
        r.checkBoxes(checkTasks)
        # plate layout
        r.plateLayouts(plates)
        # build pdf
        r.build()

    '''tube Labels'''
    def tubeLabels(self,fi='/dev/null',tags={}):
        # detects collisions, run with empty output to validate
        digests = {}  # digest -> name
        with open(fi,'w') as fh:
            print >> fh, '~SD30'  # darkness to maximum
            for n, p in enumerate(self.plates):
                for i,row in enumerate(p.M):
                    for j,cell in enumerate(row):
                        if cell:
                            d = cell.primerpairobject.uniqueid()[:10]  # uniqueid
                            if d in digests.keys():
                                try:
                                    assert cell.primerpair == digests[d]
                                except:
                                    raise Exception('BarcodeCollision')
                            else:
                                digests[d] = cell.primerpair
                            # get tag name
                            tagstring = ' '.join(set([ tags[x.tag]['name'] if x.tag in tags.keys() else x.tag \
                                for x in cell.primerpairobject ]))
                            print >> fh, "^XA"  # start label
                            print >> fh, "^FO20,25^AB^FD{}^FS".format(self.date[:self.date.rfind('.')])  # date to the second
                            print >> fh, "^FO20,45^AB,25^FD{}^FS".format(cell.primerpair)  # primer name
                            print >> fh, "^FO20,75^AB^FD{}^FS".format(tagstring)  # primer name
                            print >> fh, "^FO20,95^BY1.5^BCN,80,Y,N,N^FD{}^FS".format(d)  # uniqueid
                            print >> fh, "^XZ"  # end label


class Plate(object):
    def __init__(self,rows,columns):
        self.M = [ [ None for j in range(columns) ] for i in range(rows) ]
        self.nrow = rows
        self.ncol = columns

    def __repr__(self):
        return '\n'.join([ str(i)+': '+str(r) for i,r in enumerate(self.M)])

    def __len__(self):
        return sum([ len([ f for f in r if f]) for r in self.M ])

    def platemap(self):
        s, p = set(), set()
        for r in range(len(self.M)):
            for c in range(len(self.M[r])):
                if self.M[r][c] is not None:
                    s.add(self.M[r][c].sample)
                    p.add(self.M[r][c].primerpairobject)
        # # generate plate map
        trunc = lambda x: '<br/>'.join([ x.sample[:9]+'...' if len(x.sample)>12 else x.sample,
            x.primerpair[:9]+'...' if len(x.primerpair)>12 else x.primerpair]) if x is not None else "EMPTY"
        pm = [ map(trunc,r) for r in self.M ]
        # return values
        return sorted(list(s)), sorted(list(p),key=lambda x: x.name), pm

    def bestRowNumbers(self,val,attr):
        rowcounts = { i: len([ e for e in r if e is not None and getattr(e,attr) == val ]) \
            for i, r in enumerate(self.M) }
        return [ x[0] for x in sorted(rowcounts.items(), key=lambda x: x[1], reverse=True) ]

    def _bestRows(self,x,rvalue):
        if x.control:  # try to spread across rows
            rowcounts = { i: \
                len([ e for e in r if e is not None and getattr(e,rvalue) != getattr(x,rvalue) ]) - \
                len([ e for e in r if e is not None and getattr(e,rvalue) == getattr(x,rvalue) ]) \
                for i, r in enumerate(self.M) }
        else:  # use same row if possible
            rowcounts = { i: \
                len([ e for e in r if e is not None and getattr(e,rvalue) == getattr(x,rvalue) ]) - \
                len([ e for e in r if e is not None and getattr(e,rvalue) != getattr(x,rvalue) ]) \
                for i, r in enumerate(self.M) }
        return [ x[0] for x in sorted(rowcounts.items(), key=lambda x: x[1], reverse=True) ]

    def isfull(self):
        return True if len(self)>=len(self.M)*len(self.M[0]) else False

    # s with row priority only (corresponds to sample)
    def add(self,t,rvalue='primerpair'):
        # get best rows
        try:
            bestRows = self._bestRows(t,rvalue)  # randomises for controls
            assert len(bestRows)>0
        except AssertionError:
            raise PlateError('Plate Full')
        except:
            raise
        else:
            # get free columns
            availableWells = [ r for r in itertools.product(bestRows,range(self.ncol)) if self.M[r[0]][r[1]] is None ]
        # add to first available well
        self.M[availableWells[0][0]][availableWells[0][1]] = t

    def compress(self,rvalue):
        # find maximum row count
        maxColumn = [ int(len(self)/self.nrow) + 1 if r < len(self) % self.nrow else int(len(self)/self.nrow) \
            for r in range(self.nrow) ]
        # start with longest row (and only rows which exceed maxColumns)
        rowOrder = { i: r.count(None) for i,r in enumerate(self.M) if r.count(None) < self.ncol-maxColumn[i] }
        longRows = [ i for i,r in sorted(rowOrder.items(), key=lambda x: x[1], reverse=False) ]
        for r in longRows:
            # find best row (sum of empty and own)
            for c in range(maxColumn[r],self.ncol):
                if self.M[r][c] is not None:
                    bestRows = [ br for br in self._bestRows(self.M[r][c],rvalue) \
                        if br!=r and self.M[br].count(None) > self.ncol-maxColumn[br] ]
                    # find leftmost free position
                    for bc in range(self.ncol):
                        if self.M[bestRows[0]][bc] is None:
                            self.M[r][c], self.M[bestRows[0]][bc] = self.M[bestRows[0]][bc], self.M[r][c]  # switch position
                            break
        return
