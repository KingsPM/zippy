#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__=="""Report Generator"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "2.3.4"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import os
import sys
import time, datetime
from math import ceil
from copy import deepcopy
from random import shuffle
from itertools import product, groupby, chain
from functools import partial
from hashlib import sha1
from collections import Counter
from . import PlateError, char_range, imageDir, githash
from .primer import parsePrimerName, PrimerPair, Primer
from urllib import unquote

from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm, inch
from reportlab.platypus import SimpleDocTemplate, BaseDocTemplate, Frame, Paragraph, Spacer, Image, Table, TableStyle, PageBreak, NextPageTemplate, PageTemplate
from reportlab.platypus.flowables import KeepTogether
from reportlab.rl_config import defaultPageSize


# Page Settings
leftMargin = rightMargin = topMargin = bottomMargin = 2*cm

# MolPath version controlled document template
class MolPathTemplate(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        self.auth = ""
        self.docid = ""
        self.site = ""
        # document authorisation, header and worklist
        if 'auth' in kwargs.keys():
            self.auth = kwargs['auth']
            del kwargs['auth']
        if 'docid' in kwargs.keys():
            self.docid = kwargs['docid']
            del kwargs['docid']
        if 'site' in kwargs.keys():
            self.site = kwargs['site']
            del kwargs['site']
        if 'worklist' in kwargs.keys():
            self.worklist = kwargs['worklist']
            del kwargs['worklist']
        # build canvas
        canvas.Canvas.__init__(self, *args, **kwargs)
        self.pages = []

    def showPage(self):
        self.pages.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        page_count = len(self.pages)
        for page in self.pages:
            self.__dict__.update(page)
            self.makePage(page_count,page['_pagesize'])
            canvas.Canvas.showPage(self)
        canvas.Canvas.save(self)

    def makePage(self, page_count, pagesize):
        self.setFont("Helvetica", 9)
        # docid
        if self.site:
            self.drawRightString(pagesize[0]-rightMargin, pagesize[1]-topMargin/2.0, self.site)
        # sitestring
        if self.docid:
            self.drawString(leftMargin, pagesize[1]-topMargin/2.0, self.docid)
        # authorised by
        if self.auth:
            self.drawString(leftMargin, bottomMargin/2, "Authorised by %s" % self.auth)
        # version
        self.drawCentredString(pagesize[0]/2, bottomMargin/2, githash('zippy'))
        # page number
        page = "%sPage %s of %s" % ((self.worklist+' / ') if self.worklist else '', self._pageNumber, page_count)
        self.drawRightString(pagesize[0]-rightMargin, bottomMargin/2, page)


# Report
class Report(object):
    def __init__(self,fi,title='This is the title',logo=None,site='',auth='',docid='',worklist=''):
        # site and auth
        self.site = site
        self.auth = auth
        self.docid = docid
        self.worklist = worklist
        # get document
        self.doc = BaseDocTemplate(fi,
                      rightMargin=rightMargin,
                      leftMargin=leftMargin,
                      topMargin=topMargin,
                      bottomMargin=bottomMargin,
                      leftPadding = 0,
                      rightPadding = 0,
                      topPadding = 0,
                      bottomPadding = 0,
                      showBoundary=0)

        # style sheets
        self.styles = getSampleStyleSheet()
        self.styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
        self.styles.add(ParagraphStyle(name='Warning', backColor = '#FFFF00', borderColor = "#000000"))
        self.styles.add(ParagraphStyle(name='tiny', fontSize=4, leading=4, alignment=TA_CENTER))
        self.styles.add(ParagraphStyle(name='small', fontSize=6, leading=6))
        self.styles.add(ParagraphStyle(name='normal', fontSize=8, leading=6))
        self.styles.add(ParagraphStyle(name='big', fontSize=10, leading=6))
        self.styles.add(ParagraphStyle(name='huge', fontSize=12, leading=6))

        # add Page templates
        frame1 = Frame(self.doc.leftMargin, self.doc.bottomMargin,
            self.doc.width, self.doc.height)
        frame2 = Frame(self.doc.leftMargin, self.doc.bottomMargin,
            self.doc.height, self.doc.width)
        ptemplate = PageTemplate(id='portrait',frames =[frame1], onPage=lambda canvas, doc: canvas.setPageSize(A4))
        ltemplate = PageTemplate(id='landscape',frames =[frame2], onPage=lambda canvas, doc: canvas.setPageSize(landscape(A4)))
        self.doc.addPageTemplates([ptemplate, ltemplate])

        # flowable elements
        self.elements = []
        # Header
        logo = Image(os.path.join(imageDir,logo),width=1.32*inch,height=0.7*inch) if logo else Paragraph('LOGO', self.styles["Heading1"])
        titl = Paragraph('%s' % title, self.styles["Heading2"])
        date = Paragraph('<font size=10>Generated on %s</font>' % time.ctime(), self.styles["Normal"])
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
        data = [[ 'Date','','','Operator','','','Worklist',self.worklist]]
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
            self.elements.append(KeepTogether(t))
            self.elements.append(Spacer(1, 12))

    def setNextPageTemplate(self,template):
        self.elements.append(NextPageTemplate(template))

    # generic dataframe style table (relative column width, pagesize dependent, full width)
    def genericTable(self, data, landscape=False, rowShading=False, header=True, tableTitle=None, relativeColWidth=None, rowHeight=0.5*cm, mergeColumnFields=None):
        # global stylesheet
        stylesheet = [
            ('ALIGN',(0,0),(-1,-1),'LEFT')
        ]
        # header dependent styles
        if header:
            stylesheet += [
                ('FONTSIZE',(0,0),(-1,0),9),
                ('BACKGROUND', (0,0), (-1,0), colors.bisque),  # header color
                ('INNERGRID', (0,0), (-1,0), 0.25, colors.black), # header inner grid
                ('BOX', (0,0), (-1,0), 1, colors.black)  # header box
            ]
        # Data cell styles
        firstData = (0,1) if header else (0,0)
        stylesheet += [
            ('FONTSIZE',firstData, (-1,-1),8),
            ('VALIGN', firstData, (-1,-1), 'MIDDLE') # middle cell alignment for data rows
        ]
        # column create boxes (merge cells?)
        groupBoxes = []  # rows to box
        if mergeColumnFields is None:
            mergeColumnFields = list(range(len(data[0])))
        # get column numbers for unmerged fields (draw individual cells on those)
        for c in set(range(len(data[0]))) - set(mergeColumnFields):
            stylesheet += [
                ('INNERGRID', (c,firstData[1]), (c,-1), 0.25, colors.black),
                ('BOX', (c,firstData[1]), (c,-1), 0.25, colors.black)
            ]
        # MERGE CELLS
        for i, c in enumerate(mergeColumnFields):
            rowRanges = range(firstData[1],len(data[firstData[1]:]))  # groups that are treated seperately
            offsetRow = firstData[1]  # define first row offset
            for k,g in groupby([ data[r][c] for r in range(firstData[1],len(data)) ]):
                groupSize = len(list(g))
                startRow, endRow = offsetRow, offsetRow + groupSize - 1  # set start,end rows for group
                offsetRow = endRow + 1  # update offsetRow
                # add Box for cell groups
                stylesheet += [ ('BOX', (c,startRow), (c,endRow), 0.25, colors.black) ]
                # create box if first request column to merge cells in
                if i==0:
                    groupBoxes += [ ('BOX', (0,startRow), (-1,endRow), 1, colors.black) ]  # data box
                # merge cells
                mergeData = [ data[r][c] for r in range(startRow,offsetRow) ]
                assert len(set(mergeData))==1
                for n, r in enumerate(range(startRow,offsetRow)):
                    # remove duplicate fields
                    if n:
                        data[r][c] = ''
        stylesheet += groupBoxes  # add group boxes
        if rowShading:
            for r in range(firstData[1],len(data)):
                stylesheet += [('BACKGROUND',(0,r),(-1,r),colors.lightgrey)] if r % 2 == 0 else [('BACKGROUND',(0,r),(-1,r),colors.snow)]
        # add title
        if tableTitle:
            self.elements.append(Paragraph(tableTitle, self.styles["Heading4"]))
            self.elements.append(Spacer(1, 2))
        # define columnWidths
        frameWidth = self.doc.height if landscape else self.doc.width
        if relativeColWidth:
            # define relative columns widths
            colw = [ x*frameWidth/sum(relativeColWidth) for x in relativeColWidth ]
        else:
            colw = frameWidth/float(len(data[0]))

        # create table
        rowHeights = [rowHeight*1.2] + [rowHeight]*(len(data)-1) if header else [rowHeight] * len(data)
        t = Table(data, colWidths=colw, rowHeights=rowHeights, repeatRows=1 if header else 0)
        t.setStyle(TableStyle(stylesheet))
        self.elements.append(t)
        self.elements.append(Spacer(1, 12))

    def pageBreak(self):
        self.elements.append(PageBreak())

    def samplePrimerLists(self,s,p,counts=Counter()):
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
        centeredsmall = ParagraphStyle('locations', fontSize=6, leading=6, alignment=1)  # Location column
        data = [[str(len(s)),'Samples','',str(len(p)),'Primer Pairs', 'Suffixes', 'Locations']]
        for i in range(max(len(s),len(p))):
            if i<len(p):
                if any(p[i][2]):
                    locationString = ' '.join(map(str,p[i][2]))
                    locationParagraph = Paragraph(locationString, centered if len(locationString) < 10 else centeredsmall)
                else:
                    locationParagraph = Paragraph(' ',centered)
            data.append(([ counts[s[i]] if counts else '', s[i] ] if i<len(s) else ['','']) + [''] + \
                ([ counts[p[i][0]] if counts else '', p[i][0], Paragraph('<br/>'.join(p[i][1]),doubleLine), locationParagraph ] if i<len(p) else ['','','','']))
        self.elements.append(Spacer(1, 2))
        t = Table(data, colWidths=[0.6*cm,5*cm,0.3*cm,0.6*cm,5.3*cm,1.6*cm,2.1*cm], rowHeights=0.6*cm)
        t.setStyle(TABLE_STYLE)
        self.elements.append(t)
        self.elements.append(Spacer(1, 12))

    def volumeLists(self,reactions,mastermix,qsolution,excess,program):
        # batch mix
        data = [['Reagent','Quantity','LOT','Expiry','','Reactions', str(reactions) ],
            ['MasterMix', str((1.+excess)*reactions*mastermix)+' µl', '', '', '', 'Excess', str((excess)*100)+' %' ],
            ['Q-Solution', str((1.+excess)*reactions*qsolution)+' µl', '', '', '', 'PCR Program', program ],
            ['TOTAL', str((1.+excess)*reactions*(mastermix+qsolution))+' µl', '', '', '', 'PCR Block','']]
        t = Table(data, colWidths=[2.5*cm,2.5*cm,2.5*cm,2.5*cm,0.3*cm,2.7*cm,2.5*cm], rowHeights=0.6*cm)
        t.setStyle(TableStyle([
            ('FONTSIZE',(0,1),(0,-1),10),
            ('FONTSIZE',(1,0),(4,0),10),
            ('FONTSIZE',(1,1),(4,-1),8),
            ('FONTSIZE',(5,0),(5,-1),10),
            ('FONTSIZE',(6,0),(6,-1),8),
            ('INNERGRID', (0,0), (3,-1), 0.25, colors.black),
            ('INNERGRID', (5,0), (6,-1), 0.25, colors.black),
            ('LINEABOVE', (0,1),(3,1), 1, colors.black),
            ('BACKGROUND',(2,-1),(3,-1),colors.lightgrey),
            ('BACKGROUND', (0,0), (3,0), colors.bisque),
            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
            ('ALIGN',(0,0),(-1,-1),'RIGHT'),
            ('BOX', (0,0), (3,-1), 1, colors.black),
            ('BOX', (5,0), (6,-1), 1, colors.black)
            ]))
        self.elements.append(Spacer(1, 12))
        self.elements.append(KeepTogether(t))
        self.elements.append(Spacer(1, 12))

    def checkBoxes(self,title='Checks',table=[],tableHeader=['Task','Date','Checker'],tickbox=[],tickboxNames=['YES','NO'],textLines={}):
        # title
        if title:
            self.elements.append(Paragraph(title, self.styles["Heading4"]))
            self.elements.append(Spacer(1, 2))
        # Task table
        if table:
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
            data = [tableHeader]
            for i in range(len(table)):
                data.append([ table[i], '', '' ])
            t = Table(data, colWidths=[7.5*cm,4*cm,4*cm], rowHeights=0.6*cm)
            t.setStyle(TABLE_STYLE)
            self.elements.append(KeepTogether(t))
            self.elements.append(Spacer(1, 6))

        if tickbox:
            stylesheet = [
                ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
                ('ALIGN',(1,0),(-1,-1),'LEFT'),
                ('ALIGN',(0,0), (0,-1),'RIGHT'),
                ('FONTSIZE',(0,0),(-1,-1),9)
            ]
            # set column width
            colWidths = [3*cm]
            for t in tickboxNames:
                colWidths += [0.6*cm,1.2*cm]
            # set row heights
            rowHeights = []
            # compile data
            data = []
            for i, b in enumerate(tickbox):
                data.append([ b ])
                for j, t in enumerate(tickboxNames):
                    data[-1] += [ '', t ]
                    stylesheet += [ ('BOX', (2*j+1,i*2), (2*j+1,i*2), 1, colors.black) ]
                data.append(['']*len(data[-1]))  # empty line
                rowHeights += [0.6*cm,0.3*cm]
            t = Table(data, colWidths=colWidths, rowHeights=rowHeights)
            t.setStyle(TableStyle(stylesheet))
            self.elements.append(KeepTogether(t))
            self.elements.append(Spacer(1, 6))

        if textLines:
            stylesheet = [
                ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
                ('ALIGN',(1,0),(-1,-1),'LEFT'),
                ('ALIGN',(0,0), (0,-1),'RIGHT'),
                ('FONTSIZE',(0,0),(0,-1),12),
                ('LINEBELOW', (1,0), (1,-1), 1.0, colors.black)
            ]
            data = []
            for t,l in textLines.items():
                for i in range(l):
                    data.append([ t+':' if i==0 else '', ''])
            t = Table(data, colWidths=[4*cm,11.5*cm], rowHeights=0.6*cm)
            t.setStyle(TableStyle(stylesheet))
            self.elements.append(KeepTogether(t))
            self.elements.append(Spacer(1, 6))
        self.elements.append(Spacer(1, 6))

    def build(self):
        self.doc.build(self.elements, \
            canvasmaker=partial(MolPathTemplate, site=self.site, auth=self.auth, docid=self.docid, worklist=self.worklist))


'''PCR test (PrimerPair (with variants), samplename]'''
class Test(object):
    def __init__(self,primerpair,sample,**kwargs):
        self.sample = sample
        self.primerpair = primerpair
        self.control = False

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def __str__(self):
        return "{} {} {}".format(self.sample, self.primerpair, ','.join(map(repr,self.primerpairobject.variants)))

    def __hash__(self):
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
                for col in range(p.ncol):
                    C = str(col+1)
                    for row in range(p.nrow):
                        R = chr(ord('A')+row)
                        cell = p.M[row][col]
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

    ''' count reactions '''
    def reactionCount(self):
        reactions = Counter()
        for plate in self.plates:
            for test in plate.testList():
                reactions[test.sample] += 1
                reactions[test.primerpair] += 1
        return reactions

    '''PDF worksheet'''
    def createWorkSheet(self,fi,primertest=False,worklist='',**kwargs):
        logo = kwargs['logo'] if 'logo' in kwargs.keys() and kwargs['logo'] else None
        site = kwargs['site'] if 'site' in kwargs.keys() and kwargs['site'] else None
        auth = kwargs['auth'] if 'auth' in kwargs.keys() and kwargs['auth'] else None
        docid = kwargs['docid'] if 'docid' in kwargs.keys() and kwargs['docid'] else None
        r = Report(fi,title=self.name,logo=logo,site=site,auth=auth,docid=docid,worklist=worklist)
        # add plates
        samples, primers, plates = [], [], []
        for plate in self.plates:
            s, p, m = plate.platemap()  # gets samples, (pairname, (primersuffixes), (locations)), platemap
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
        r.samplePrimerLists(orderedSamples,orderedPrimers,counts=self.reactionCount())
        # reaction volume list
        r.volumeLists(sum([len(p) for p in self.plates]),kwargs['volumes']['mastermix'],kwargs['volumes']['qsolution'],kwargs['volumes']['excess'],kwargs['volumes']['program'])
        # add checkboxes
        checkTasks = ['New primers ordered', 'Plate orientation checked', 'Primer checked and storage assigned'] if primertest \
            else ['Plate orientation checked', 'DNA barcodes relabeled']
        r.checkBoxes(title='',table=checkTasks)
        # plate layout
        r.plateLayouts(plates)
        # print result table
        if primertest:
            fields = [[ 'Primer Pair', 'Amplicon Size', 'Result']]
            for i,t in enumerate([ x for x in sorted(self,key=lambda x: (x.sample,x.primerpair)) if not x.control ]):
                fields += [[ t.primerpair, str(t.primerpairobject.targetLength(includePrimers=True))+' bp', '']]
            # create result table
            r.pageBreak()
            r.genericTable(fields,tableTitle='Results',landscape=False, mergeColumnFields=[], relativeColWidth=[2,1,3])
            # add checkboxes
            checkTasks = ['Primer checked and storage assigned']
            r.checkBoxes(title='',table=checkTasks)
        else:
            fields = [['DNA #', 'Patient Name', 'Primer Pair', 'Variant', 'Zygosity', 'Result', 'Check']]
            for i,t in enumerate([ x for x in sorted(self,key=lambda x: (x.sample,x.primerpair)) if not x.control ]):
                for v in t.primerpairobject.variants:
                    fields += [[ t.sample, '', t.primerpair, ' '.join(unquote(v.name).split(',')[:-1]), unquote(v.name).split(',')[-1], '', '']]
            # create result table
            r.setNextPageTemplate('landscape')
            r.pageBreak()
            r.genericTable(fields,tableTitle='Results',landscape=True, mergeColumnFields=[0,1],relativeColWidth=[0.8,0.8,0.8,2.6,0.5,2.2,0.4])
            # add checkboxes
            r.checkBoxes(title='',table=['Primary Reporter', 'Secondary Reporter'],tableHeader=['Reporter','Date','Initial'],
                tickbox=['Unmatched Sample Check', 'Control Check'], tickboxNames=['YES','NO'],
                textLines={'Comments': 3})
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
                            # barcode id (with collision check, as trucated 32 byte string)
                            d = cell.primerpairobject.uniqueid()[:10]  # truncated uniqueid (1,099,511,627,776)
                            if d in digests.keys():
                                try:
                                    assert cell.primerpair == digests[d]
                                except:
                                    raise Exception('BarcodeCollision')
                                else:
                                    continue  # dont print same barcode multiple times
                            else:
                                digests[d] = cell.primerpair
                            # Location string
                            locations = ' '.join([ str(l) if l else '' for l in cell.primerpairobject.locations() ])
                            # get tag name
                            tagstring = '/'.join(set([ x.tag for x in cell.primerpairobject ]))
                            print >> fh, "^XA"  # start label
                            print >> fh, "^PR1,A,A"  # slower print speed
                            print >> fh, "^FO20,50^AB^FD{}  {}^FS".format(self.date[:self.date.rfind('.')],locations)  # date and location
                            print >> fh, "^FO20,70^AB,25^FD{}^FS".format(cell.primerpair)  # primer name
                            print >> fh, "^FO20,100^BY1.5^BCN,80,Y,N,N^FD{}^FS".format(d)  # Barcode uniqueid
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

    def testList(self):
        return [ t for t in chain.from_iterable(zip(*self.M)) if t is not None ]

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
            availableWells = [ r for r in product(bestRows,range(self.ncol)) if self.M[r[0]][r[1]] is None ]
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
