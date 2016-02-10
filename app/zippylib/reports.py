#!/usr/bin/env python

__doc__=="""
####################
# Report generator #
####################
"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Prototype"

import sys
import time
from math import ceil
from copy import deepcopy
from random import shuffle
import itertools
from collections import Counter
from . import PlateError, char_range

from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm, inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.rl_config import defaultPageSize


class Report(object):
    def __init__(self,fi,title='This is the title'):
        # get document
        self.doc = SimpleDocTemplate(fi,pagesize=A4,
                        rightMargin=2.54*cm,leftMargin=2.54*cm,
                        topMargin=2*cm,bottomMargin=cm)
        # style sheets
        self.styles = getSampleStyleSheet()
        self.styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
        self.styles.add(ParagraphStyle(name='Warning', backColor = '#FFFF00', borderColor = "#000000"))
        self.styles.add(ParagraphStyle(name='tiny', fontSize=4, leading=4))
        self.styles.add(ParagraphStyle(name='small', fontSize=6, leading=6))
        # flowable elements
        self.elements = []
        # Header
        logo = Image("viapath.png",width=1.32*inch,height=0.7*inch)
        titl = Paragraph('%s' % title, self.styles["Heading1"])
        date = Paragraph('<font size=12>Generated on %s</font>' % time.ctime(), self.styles["Normal"])
        self.elements.append(Table([[logo,titl],['',date]],
            colWidths=[1.5*inch,4.5*inch],
            style=[ ('SPAN',(0,0),(0,1)), ('VALIGN',(0,0),(-1,-1),'BOTTOM'),
                ('VALIGN',(1,1),(-1,-1),'TOP'), ('ALIGN',(1,0),(-1,-1),'LEFT') ]))
        self.elements.append(Spacer(1, 12))
        # Warning
        # self.elements.append(Paragraph('<font size=14>Always double check everything!</font>', self.styles["Warning"]))
        # self.elements.append(Spacer(1, 24))

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
            ('BACKGROUND',(1,0),(-1,0),colors.snow)
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
            t = Table(data, colWidths=cm, rowHeights=0.5*cm)
            t.setStyle(PLATE_STYLE)
            self.elements.append(t)
            self.elements.append(Spacer(1, 12))

    def samplePrimerLists(self,s,p):
        TABLE_STYLE = TableStyle([
            ('FONTSIZE',(0,1),(-1,-1),8),
            ('FONTSIZE',(0,0),(-1,0),10),
            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
            ('ALIGN',(0,0),(-1,-1),'LEFT'),
            ('INNERGRID', (0,1), (1,len(s)), 0.25, colors.black),
            ('LINEABOVE', (0,1),(1,1),1,colors.black),
            ('BOX', (0,0), (1,len(s)), 1, colors.black),
            ('INNERGRID', (3,1), (4,len(p)), 0.25, colors.black),
            ('LINEABOVE', (3,1),(4,1),1,colors.black),
            ('BOX', (3,0), (4,len(p)), 1, colors.black),
            ])

        # self.elements.append(Paragraph('Sample - Primers', self.styles["Heading2"]))
        self.elements.append(Spacer(1, 2))
        data = [[str(len(s)),'Samples','',str(len(p)),'Primer Pairs']]
        for i in range(max(len(s),len(p))):
            data.append([ '', s[i] if i<len(s) else '', '', '', p[i] if i<len(p) else '' ])
        t = Table(data, colWidths=[0.6*cm,7*cm,0.3*cm,0.6*cm,7*cm], rowHeights=0.6*cm)
        t.setStyle(TABLE_STYLE)
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
            ('BACKGROUND', (0,1), (0,-1), colors.bisque)
            ])

        self.elements.append(Paragraph('Checks', self.styles["Heading4"]))
        self.elements.append(Spacer(1, 2))
        data = [[ 'Task', 'Date', 'Checker']]
        for i in range(len(titles)):
            data.append([ titles[i], '', '' ])
        t = Table(data, colWidths=[7.5*cm,4*cm,4*cm], rowHeights=0.6*cm)
        t.setStyle(TABLE_STYLE)
        self.elements.append(t)
        self.elements.append(Spacer(1, 12))


    def build(self):
        self.doc.build(self.elements)


'''PCR test [PrimerPair X Sample] (multiple variants)'''
class Test(object):
    def __init__(self,sample,primerpair):
        self.sample = sample
        self.primerpair = primerpair
        self.control = False
        self.variants = []  # variants to test for (name, offset)

    def __str__(self):
        return '<'+'|'.join([self.sample,self.primerpair,str(self.control),'Variants:'+str(len(self.variants))])+'>'

    def __repr__(self):
        return self.sample[:2]+("-" if self.control else "+")+self.primerpair[:2]

    def __hash__(self):
        #hash((self.sample,self.primerpair,self.control))
        return hash(self.sample) ^ hash(self.primerpair) ^ hash(self.control)

    def __eq__(self,other):
        return hash(self) == hash(other)

'''Worksheet data'''
class Worksheet(list):
    def __init__(self,elements,fieldNames=['sample', 'variant', 'primerpair', \
        'fwd', 'fwdTm', 'fwdGC', 'rev', 'revTm', 'revGC', 'chrom', \
        'chromStart', 'chromEnd']):
        list.__init__(self, [])
        self.plates = []
        for i, e in enumerate(elements):
            data = dict(zip(fieldNames,e))
            t = Test(sample=data['sample'],primerpair=data['primerpair'])
            if t in self:
                for test in self[::-1]:  # go backwards as variants are sorted by sample and coordinate
                    if test == t:
                        test.variants.append(data['variant'])
                        break
            else:
                t.variants.append(data['variant'])
                self.append(t)

    def __str__(self):
        return "<Worksheet (%d tests, %d variants)> " % (len(self),sum([len(e.variants) for e in self]))

    def count(self,attr,x):
        return len([ e for e in self if getattr(e,attr)==x ])

    '''print robot csv, (DestinationPlate,DestinationWell,SampleID,PrimerID)'''
    def robotCsv(self,fi,sep=','):
        with open(fi,'w') as fh:
            print >> fh, sep.join(['DestinationPlate','DestinationWell','PrimerID','SampleID','Control','Confirmations'])
            for n, p in enumerate(self.plates):
                P = 'Plate'+str(n+1)
                for i,row in enumerate(p.M):
                    R = chr(ord('A')+i)
                    for j,cell in enumerate(row):
                        C = str(j+1)
                        if cell:
                            print >> fh, sep.join([P,R+C,cell.primerpair,cell.sample,str(cell.control),str(len(cell.variants))])

    '''assign tests to plate (smart work better only for big test sets)'''
    def fillPlates(self,size=[8,12],randomize=True):
        if randomize:
            shuffle(self)  # shuffle filling
        else:
            self = sorted(self, key=lambda x: x.sample)  # sort by sample
        # fill plates
        self.plates.append(Plate(*size))  # start plate
        for t in self:
            if self.plates[-1].isfull():
                self.plates.append(Plate(*size))
            # add sample
            self.plates[-1].add(t,'sample')
        # compress plates
        for p in self.plates:
            p.compress('sample')
        return

    '''add control samples'''
    def addControls(self):
        controlsamples = []
        for e in self:
            if not e.control:
                controlsamples.append(deepcopy(e))
                controlsamples[-1].control = True
        self += controlsamples
        return

    # print worksheet
    def workSheet(self,fi):
        r = Report(fi,'Batch PCR')
        # add plates
        samples, primers, plates = [], [], []
        for p in self.plates:
            s, p, m = p.platemap()
            samples += s
            primers += p
            plates.append(m)
        # sample and primer list
        r.samplePrimerLists(list(set(samples)),list(set(primers)))
        # add checkboxes
        r.checkBoxes(['Plate orientation'])
        # plate layout
        r.plateLayouts(plates)
        # build pdf
        r.build()


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
        s, p, c = set(), set(), set()
        for r in range(len(self.M)):
            for c in range(len(self.M[r])):
                if self.M[r][c] is not None:
                    s.add(self.M[r][c].sample)
                    p.add(self.M[r][c].primerpair)
        # generate plate map
        trunc = lambda x: str(x)[1:].split('|')[0]+'\n'+str(x)[1:].split('|')[1].split(':')[0] if x is not None else "EMPTY"
        pm = [ map(trunc,r) for r in self.M ]
        # return values
        return sorted(list(s)), sorted(list(p)), pm

    def _bestRows(self,x,rvalue):
        rowcounts = { i: len([ e for e in r if e is not None and getattr(e,rvalue) == getattr(x,rvalue) ]) - \
            len([ e for e in r if e is not None and getattr(e,rvalue) != getattr(x,rvalue) ]) \
            for i, r in enumerate(self.M) }
        return [ x[0] for x in sorted(rowcounts.items(), key=lambda x: x[1], reverse=True) ]

    def isfull(self):
        return True if len(self)>=len(self.M)*len(self.M[0]) else False

    # fills with row priority only (corresponds to sample)
    def add(self,t,rvalue):
        # get best rows
        try:
            bestRows = self._bestRows(t,rvalue)
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
        maxColumn = int(len(self)/float(self.nrow))+1
        # start with longest row (and only rows wich exceed maxColumns)
        rowOrder = { i: r.count(None) for i,r in enumerate(self.M) if r.count(None) < self.ncol-maxColumn }
        longRows = [ i for i,r in sorted(rowOrder.items(), key=lambda x: x[1], reverse=False) ]
        for r in longRows:
            # find best row (sum of empty and own)
            for c in range(maxColumn,self.ncol):
                if self.M[r][c] is not None:
                    bestRows = [ br for br in self._bestRows(self.M[r][c],rvalue) \
                        if br!=r and self.M[br].count(None) > self.ncol-maxColumn ]
                    # find leftmost free position
                    for bc in range(self.ncol):
                        if self.M[bestRows[0]][bc] is None:
                            self.M[r][c], self.M[bestRows[0]][bc] = self.M[bestRows[0]][bc], self.M[r][c]  # switch position
                            break
        return
