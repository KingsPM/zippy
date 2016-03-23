#!/usr/bin/env python

__doc__=="""SQLITE Database API"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "1.2"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys, os, re, ast
import datetime
import json
import hashlib
import sqlite3
import fnmatch
import copy
from collections import defaultdict
from . import flatten
from .primer import Primer, Locus, PrimerPair

class PrimerDB(object):
    def __init__(self, database, user='unkown'):
        # open database and get a cursor
        self.sqlite = database
        self.db = sqlite3.connect(self.sqlite)
        self.user = user
        # create file table if not exists
        cursor = self.db.cursor()
        try:
            # TABLE
            cursor.execute('PRAGMA foreign_keys = ON')
            cursor.execute('CREATE TABLE IF NOT EXISTS primer(seq TEXT PRIMARY KEY, tm REAL, gc REAL, FOREIGN KEY(seq) REFERENCES target(seq));')
            cursor.execute('CREATE TABLE IF NOT EXISTS target(seq TEXT, chrom TEXT, position INT, reverse BOOLEAN, FOREIGN KEY(seq) REFERENCES primer(seq));')
            cursor.execute('CREATE TABLE IF NOT EXISTS pairs(pairid TEXT, uniqueid TEXT, left TEXT, right TEXT, chrom TEXT, start INT, end INT, FOREIGN KEY(left) REFERENCES primer(seq), FOREIGN KEY(right) REFERENCES primer(seq), FOREIGN KEY(pairid, uniqueid) REFERENCES status(pairid, uniqueid) ON DELETE CASCADE, UNIQUE (pairid, uniqueid) ON CONFLICT REPLACE);')
            cursor.execute('CREATE TABLE IF NOT EXISTS status(pairid TEXT NOT NULL, uniqueid TEXT NOT NULL, dateadded TEXT, vessel INT, well TEXT, PRIMARY KEY (pairid, uniqueid) ON CONFLICT REPLACE, UNIQUE (vessel, well));')
            cursor.execute('CREATE INDEX IF NOT EXISTS seq_index_in_target ON target(seq);')
            cursor.execute('CREATE TABLE IF NOT EXISTS blacklist(uniqueid TEXT PRIMARY KEY, blacklistdate TEXT );')
            self.db.commit()
        except:
            print >> sys.stderr, self.sqlite
            print >> sys.stderr, self.db
            print >> sys.stderr, cursor
            raise
        finally:
            self.db.close()
        return

    def __str__(self):
        return '<ZippyDB at %s>' % self.sqlite

    def __repr__(self):
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            cursor = self.db.cursor()
            cursor.execute('''SELECT DISTINCT
                p.pairid, p.uniqueid, p.left, p.right, p.chrom, p.start, p.end, s.dateadded
                FROM pairs as p, status as s
                WHERE p.pairid = s.pairid AND p.uniqueid = s.uniqueid;''')
            rows = cursor.fetchall()
        finally:
            self.db.close()
        return "\n".join([ '{:<20} {:40} {:<25} {:<25} {:<8} {:>9d} {:>9d} {}'.format(*row) for row in rows ])

    '''show/update blacklist'''
    def blacklist(self,add=None):
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            if add:
                # get uniqueid from status table for pairid
                blacklisttime = datetime.datetime.now()
                cursor = self.db.cursor()
                cursor.execute('PRAGMA foreign_keys = ON')
                cursor.execute('''SELECT DISTINCT s.uniqueid
                    FROM status AS s WHERE s.pairid = ?;''', (add,))
                bl_uniqueid = [ row[0] for row in cursor.fetchall() ]
                # get list of pairs from pairs table with uniqueid
                second_cursor = self.db.cursor()
                second_cursor.execute('PRAGMA foreign_keys = ON')
                pairlist = []
                for uid in bl_uniqueid:
                    # add uniqueid to blacklist
                    second_cursor.execute('''INSERT INTO blacklist(uniqueid, blacklistdate) VALUES(?,?);''', \
                    (uid, blacklisttime))
                    # get list of pairs to be deleted
                    second_cursor.execute('''SELECT DISTINCT p.pairid
                    FROM pairs AS p
                    WHERE p.uniqueid = ?;''', (uid,))
                    pairlist += [ x[0] for x in second_cursor.fetchall() ]
                    # delete all those pairs from status table
                    second_cursor.execute('''DELETE FROM status
                    WHERE uniqueid = ?;''', (uid,))
                    self.db.commit()
                return pairlist
            else: #return list of uniqueids from blacklist
                cursor = self.db.cursor()
                cursor.execute('''SELECT DISTINCT uniqueid FROM blacklist;''')
                rows = cursor.fetchall()
                return [ row[0] for row in rows ]
        finally:
            self.db.close()

    '''adds list of primers to database'''
    def addPrimer(self, *primers):
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            for p in primers:
                cursor = self.db.cursor()
                cursor.execute('''INSERT OR REPLACE INTO primer(seq,tm,gc) VALUES(?,?,?)''', \
                    (p.seq, p.tm, p.gc))
                for l in p.loci:
                    cursor.execute('''INSERT OR REPLACE INTO target(seq,chrom,position,reverse) VALUES(?,?,?,?)''', \
                        (p.seq, l.chrom, l.offset, l.reverse))
            self.db.commit()
        finally:
            self.db.close()
        return

    def addPair(self, *pairs):
        '''adds primer pairs (and individual primers)'''
        # add primers
        flat = []
        for p in pairs:
            flat.append(p[0])
            flat.append(p[1])

        self.addPrimer(*flat)
        # add pairs
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            current_time = datetime.datetime.now()
            cursor = self.db.cursor()
            # add pairs
            for p in pairs:
                left, right = p[0],p[1]
                # find common substring in name for automatic naming
                chrom = left.targetposition.chrom
                start = left.targetposition.offset
                end = right.targetposition.offset+right.targetposition.length
                cursor.execute('''INSERT OR REPLACE INTO pairs(pairid,uniqueid,left,right,chrom,start,end) VALUES(?,?,?,?,?,?,?)''', \
                    (p.name(), p.uniqueid(), left.seq, right.seq, chrom, start, end))
                cursor.execute('''INSERT OR IGNORE INTO status(pairid,uniqueid,dateadded) VALUES(?,?,?)''', \
                    (p.name(), p.uniqueid(), current_time))
            self.db.commit()
        finally:
            self.db.close()
        return

    def query(self, variant):
        '''returns suitable primer pairs for the specified interval'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            cursor = self.db.cursor()
            cursor.execute('''SELECT DISTINCT p.pairid, p.left, p.right, p.chrom, p.start, p.end, s.vessel, s.well,
                abs(p.start+((p.end-p.start)/2) - ?) as midpointdistance
                FROM pairs AS p, status AS s
                WHERE p.pairid = s.pairid AND p.uniqueid = s.uniqueid
                AND p.chrom = ?
                AND p.start + length(p.left) <= ?
                AND p.end - length(p.right) >= ?
                ORDER BY midpointdistance;''', \
                (int(variant.chromStart+int(variant.chromEnd-variant.chromStart)/2.0), variant.chrom, variant.chromStart, variant.chromEnd))
            rows = cursor.fetchall()
        finally:
            self.db.close()
        # return primer pairs that would match
        primerPairs = []
        for row in rows:
            name = row[0]
            leftSeq = row[1]
            rightSeq = row[2]
            leftTargetposition = Locus(row[3], row[4], len(row[1]), False)
            rightTargetposition = Locus(row[3], row[5]-len(row[2]), len(row[2]), True)
            leftPrimer = Primer(name+'_left', leftSeq, leftTargetposition)
            rightPrimer = Primer(name+'_right', rightSeq, rightTargetposition)
            leftPrimer.calcProperties()
            rightPrimer.calcProperties()
            primerPairs.append(PrimerPair([leftPrimer, rightPrimer], location=row[6:8]))
        return primerPairs  # ordered by midpoint distance

    def storePrimer(self,pairid,vessel,well):
        '''updates the location in which primer pairs are stored in the status table'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            # update
            try:
                cursor = self.db.cursor()
                cursor.execute('''UPDATE OR ABORT status SET vessel = ?, well = ?
                    WHERE pairid = ?''', (vessel, well, pairid))
                self.db.commit()
            except sqlite3.IntegrityError:
                return False
            except:
                raise
            # check if updated
            try:
                cursor.execute('''SELECT DISTINCT vessel, well, pairid
                    FROM status WHERE pairid = ?;''', (pairid,) )
                rows = cursor.fetchall()
                assert len(rows)==1 and rows[0][0] == vessel and rows[0][1] == well
            except AssertionError:
                return False
            except:
                raise
        finally:
            self.db.close()
        return True

    def dump(self,what,**kwargs):
        if what=='amplicons':
            # dump amplicons (all possible)
            try:
                self.db = sqlite3.connect(self.sqlite)
            except:
                raise
            else:
                cursor = self.db.cursor()
                # FWD strand
                cursor.execute('''SELECT DISTINCT p.chrom, p.start, p.end, p.pairid
                    FROM pairs AS p
                    ORDER BY p.chrom, p.start;''')
                rows = cursor.fetchall()
            finally:
                self.db.close()
            return rows, ('chrom','chromStart','chromEnd','name')  # rows and colnames
            # return [ '{}\t{}\t{}\t{}'.format(*row) for row in rows ]
        elif what=='ordersheet':
            # dump pending orders
            try:
                self.db = sqlite3.connect(self.sqlite)
            except:
                raise
            else:
                cursor = self.db.cursor()
                # PAIRS
                cursor.execute('''SELECT DISTINCT p.pairid, p.left, p.right
                    FROM pairs AS p, status AS s
                    WHERE p.pairid = s.pairid AND p.uniqueid = s.uniqueid
                    AND s.well IS NULL AND s.vessel IS NULL
                    ORDER BY s.dateadded, p.pairid;''')
                rows = cursor.fetchall()
            finally:
                self.db.close()
            # define columns
            columns = ['name','forward','reverse']
            # add tags and extra columns
            if "extracolumns" in kwargs.keys():
                extracolumns = kwargs['extracolumns'].keys()
                columns += extracolumns
                for i in range(len(rows)):
                    rows[i] = list(rows[i]) + [ kwargs['extracolumns'][c] for c in extracolumns ]
            # add sequence tag
            if "sequencetags" in kwargs.keys():
                for row in rows:
                    row[1] = kwargs['sequencetags']['left'] + row[1]
                    row[2] = kwargs['sequencetags']['right'] + row[2]
            # expand fwd reverse
            expandedrows = []
            columns = tuple(columns[:1] + [ 'sequence' ] + columns[3:])
            for row in rows:
                expandedrows.append([row[0]+'_fwd', row[1]] + row[3:])
                expandedrows.append([row[0]+'_rev', row[2]] + row[3:])
            return expandedrows, columns
        elif what=='locations':
            # dump locations (all possible)
            try:
                self.db = sqlite3.connect(self.sqlite)
            except:
                raise
            else:
                cursor = self.db.cursor()
                # PAIRS
                cursor.execute('''SELECT DISTINCT p.pairid, s.vessel, s.well
                    FROM pairs AS p, status AS s
                    WHERE p.pairid = s.pairid AND p.uniqueid = s.uniqueid;''')
                rows = cursor.fetchall()
            finally:
                self.db.close()
            return rows, ['name','vesel','well']
