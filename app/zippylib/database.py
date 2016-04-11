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
from .primer import Primer, Locus, PrimerPair, Location

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
            cursor.execute('CREATE TABLE IF NOT EXISTS primer(name TEXT, seq TEXT, tag TEXT, tm REAL, gc REAL, vessel INT, well TEXT, dateadded TEXT, FOREIGN KEY(seq) REFERENCES target(seq), UNIQUE (name), UNIQUE (vessel, well));')
            cursor.execute('CREATE TABLE IF NOT EXISTS target(seq TEXT PRIMARY KEY, chrom TEXT, position INT, reverse BOOLEAN, FOREIGN KEY(seq) REFERENCES primer(seq));')
            cursor.execute('CREATE INDEX IF NOT EXISTS seq_index_in_target ON target(seq);')
            cursor.execute('CREATE TABLE IF NOT EXISTS pairs(pairid TEXT, uniqueid TEXT, left TEXT, right TEXT, chrom TEXT, start INT, end INT, dateadded TEXT, FOREIGN KEY(left) REFERENCES primer(name), FOREIGN KEY(right) REFERENCES primer(name), FOREIGN KEY(pairid, uniqueid) REFERENCES status(pairid, uniqueid) ON DELETE CASCADE, UNIQUE (pairid, uniqueid) ON CONFLICT REPLACE);')
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
                p.pairid, p.uniqueid, p.left, l.seq, p.right, r.seq, p.chrom, p.start, p.end, p.dateadded
                FROM pairs as p
                LEFT JOIN primer as l ON l.name = p.left
                LEFT JOIN primer as r ON r.name = p.right;''')
            rows = cursor.fetchall()
        finally:
            self.db.close()
        return "\n".join([ '{:<20} {:40} {:>20} {:<25} {:>20} {:<25} {:>8} {:>9d} {:>9d} {}'.format(*row) for row in rows ])

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
                    # delete all those pairs from pairs table
                    second_cursor.execute('''DELETE FROM pairs
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
            current_time = datetime.datetime.now()
            for p in primers:
                cursor = self.db.cursor()
                cursor.execute('''INSERT OR REPLACE INTO primer(name,seq,tag,tm,gc,dateadded) VALUES(?,?,?,?,?,?)''', \
                    (p.name, p.seq, p.tag, p.tm, p.gc, current_time))
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
                cursor.execute('''INSERT OR REPLACE INTO pairs(pairid,uniqueid,left,right,chrom,start,end,dateadded) VALUES(?,?,?,?,?,?,?,?)''', \
                    (p.name, p.uniqueid(), left.name, right.name, chrom, start, end, current_time))
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
            cursor.execute('''SELECT DISTINCT p.pairid, p.left, p.right, p.chrom, p.start, p.end, l.vessel, l.well, r.vessel, r.well,
                abs(p.start+((p.end-p.start)/2) - ?) as midpointdistance
                FROM pairs AS p
                LEFT JOIN primer as l ON p.left = l.name
                LEFT JOIN primer as r ON p.right = r.name
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
            primerPairs.append(PrimerPair([leftPrimer, rightPrimer], locations=[Location(row[6:8]),Location(row[8:10])]))
        return primerPairs  # ordered by midpoint distance

    def getLocation(self,loc):
        '''returns whats stored at location'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            cursor.execute('''SELECT DISTINCT name FROM primer
                WHERE vessel = ? AND well LIKE ?;''', (loc.vessel(),'%'+loc.well()+'%') )
            return [ x[0] for x in cursor.fetchall() ]
        finally:
            self.db.close()

    def addLocations(self, *locations):
        '''updates location for a batch of primers'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            # update
            try:
                cursor = self.db.cursor()
                cursor.executemany('''UPDATE OR IGNORE primer
                    SET vessel = ?, well = ? WHERE name = ?''', \
                    ((loc.vessel(), loc.well(), primerid) for primerid, loc in locations))
                self.db.commit()
            except:
                raise
        finally:
            self.db.close()
        return

    def storePrimer(self,primerid,loc):
        '''updates the location in which primer pairs are stored in the status table'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            # update
            try:
                cursor = self.db.cursor()
                cursor.execute('''UPDATE OR ABORT primer SET vessel = ?, well = ?
                    WHERE name = ?''', (loc.vessel(), loc.well(), primerid))
                self.db.commit()
            except sqlite3.IntegrityError:
                return False
            except:
                raise
            # check if updated
            try:
                cursor.execute('''SELECT DISTINCT vessel, well, pairid
                    FROM status WHERE name = ?;''', (primerid,) )
                rows = cursor.fetchall()
                assert len(rows)==1 and Location(rows[0][0],rows[0][1]) == loc
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
                if 'size' in kwargs.keys() and len(kwargs['size'])==2:
                    cursor.execute('''SELECT DISTINCT p.chrom, p.start, p.end, p.pairid
                        FROM pairs AS p
                        WHERE p.end - p.start >= ? AND p.end - p.start <= ?
                        ORDER BY p.chrom, p.start;''', (kwargs['size'][0],kwargs['size'][1]) )
                else:
                    cursor.execute('''SELECT DISTINCT p.chrom, p.start, p.end, p.pairid
                        FROM pairs AS p
                        ORDER BY p.chrom, p.start;''')
                rows = cursor.fetchall()
            finally:
                self.db.close()
            return rows, ('chrom','chromStart','chromEnd','name')  # rows and colnames
            # return [ '{}\t{}\t{}\t{}'.format(*row) for row in rows ]
        elif what=='ordersheet':
            try:
                self.db = sqlite3.connect(self.sqlite)
            except:
                raise
            else:
                cursor = self.db.cursor()
                # PAIRS
                cursor.execute('''SELECT DISTINCT
                    p.pairid AS pairname, l.name AS primername, l.seq AS sequence, l.tag as seqtag, 'fwd' AS direction
                    FROM pairs AS p
                    LEFT JOIN primer as l ON p.left = l.name
                    WHERE l.well IS NULL AND l.vessel IS NULL
                    UNION
                    SELECT DISTINCT
                    p.pairid AS pairname, r.name AS primername, r.seq AS sequence, r.tag AS seqtag, 'rev' AS direction
                    FROM pairs AS p
                    LEFT JOIN primer as r ON p.right = r.name
                    WHERE r.well IS NULL AND r.vessel IS NULL
                    ORDER BY pairname, direction;''')
                rows = cursor.fetchall()
            finally:
                self.db.close()
            # define columns
            columns = ['pairname','primername','sequence','seqtag','direction']
            # add tags and extra columns
            if 'extra' in kwargs.keys() and kwargs['extra']:
                columns += [ c[0] for c in extra ]
                for i in range(len(rows)):
                    rows[i] = list(rows[i]) + [ c[1] for c in kwargs['extra'] ]
            # add sequence tag
            if 'tags' in kwargs.keys() and kwargs['tags']:
                for row in rows:
                    # get correct tag
                    try:
                        if row[4] == 'fwd':
                            prepend = tags[row[3]]['tags'][0]
                        elif row[4] == 'rev':
                            prepend = tags[row[3]]['tags'][1]
                        else:
                            prepend = tags[row[3]]
                    except:
                        row[2] = row[3] + '-' + row[2]  # prepend tag sequence
                    else:
                        row[2] = prepend + row[2]  # prepend tag sequence
            # return rows (list of list) and column names (headers)
            return rows, columns
        elif what=='locations':
            # dump locations (all possible)
            try:
                self.db = sqlite3.connect(self.sqlite)
            except:
                raise
            else:
                cursor = self.db.cursor()
                # PAIRS
                cursor.execute('''SELECT DISTINCT pp.pairid, p.name, p.vessel, p.well
                        FROM pairs AS pp, primer as p
                        WHERE pp.left = p.name
                    UNION
                    SELECT DISTINCT pp.pairid, p.name, p.vessel, p.well
                        FROM pairs AS pp, primer as p
                        WHERE pp.right = p.name;''')
                rows = cursor.fetchall()
            finally:
                self.db.close()
            return rows, ['pair','primer','vessel','well']
