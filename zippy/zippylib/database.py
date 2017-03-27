#!/usr/bin/env python

__doc__=="""SQLITE Database API"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "2.3.4"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys, os, re, ast
import datetime
import json
import hashlib
import sqlite3
import fnmatch
import primer3
from copy import deepcopy
from collections import defaultdict
from . import flatten
from .primer import Primer, Locus, PrimerPair, Location, parsePrimerName

# changes conflicting name
def changeConflictingName(n):
    f = n.split('_')
    # find suffix to increment
    if len(f)==4:  # try to increase suffix
        try:
            f[2] = str(int(f[2])+1)
            newname = '_'.join(f)
        except:
            raise Exception('PrimerNameChangeError')
        return newname
    elif len(f)==3:  # add number suffix after exon
        return '_'.join(f[:2]+[str(1)]+f[2:])
    raise Exception('PrimerNameChangeError')

# Primer Database
class PrimerDB(object):
    def __init__(self, database, dump=None):
        # open database and get a cursor
        self.sqlite = database
        self.db = sqlite3.connect(self.sqlite)
        self.dump = dump  # Primer BED file created by destructor
        # create file table if not exists
        cursor = self.db.cursor()
        try:
            # TABLE
            cursor.execute('''PRAGMA foreign_keys = ON''')
            # names unique, reinsertion if identical primers ignored
            cursor.execute('''CREATE TABLE IF NOT EXISTS primer(
                name TEXT, seq TEXT, tag TEXT, tm REAL, gc REAL, vessel INT,
                well TEXT, dateadded TEXT,
                FOREIGN KEY(seq) REFERENCES target(seq),
                CONSTRAINT uniquenames PRIMARY KEY (name) ON CONFLICT ABORT,
                UNIQUE (name,seq,tag) ON CONFLICT IGNORE,
                UNIQUE (vessel, well));''')
            cursor.execute('''CREATE TABLE IF NOT EXISTS pairs(
                pairid TEXT PRIMARY KEY, uniqueid TEXT, left TEXT, right TEXT,
                chrom TEXT, start INT, end INT, dateadded TEXT,
                FOREIGN KEY(left) REFERENCES primer(name) ON UPDATE CASCADE,
                FOREIGN KEY(right) REFERENCES primer(name) ON UPDATE CASCADE,
                UNIQUE (pairid, uniqueid) ON CONFLICT IGNORE);''')
            cursor.execute('''CREATE TABLE IF NOT EXISTS target(
                seq TEXT, chrom TEXT, position INT, reverse BOOLEAN, tm REAL,
                UNIQUE (seq,chrom,position,reverse),
                FOREIGN KEY(seq) REFERENCES primer(seq) ON DELETE CASCADE);''')
            cursor.execute('''CREATE TABLE IF NOT EXISTS blacklist(
                uniqueid TEXT PRIMARY KEY, blacklistdate TEXT);''')
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

    def writeAmpliconDump(self):
        ## dump amplicons to bed file
        if self.dump:
            try:
                self.db = sqlite3.connect(self.sqlite)
            except:
                raise
            else:
                # get amplicons
                cursor = self.db.cursor()
                cursor.execute('''SELECT DISTINCT p.chrom, p.start, p.end, p.pairid
                    FROM pairs AS p
                    ORDER BY p.chrom, p.start;''')
                rows = cursor.fetchall()
                # write bed file
                try:
                    with open(self.dump,'w') as fh:
                        for row in rows:
                            print >> fh, '\t'.join(map(str,row))
                except IOError:
                    print >> sys.stderr, "cannot write to %s" % self.dump
                    pass  # fail silently (eg if data cannot be written)
                except:
                    raise
            finally:
                self.db.close()

    def removeOrphans(self):
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            cursor = self.db.cursor()
            cursor.execute('''
            SELECT p.name FROM primer AS p WHERE NOT EXISTS (
            SELECT left FROM pairs AS pp WHERE pp.left = p.name
            UNION
            SELECT right FROM pairs AS pp WHERE pp.right = p.name);''')
            orphans = cursor.fetchall()
            cursor.executemany('''DELETE FROM primer
                WHERE name = ?''', (orphans))
            self.db.commit()
            return [ x[0] for x in orphans ]
        finally:
            self.db.close()

    '''show/update blacklist'''
    def blacklist(self,add=None,justdelete=False):
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
                cursor.execute('''SELECT DISTINCT p.uniqueid
                    FROM pairs AS p WHERE p.pairid = ?;''', (add,))
                bl_uniqueid = [ row[0] for row in cursor.fetchall() ]
                # get list of pairs from pairs table with uniqueid
                second_cursor = self.db.cursor()
                second_cursor.execute('PRAGMA foreign_keys = OFF')
                pairlist = []
                for uid in bl_uniqueid:
                    # add uniqueid to blacklist
                    if not justdelete:
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
            self.removeOrphans()
            self.writeAmpliconDump()


    '''adds list of primers to database and automatically renames'''
    def addPrimer(self, *primers):
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            current_time = datetime.datetime.now()
            cursor = self.db.cursor()
            for p in primers:
                # store primers and modify names if necessary
                while True:
                    originalName = deepcopy(p.name)
                    try:
                        cursor.execute('''INSERT INTO primer(name,seq,tag,tm,gc,dateadded) VALUES(?,?,?,?,?,?)''', \
                            (p.name, p.seq, p.tag, p.tm, p.gc, current_time))
                    except sqlite3.IntegrityError:
                        try:
                            p.name = changeConflictingName(p.name)
                        except Exception as e:
                            raise e
                    except:
                        raise
                    else:
                        if originalName != p.name:
                            print >> sys.stderr, "WARNING: renamed primer {} -> {} in database".format(originalName, p.name)
                        break  # sucessfully stored
                # store mapping loci
                for l in p.loci:
                    cursor.execute('''INSERT OR IGNORE INTO target(seq,chrom,position,reverse,tm) VALUES(?,?,?,?,?)''', \
                        (p.seq, l.chrom, l.offset, l.reverse, l.tm))
            self.db.commit()
        finally:
            self.db.close()
            self.writeAmpliconDump()
        return

    def addPair(self, *pairs):
        '''adds primer pairs (and individual primers)'''
        # add primers (and rename if necessary)
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
                p.fixName()  # changes name if there is a longer common name (catches primer renaming)
                # find common substring in name for automatic naming
                chrom = p[0].targetposition.chrom
                start = p[0].targetposition.offset
                end = p[1].targetposition.offset+p[1].targetposition.length
                cursor.execute('''INSERT OR IGNORE INTO pairs(pairid,uniqueid,left,right,chrom,start,end,dateadded) VALUES(?,?,?,?,?,?,?,?)''', \
                    (p.name, p.uniqueid(), p[0].name, p[1].name, chrom, start, end, current_time))
            self.db.commit()
        finally:
            self.db.close()
            self.writeAmpliconDump()
        return

    '''query for interval or name'''
    def query(self, query,opendb=None):
        '''returns suitable primer pairs for the specified interval'''
        try:
            self.db = opendb if opendb else sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            cursor = self.db.cursor()
            datematch = re.compile("([0-9\s-]+)$")
            if datematch.match(str(query)): # query date
                subSearchName = '%'+query+'%'
                cursor.execute('''SELECT DISTINCT p.pairid, l.tag, r.tag, l.seq, r.seq, p.left, p.right,
                    p.chrom, p.start, p.end, l.vessel, l.well, r.vessel, r.well, 0
                    FROM pairs AS p
                    LEFT JOIN primer as l ON p.left = l.name
                    LEFT JOIN primer as r ON p.right = r.name
                    where p.dateadded LIKE ?
                    ORDER BY p.pairid;''', \
                    (subSearchName,))
            elif type(query) in [str,unicode]:  # use primerpair name
                subSearchName = '%'+query+'%'
                cursor.execute('''SELECT DISTINCT p.pairid, l.tag, r.tag, l.seq, r.seq, p.left, p.right,
                    p.chrom, p.start, p.end, l.vessel, l.well, r.vessel, r.well, 0
                    FROM pairs AS p
                    LEFT JOIN primer as l ON p.left = l.name
                    LEFT JOIN primer as r ON p.right = r.name
                    WHERE p.pairid LIKE ?
                    ORDER BY p.pairid;''', \
                    (subSearchName,))
            else:  # is interval
                cursor.execute('''SELECT DISTINCT p.pairid, l.tag, r.tag, l.seq, r.seq, p.left, p.right,
                    p.chrom, p.start, p.end, l.vessel, l.well, r.vessel, r.well,
                    abs(p.start+((p.end-p.start)/2) - ?) as midpointdistance
                    FROM pairs AS p
                    LEFT JOIN primer as l ON p.left = l.name
                    LEFT JOIN primer as r ON p.right = r.name
                    WHERE p.chrom = ?
                    AND p.start + length(l.seq) <= ?
                    AND p.end - length(r.seq) >= ?
                    ORDER BY midpointdistance;''', \
                    (int(query.chromStart+int(query.chromEnd-query.chromStart)/2.0), query.chrom, query.chromStart, query.chromEnd))
            rows = cursor.fetchall()
        finally:
            if not opendb:
                self.db.close()
        # return primer pairs that would match
        primerPairs = []
        for row in rows:
            # build targets
            leftTargetposition = Locus(row[7], row[8], len(row[3]), False, primer3.calcTm(str(row[3])))
            rightTargetposition = Locus(row[7], row[9]-len(row[4]), len(row[4]), True, primer3.calcTm(str(row[4])))
            # build storage locations (if available)
            leftLocation = Location(*row[10:12]) if all(row[10:12]) else None
            rightLocation = Location(*row[12:14]) if all(row[12:14]) else None
            # Build primers
            leftPrimer = Primer(row[5], row[3], targetposition=leftTargetposition, tag=row[1], location=leftLocation)
            rightPrimer = Primer(row[6], row[4], targetposition=rightTargetposition, tag=row[2], location=rightLocation)
            # get reverse status (from name)
            orientations = [ x[1] for x in map(parsePrimerName,row[5:7]) ]
            if not any(orientations) or len(set(orientations))==1:
                print >> sys.stderr, '\rWARNING: {} orientation is ambiguous ({},{}){}\r'.format(row[0],\
                    '???' if orientations[0]==0 else 'rev' if orientations[0]<0 else 'fwd', \
                    '???' if orientations[0]==0 else 'rev' if orientations[1]<0 else 'fwd'," "*20)
                reverse = False
            elif orientations[0]>0 or orientations[1]<0:
                reverse = False
            elif orientations[1]>0 or orientations[0]<0:
                reverse = True
            else:
                raise Exception('PrimerPairStrandError')
            # Build pair
            primerPairs.append(PrimerPair([leftPrimer, rightPrimer],name=row[0],reverse=reverse))
        return primerPairs  # ordered by midpoint distance

    def getLocation(self,loc):
        '''returns whats stored at location'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            cursor = self.db.cursor()
            cursor.execute('''SELECT DISTINCT name FROM primer
                WHERE vessel = ? AND well LIKE ?;''', (loc.vessel(),'%'+loc.well()+'%') )
            return [ x[0] for x in cursor.fetchall() ]
        finally:
            self.db.close()

    def getRedundantPrimers(self):
        '''returns redundant primer (same tag and sequence)'''
        redundant = defaultdict(list)
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            cursor = self.db.cursor()
            cursor.execute('''SELECT p1.seq, p1.tag, p1.name
            FROM primer AS p1 LEFT JOIN primer AS p2 ON p1.seq = p2.seq AND p1.tag = p2.tag
            GROUP BY p1.name,p1.tag HAVING COUNT(p1.seq) > 1;''')
            for l in cursor.fetchall():
                redundant[(l[0],l[1])].append(l[2])
            # return list of list
            return [ [k[0], k[1], ','.join(v)] for k,v in redundant.items() ], ['seq','tag','synonyms']
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

    def storePrimer(self,primerid,loc,force=False):
        '''updates the location in which primers are stored'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            # reset storage location
            if force:
                try:
                    cursor = self.db.cursor()
                    cursor.execute('''UPDATE OR IGNORE primer SET vessel = NULL, well = NULL
                        WHERE vessel = ? AND well = ?''', (loc.vessel(), loc.well()))
                    self.db.commit()
                except:
                    raise
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
                cursor.execute('''SELECT DISTINCT vessel, well
                    FROM primer WHERE name = ?;''', (primerid,) )
                rows = cursor.fetchall()
                assert len(rows)==1 and Location(str(rows[0][0]),str(rows[0][1])) == loc
            except AssertionError:
                return False
            except:
                raise
        finally:
            self.db.close()
        return True

    def updateName(self,primerName,newName):
        '''changes the name of a primer stored in the database'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            # update primer name in primer and pairs table
            try:
                cursor = self.db.cursor()
                cursor.execute('''UPDATE OR ABORT pairs SET left = ?
                    WHERE left = ?;''', (newName, primerName))
                cursor.execute('''UPDATE OR ABORT pairs SET right = ?
                    WHERE right = ?;''', (newName, primerName))
                cursor.execute('''UPDATE OR ABORT primer SET name = ?
                    WHERE name = ?;''', (newName, primerName))
                self.db.commit()
            except sqlite3.IntegrityError:
                return False
            except:
                raise
            # check if updated
            try:
                cursor.execute('''SELECT DISTINCT name
                    FROM primer WHERE name = ?;''', (newName,) )
                rows = cursor.fetchall()
                assert rows[0][0] == newName
                return True
            except AssertionError:
                return False
            except:
                raise
        finally:
            self.db.close()
            self.writeAmpliconDump()

    def updatePairName(self,pairName,newName):
        '''changes the name of a primer stored in the database'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            # update primer name in primer and pairs table
            try:
                cursor = self.db.cursor()
                cursor.execute('''UPDATE OR ABORT pairs SET pairid = ?
                    WHERE pairid = ?;''', (newName, pairName))
                self.db.commit()
            except sqlite3.IntegrityError:
                return False
            except:
                raise
            # check if updated
            try:
                cursor.execute('''SELECT DISTINCT pairid
                    FROM pairs WHERE pairid = ?;''', (newName,) )
                rows = cursor.fetchall()
                assert rows[0][0] == newName
                return True
            except AssertionError:
                return False
            except:
                raise
        finally:
            self.db.close()
            self.writeAmpliconDump()

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
            if 'extracolumns' in kwargs.keys() and kwargs['extracolumns']:
                columns += [ c[0] for c in kwargs['extracolumns'] ]
                for i in range(len(rows)):
                    rows[i] = list(rows[i]) + [ c[1] for c in kwargs['extracolumns'] ]
            # add sequence tag
            if 'sequencetags' in kwargs.keys() and kwargs['sequencetags']:
                for row in rows:
                    # get correct tag
                    try:
                        p = parsePrimerName(row[1])
                        if p[1] > 0:
                            prepend = kwargs['sequencetags'][row[3]]['tags'][0]
                        elif p[1] < 0:
                            prepend = kwargs['sequencetags'][row[3]]['tags'][1]
                        else:
                            raise Exception('PrimerNameParseError')
                        print >> sys.stderr, prepend
                    except AssertionError:
                        raise
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
        elif what=='table':
            # dump table with pairs primers and locations (which can be reimported)
            try:
                self.db = sqlite3.connect(self.sqlite)
            except:
                raise
            else:
                cursor = self.db.cursor()
                cursor.execute('''SELECT * FROM (
                    SELECT p.name, pp.pairid, p.tag, p.seq, p.vessel, p.well
                    FROM pairs AS pp, primer AS p
                    WHERE pp.left = p.name
                    UNION
                    SELECT p.name, pp.pairid, p.tag, p.seq, p.vessel, p.well
                    FROM pairs AS pp, primer AS p
                    WHERE pp.right = p.name)
                    ORDER BY pairid;''')
                rows = cursor.fetchall()
            finally:
                self.db.close()
            return rows, ['primername', 'primerset', 'tag', 'sequence', 'vessel', 'well']
