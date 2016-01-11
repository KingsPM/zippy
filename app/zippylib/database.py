#!/usr/bin/env python

import sys, os, re, ast
import datetime
import json
import hashlib
import sqlite3
import fnmatch
import copy
from collections import defaultdict
from zippylib import flatten, commonPrefix
from zippylib.primer import Primer, Locus

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
            cursor.execute('CREATE TABLE IF NOT EXISTS primer(name TEXT PRIMARY KEY, seq TEXT, tm REAL, gc REAL, FOREIGN KEY(seq) REFERENCES target(seq));')
            cursor.execute('CREATE TABLE IF NOT EXISTS target(seq TEXT, chrom TEXT, position INT, reverse BOOLEAN);')
            cursor.execute('CREATE TABLE IF NOT EXISTS pairs(pairid TEXT PRIMARY KEY, left TEXT, right TEXT, chrom TEXT, start INT, end INT, FOREIGN KEY(left) REFERENCES primer(seq), FOREIGN KEY(right) REFERENCES primer(seq));')
            cursor.execute('CREATE TABLE IF NOT EXISTS status(pairid TEXT PRIMARY KEY, status INT, dateadded TEXT, FOREIGN KEY(pairid) REFERENCES pairs(pairid));')
            cursor.execute('CREATE INDEX IF NOT EXISTS seq_index_in_target ON target(seq);')
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
            cursor.execute('SELECT pairs.*, p1.tm, p2.tm FROM pairs LEFT JOIN primer as p1 ON pairs.left = p1.seq LEFT JOIN primer as p2 ON pairs.right = p2.seq;')
            rows = cursor.fetchall()
        finally:
            self.db.close()
        return "\n".join([ '{:<10} {:<30} {:<30} {:.1f} {:.2f}'.format(*row) for row in rows ])

    def getPrimers(self, paired=True):
        raise NotImplementedError
        # try:
        #     self.db = sqlite3.connect(self.sqlite)
        # except:
        #     raise
        # else:
        #     cursor = self.db.cursor()
        #     cursor.execute('''SELCECT FROM primer''')



        #     self.db.commit()
        # finally:
        #     self.db.close()
        # return

    '''adds list of primers to database'''
    def addPrimer(self, *primers):
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            for p in primers:
                cursor = self.db.cursor()
                cursor.execute('''INSERT OR REPLACE INTO primer(name,seq,tm,gc) VALUES(?,?,?,?)''', \
                    (p.name, p.seq, p.tm, p.gc))
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
            cursor = self.db.cursor()
            # add pairs
            for p in pairs:
                left, right = p[0],p[1]
                # find common substring in name for automatic naming
                pairid = commonPrefix(left.name,right.name)
                chrom = left.targetposition.chrom
                start = left.targetposition.offset
                end = right.targetposition.offset+right.targetposition.length
                if not pairid:  # fallback
                    pairid = "X"+md5(left.seq+'_'+right.seq).hexdigest()
                cursor.execute('''INSERT OR REPLACE INTO pairs(pairid,left,right,chrom,start,end) VALUES(?,?,?,?,?,?)''', \
                    (pairid, left.seq, right.seq, chrom, start, end))
                cursor.execute('''INSERT OR REPLACE INTO status(pairID,dateadded) VALUES(?,?)''', \
                    (pairid, datetime.datetime.now()))
            self.db.commit()
        finally:
            self.db.close()
        return

    def query(self, variant, flank):
        '''returns suitable primer pairs for the specified loci'''
        try:
            self.db = sqlite3.connect(self.sqlite)
        except:
            raise
        else:
            chrom = variant.chrom
            chromStart = variant.chromStart
            chromEnd = variant.chromEnd
            # print chrom
            # print chromStart
            # print chromEnd
            # print flank
            cursor = self.db.cursor()
            cursor.execute('''SELECT *
                FROM pairs AS p
                LEFT JOIN target AS t1 ON p.left = t1.seq
                LEFT JOIN target AS t2 ON p.right = t2.seq
                LEFT JOIN status AS s ON p.pairID = s.pairID
                WHERE t1.chrom = t2.chrom
                AND t1.chrom = ?
                AND t1.position + length(t1.seq) + ? <= ?
                AND t2.position - ? >= ?;''', (chrom, flank, chromStart, flank, chromEnd))

            rows = cursor.fetchall()
            # cursor.execute('''SELECT COUNT *
                # FROM target
                # WHERE seq = seq;''')


                # AS t1, t2
                # LEFT JOIN pairs AS p ON p.left = t1.seq
                # LEFT JOIN pairs AS p ON p.right = t2.seq
                # WHERE t1.chrom = t2.chrom
                # AND t1.chrom = ?
                # AND t1.position + length(t1.seq) + ? <= ?
                # AND t2.position - ? >= ?                ;''', ())
            leftLoci = cursor.fetchall()
            # print rows

        finally:
            self.db.close()
        # return primer pairs that would match
        primerPairs = []
        for row in rows:
            # print row,'\n'
            name = row[0]
            leftSeq = row[6]
            rightSeq = row[10]
            leftTargetposition = Locus(row[7], row[8], len(row[6]), True if row[9] else False) #(chrom,offset,length,reverse)
            rightTargetposition = Locus(row[11], row[12], len(row[10]),True if row[13] else False)
            # print name
            # print leftSeq
            # print leftTargetposition
            # print rightSeq
            # print rightTargetposition
            leftPrimer = Primer(name+'_left', leftSeq, leftTargetposition) #(name,seq,targetposition(locus),tm) {, calcProperties(leftSeq)}
            rightPrimer = Primer(name+'_right', rightSeq, rightTargetposition)
            leftPrimer.calcProperties()
            rightPrimer.calcProperties()

            status = row[15]
            dateAdded = row[16]
            # print leftPrimer
            # print rightPrimer
            primerPairs.append([[leftPrimer, rightPrimer], status])
        return primerPairs


    def dump(self,what,**kwargs):
        raise NotImplementedError
        if what=='amplicons':
            # dump amplicons (all possible)
            try:
                self.db = sqlite3.connect(self.sqlite)
            except:
                raise
            else:
                cursor = self.db.cursor()
                # FWD strand
                cursor.execute('''SELECT DISTINCT t1.chrom, t1.pos, t2.pos+length(t2.seq), p.pairid
                    FROM pairs AS p
                    LEFT JOIN target AS t1 ON p.left = t1.seq
                    LEFT JOIN target AS t2 ON p.right = t2.seq
                    WHERE t1.chrom = t2.chrom
                    AND t2.pos+length(t2.seq)-t1.pos >= ?
                    AND t2.pos+length(t2.seq)-t1.pos <= ?
                    AND t1.reverse = 0
                    AND t2.reverse = 1
                    ORDER BY t1.chrom, t1.pos;''', (kwargs['size'][0],kwargs['size'][1]))
                rows = cursor.fetchall()
                # REV strand
                cursor.execute('''SELECT DISTINCT t1.chrom, t1.pos, t2.pos+length(t2.seq), p.pairid
                    FROM pairs AS p
                    LEFT JOIN target AS t1 ON p.right = t1.seq
                    LEFT JOIN target AS t2 ON p.left = t2.seq
                    WHERE t1.chrom = t2.chrom
                    AND t2.pos+length(t2.seq)-t1.pos >= ?
                    AND t2.pos+length(t2.seq)-t1.pos <= ?
                    AND t1.reverse = 0
                    AND t2.reverse = 1
                    ORDER BY t1.chrom, t1.pos;''', (kwargs['size'][0],kwargs['size'][1]))
                rows += cursor.fetchall()
            finally:
                self.db.close()
            return rows, ('chrom','chromStart','chromEnd','name')  # rows and colnames
            # return [ '{}\t{}\t{}\t{}'.format(*row) for row in rows ]

''' SQLite based file and checkpoint manager
    def getTasks(self,status):
        # returns a list of files in a bucket, ordered by timestamp
        self.db = sqlite3.connect(self.sqlite)
        cursor = self.db.cursor()
        cursor.execute(SELECT task FROM taskregistry WHERE status = ?', (status,))
        rows = cursor.fetchall()
        self.db.close()
        return [ x[0] for x in rows ]

    ####################
    ### FILE BUCKETS ###
    ####################
    def update(self,filename,bucket='default'):
        # parse filenames
        if type(filename) is not list and type(filename) is not tuple:
            inserts = [ (filename, bucket, datetime.datetime.now()) ]
        else:
            inserts = [ (fn, bucket, datetime.datetime.now()) for fn in filename ]
        # execute query
        self.db = sqlite3.connect(self.sqlite)
        cursor = self.db.cursor()
        try:
            cursor.executemany('INSERT OR REPLACE INTO files(filename,bucket,created) VALUES(?,?,?)', inserts)
            self.db.commit()
        except:
            print >> sys.stderr, filename
            print >> sys.stderr, bucket
            print >> sys.stderr, datetime.date.today()
            raise
        finally:
            self.db.close()
        return


    #####################
    ### CHECKPOINTING ###
    #####################
    def getMetricgroup(self, fun):
        # returns metric group from function or fun if not available
        self.db = sqlite3.connect(self.sqlite)
        cursor = self.db.cursor()
        cursor.execute('SELECT DISTINCT metricgroup FROM checkpoints where task = ?', (fun,))
        rows = cursor.fetchall()
        self.db.close()
        try:
            assert len(rows)<2
        except AssertionError:
            print rows
            raise Exception('MultipleMetricgroups')
        return rows[0][0]  # return first field of first row

    def failedCheckpoints(self, exclude=''):
        # get passed = 0
        self.db = sqlite3.connect(self.sqlite)
        cursor = self.db.cursor()
        cursor.execute('SELECT DISTINCT metric FROM checkpoints where passed = ? and metricgroup != ? and task != ?', (0,exclude,exclude,))
        rows = cursor.fetchall()
        self.db.close()
        return [ x for x in rows ]

    def checkpoints(self):
        self.db = sqlite3.connect(self.sqlite)
        cursor = self.db.cursor()
        cursor.execute('SELECT DISTINCT assessed, token, passed, task, metric, comment FROM checkpoints WHERE token != "NULL" ORDER BY assessed')
        rows = cursor.fetchall()
        self.db.close()
        return rows

    def setCheckpoint(self,fun,met,state,comment='NULL',metgrp=None,token='NULL'):
        if not metgrp:
            metgrp = fun
        self.db = sqlite3.connect(self.sqlite)
        state = 1 if state else 0
        try:
            cursor = self.db.cursor()
            cursor.execute('INSERT OR REPLACE INTO checkpoints(task,token,metricgroup,metric,passed,assessed,comment) VALUES(?,?,?,?,?,?,?)', \
                (fun, token, metgrp, met, state, datetime.datetime.now(), comment))
            self.db.commit()
        except:
            print >> sys.stderr, (fun, metgrp, state, datetime.datetime.now(), comment)
            raise
        finally:
            self.db.close()
        return

    def resetCheckpoints(self):
        self.db = sqlite3.connect(self.sqlite)
        try:
            cursor = self.db.cursor()
            cursor.execute('DELETE FROM checkpoints')
            self.db.commit()
        except:
            raise
        finally:
            self.db.close()
        return

    # get failed checkpoints
    # gets output files and deletes/unlinks them
    # resets checkpoints for task
    def resetQCfail(self):
        token = re.compile(r'\w{8}$')
        self.db = sqlite3.connect(self.sqlite)
        # get token of failed taks
        cursor = self.db.cursor()
        cursor.execute('SELECT DISTINCT token FROM checkpoints WHERE passed != 1')
        rows = cursor.fetchall()
        failedToken = [ x for x in snappylib.utils.flatten(rows) if token.match(x) ]
        for t in failedToken:
            # get args and delete
            cursor = self.db.cursor()
            cursor.execute('SELECT DISTINCT args FROM runtoken where token = ?', (t,))
            rows = cursor.fetchall()
            for argstring in snappylib.utils.flatten(rows):
                args = ast.literal_eval(argstring)
                for outfile in args[1]:
                    try:
                        os.unlink(outfile)
                    except OSError:
                        pass
                    except:
                        raise
            # remove from registry, runtoken and checkpoints
            cursor = self.db.cursor()
            cursor.execute('DELETE FROM taskregistry WHERE token = ?', (t,))
            cursor.execute('DELETE FROM runtoken WHERE token = ?', (t,))
            cursor.execute('DELETE FROM checkpoints WHERE token = ?', (t,))
            self.db.commit()
        self.db.close()
        return
'''
