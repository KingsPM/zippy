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
from zippylib.primer import Primer

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
            cursor.execute('CREATE TABLE IF NOT EXISTS primer(seq TEXT PRIMARY KEY, tm REAL, gc REAL)')
            cursor.execute('CREATE TABLE IF NOT EXISTS design(seq TEXT, user TEXT, dateadded DATE, FOREIGN KEY(seq) REFERENCES primer(seq) ON UPDATE CASCADE)')
            cursor.execute('CREATE TABLE IF NOT EXISTS target(seq TEXT, chrom TEXT, pos TEXT, reverse BOOLEAN, FOREIGN KEY(seq) REFERENCES primer(seq) ON UPDATE CASCADE)')
            cursor.execute('CREATE TABLE IF NOT EXISTS pairs(pairid TEXT PRIMARY KEY, left TEXT, right TEXT, FOREIGN KEY(left) REFERENCES primer(seq) ON UPDATE CASCADE, FOREIGN KEY(right) REFERENCES primer(seq) ON UPDATE CASCADE)')
            cursor.execute('CREATE TABLE IF NOT EXISTS status(pairid TEXT, status INT, stock TEXT, FOREIGN KEY(pairid) REFERENCES pairs(pairid) ON UPDATE CASCADE)')
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
            cursor.execute('SELECT pairs.*, p1.tm, p2.tm FROM pairs LEFT JOIN primer as p1 ON pairs.left = p1.seq LEFT JOIN primer as p2 ON pairs.right = p2.seq')
            rows = cursor.fetchall()
        finally:
            self.db.close()
        return "\n".join([ '{:<10} {:<30} {:<30} {:.1f} {:.2f}'.format(*row) for row in rows ])

    def getPrimers(self, paired=True):
        raise NotImplementedError

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
                    cursor.execute('''INSERT OR REPLACE INTO target(seq,chrom,pos,reverse) VALUES(?,?,?,?)''', \
                        (p.seq, l[0], l[1], l[2]))
                cursor.execute('''INSERT OR REPLACE INTO design(seq,user,dateadded) VALUES(?,?,?)''', \
                    (p.seq, self.user, datetime.datetime.now()))
            self.db.commit()
        finally:
            self.db.close()
        return

    def addPair(self, *pairs):
        '''adds primer pairs (and individual primers)'''
        # add primers
        flat = list(set(list(flatten(pairs))))
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
                if not pairid:  # fallback
                    pairid = "X"+md5(left.seq+'_'+right.seq).hexdigest()
                cursor.execute('''INSERT OR REPLACE INTO pairs(pairid,left,right) VALUES(?,?,?)''', \
                    (pairid, left.seq, right.seq))
            self.db.commit()
        finally:
            self.db.close()
        return

    def query(self, variant, cfg):
        raise NotImplementedError
        # return primer pairs that would match

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
