#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# scan_region.py
#  
#  

import commands
import optparse
import os
import shutil
import subprocess
import sys
import re

BASE_DIR = ''

# Dictionary for databases with respective file path and output headers

def get_params(key):
    return {
        'phastconsele_100': (
            '--mce', BASE_DIR + 'phastConsElements100way.txt', 'pce100name',
            'pce100score'),
        'phastconsele_46': (
            '--mce', BASE_DIR + 'phastConsElements46way.txt', 'pce46name',
            'pce46score'),
        'phastconsele_46_pla': ('--mce', BASE_DIR + 'phastConsElements46wayPlacental.txt',
                                'pce46plname', 'pce46plscore'),
        'phastconsele_46_pri': ('--mce', BASE_DIR + 'phastConsElements46wayPrimates.txt', 'pce46prname',
                                'pce46prscore'),
        #'knowngene': ('--knowngene', BASE_DIR + 'knownGene.txt', 'knownGene gene_name', 'knownGene distance'),
        'kgxref': ('--knowngene --kgxref ' + BASE_DIR + 'kgXref.txt', BASE_DIR + 'knownGene.txt',
                   'kgsym', 'kgdist'),
        #'refgene': ('--refgene', BASE_DIR + 'refGene.txt', 'refGene gene_name', 'refGene distance'),
        'refcds': ('--refgene --refcds --reflink ' + BASE_DIR + 'refLink.txt', BASE_DIR + 'refGene.txt', 'refgsym',
                   'cdsdist'),
        'refexon': ('--refgene --refexon --reflink ' + BASE_DIR + 'refLink.txt', BASE_DIR + 'refGene.txt', 'refgsym',
                    'exondist'),
        'reflink': ('--refgene --reflink ' + BASE_DIR + 'refLink.txt', BASE_DIR + 'refGene.txt', 'refgsym',
                    'refgdist'),
        'evofold': ('--evofold', BASE_DIR + 'evofold.txt', 'evoscore'), #, 'evofold.norm_score'),
        'tfbs': ('--tfbs', BASE_DIR + 'tfbsConsSites.txt', 'tfbsscorez', 'tfbsscore'),
        'wgrna': ('--wgrna', BASE_DIR + 'wgRna.txt', 'miRNA'), #, 'wgRna.norm_score'),
        'segdup': ('--segdup', BASE_DIR + 'genomicSuperDups.txt', 'dupreg'), #, 'segdup.norm_score'),
        'dgv': ('--dgv', BASE_DIR + 'dgv.txt', 'dgvgain', 'dgvloss' ),
        'isca': ('--isca', BASE_DIR + 'isca.txt', 'annotation', 'pathscore')#Placental.txt'),
        #'phastcons_46_pri': ('--phastcons', BASE_DIR+'phastCons46wayPrimates.txt')
        #'phastcons_100': ('--phastcons', BASE_DIR+'phastCons100way.txt'),
        #'phastcons_46': ('--phastcons', BASE_DIR+'phastCons46way.txt'),
        #'phastcons_46_pla': ('--phastcons', BASE_DIR+'phastCons46way
    }.get(key)


def __main__():
    #print 'Parsing scan input options...'
    parser = optparse.OptionParser()                                # options parsing
    parser.add_option('--query', dest='queryfile')
    parser.add_option('--db', dest='db')
    parser.add_option('--resources', dest='resources')
    parser.add_option('--ex_path', dest='ex_path')
    #parser.add_option('--condense', action="store_true", dest='condense_query')
    parser.add_option('-p', '--pass_through', dest='pass_through_options', action='append', type="string")
    parser.add_option('--out', dest='outputFile')
    (options, args) = parser.parse_args()

    ex_path = options.ex_path + '/scan_region.pl'
    global BASE_DIR
    BASE_DIR = options.resources + '/'

    if len(args) > 0:
        parser.error('Wrong number of arguments')

    db = options.db.split(',')                                      # creating db list
    header = []
    #condense_flag = '--condense_query' if options.condense_query else ''
    query = open(options.queryfile, 'r')                            # open queryfile
    line = query.readline()                                         # read the queryfile's first line
    line = line.replace('\n', ' ')                                  # replace th \n
    obj = re.match('^(chr)?(\w+)\t(\d+)(\t(\d+))?(.*)', line, flags=0)  # cerco un match tra la regexp e la prima riga del file
    #if options.condense_query:                                          # se viene selezionata condense query elimino le colonne dopo l'end dall'header
    #   header.append('chromosome\tstart\tend')
    #else:
    if obj:                                                         # if the queryfile doesn't have a header create
        header.append('chromosome')
        a = line.count('\t')
        for i in range(a):
            header.append(' ')
    else:  #altrimenti inserisco la riga di intestazione della query
        header.append(line)
    query.close()
    for i, j in enumerate(db):
        params = get_params(j)
        if options.pass_through_options:
            common_options = ' '.join(options.pass_through_options)
            common_options = common_options.replace(",", " ")
        else:
            common_options = ''
        if len(params)== 4:
            line = params[2]
            header.append(line)
            line2 = params[3]
            header.append(line2)
        else:
            line = params[2]
            header.append(line)
        out = open('penn.txt', 'w')
        #if i > 0:
        #condense_flag = ''
        cmd = '%s --append --quiet %s %s %s %s' % (
            ex_path, params[0], options.queryfile, params[1], common_options)
        #print '%s' % (cmd)
        subprocess.check_call(cmd, stdout=out, shell=True)
        out.close()
        shutil.copy('penn.txt', 'temp')
        options.queryfile = 'temp'
    with open(options.outputFile, 'w')as out:
        with open('penn.txt', 'r') as dati:
            for flag in header:
                out.write(flag + '\t')
            out.write('\n')
            shutil.copyfileobj(dati, out)


if __name__ == '__main__':
    __main__()
