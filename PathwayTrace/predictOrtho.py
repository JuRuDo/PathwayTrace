#!/bin/env python

#######################################################################
# Copyright (C) 2021 Julian Dosch
#
# This file is part of PathwayTrace.
#
#  PathwayTrace is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PathwayTrace is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ProofReader.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import os
import sys
import subprocess
from Bio import SeqIO
import multiprocessing as mp
from pathlib import Path
from fcat import showTaxa
from fcat import searchOrtho
import fcat.functions as fcatfn


def checkRefspec(coreDir, coreSet, refspecList):
    taxaFile = '%s/core_orthologs/%s/done.txt' % (coreDir, coreSet)
    fcatfn.checkFileExist(taxaFile, '')
    allTaxa = fcatfn.readFile(taxaFile)
    invalid = []
    for refspec in refspecList:
        if not refspec in allTaxa:
            invalid.append(refspec)
    if len(invalid) > 0:
        print('ERROR: Invalid refspec found: %s' % '; '.join(invalid))
        print('Please use fcat.showTaxa for a list of valid reference taxa!')
        sys.exit()


def fcat(args):
    outDir = args.outDir
    if outDir == '':
        outDir = os.getcwd()
    else:
        Path(outDir).mkdir(parents=True, exist_ok=True)
    annoDir = args.annoDir
    if (annoDir == ''
            or annoDir == '%s/weight_dir' % os.path.abspath(args.coreDir)
            or annoDir == '%s/weight_dir/' % os.path.abspath(args.coreDir)):
        annoDir = '%s/fcatOutput/%s/weight_dir' % (outDir, args.coreSet)
        Path(annoDir).mkdir(parents=True, exist_ok=True)
        # annoDir = '%s/weight_dir' % args.coreDir
    annoDir = os.path.abspath(annoDir)
#    for annoFile in glob.glob('%s/weight_dir/*.json' % args.coreDir):
#        annoFileName = annoFile.split('/')[-1]
#        if not os.path.exists('%s/%s' % (annoDir, annoFileName)):
#            try:
#                os.symlink('%s/weight_dir/%s' % (args.coreDir, annoFileName), '%s/%s' % (annoDir, annoFileName))
#            except FileExistsError:
#                os.remove('%s/%s' % (annoDir, annoFileName))
#                os.symlink('%s/weight_dir/%s' % (args.coreDir, annoFileName), '%s/%s' % (annoDir, annoFileName))
    cpus = args.cpus
    if cpus >= mp.cpu_count():
        cpus = mp.cpu_count()-1

    # get queryID
    if args.annoQuery == '':
        print('Annotation for %s not given! It might take a while for annotating...' % args.querySpecies)
    (doAnno, queryTaxId) = searchOrtho.checkQueryAnno(args.annoQuery, annoDir, args.taxid, args.querySpecies)
    args.queryID = searchOrtho.parseQueryFa(args.coreSet, os.path.abspath(args.querySpecies), args.annoQuery,
                                            str(args.taxid), outDir, doAnno, annoDir, cpus)
    print('Query ID used by fCAT and fDOG: %s' % args.queryID)
    if doAnno == False:
        if os.path.exists( '%s/query_%s.json' % (annoDir, queryTaxId)):
            os.remove('%s/query_%s.json' % (annoDir, queryTaxId))

    # calculate group specific cutoffs
    # print('##### Calculating group specific cutoffs...')
    # fcatC.calcGroupCutoff(args)

    # check for valid refspec
    checkRefspec(args.coreDir, args.coreSet, str(args.refspecList).split(","))

    # search for orthologs and create phylognetic profile files
    print('##### Searching for orthologs...')
    searchOrtho.searchOrtho(args)





def


def main():
    pass
