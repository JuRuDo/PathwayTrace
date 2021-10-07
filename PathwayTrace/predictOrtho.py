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
#  along with PathwayTrace.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
import datetime
import os
import shutil
import sys
import subprocess
import argparse
import time

from Bio import SeqIO
import multiprocessing as mp
from pathlib import Path
from fcat.searchOrtho import isInt


#def checkRefspec(coreDir, coreSet, refspecList):
#    taxaFile = '%s/core_orthologs/%s/done.txt' % (coreDir, coreSet)
#    fcatfn.checkFileExist(taxaFile, '')
#    allTaxa = fcatfn.readFile(taxaFile)
#    invalid = []
#    for refspec in refspecList:
#        if not refspec in allTaxa:
#            invalid.append(refspec)
#    if len(invalid) > 0:
#        print('ERROR: Invalid refspec found: %s' % '; '.join(invalid))
#        print('Please use fcat.showTaxa for a list of valid reference taxa!')
#        sys.exit()


def get_refspec(refspecList, groupFa):
    coreSpec = []
    for s in SeqIO.parse(groupFa, 'fasta'):
        ref = s.id.split('|')[1]
        coreSpec.append(ref)
    for r in refspecList:
        if r in coreSpec:
            return r
    return None


def parseQueryFa(query, taxid, outDir):
    queryID = query.split('/')[-1].split('.')[0]
    queryIDtmp = queryID.split('@')
    if not (len(queryIDtmp) == 3 and isInt(queryIDtmp[1])):
        if taxid == '0':
            sys.exit('Query taxon does not have suitable ID format (e.g. HUMAN@9606@3). '
                     'Please provide its taxonomy ID additionaly using --taxid option.')
        else:
            addTaxon = 'fdog.addTaxon -f %s -i %s -o %s --replace' % (query, taxid, outDir)
            addTaxon = addTaxon + ' --noAnno'
            try:
                addTaxonOut = subprocess.run([addTaxon], shell=True, capture_output=True, check=True)
            except:
                sys.exit('Problem occurred while parsing query fasta file\n%s' % addTaxon)
            queryID = ''
            for line in addTaxonOut.stdout.decode().split('\n'):
                if "Species name" in line:
                    queryID = line.split('\t')[1]
                    # print('Query ID used by fCAT and fDOG: %s' % queryID)
            if len(queryID) == 0:
                sys.exit('Cannot identidy queryID!')
    else:
        Path('%s/%s/genome_dir/%s' % (outDir, queryID, queryID)).mkdir(parents=True, exist_ok=True)
        shutil.copy(query, '%s/%s/genome_dir/%s/%s.fa' % (outDir, queryID, queryID, queryID))
        os.chmod('%s/%s/genome_dir/%s/%s.fa' % (outDir, queryID, queryID, queryID), 0o755)
        checkedFile = open('%s/%s/genome_dir/%s/%s.fa.checked' % (outDir, queryID, queryID, queryID), 'w')
        now = datetime.datetime.now()
        checkedFile.write(now.strftime("%Y-%m-%d %H:%M:%S"))
        checkedFile.close()
    return queryID


def prepareJob(coreDir, coreSet, queryID, refspecList, outDir, blastDir, force):
    fdogJobs = []
    ignored = []
    groupRefspec = {}
    hmmPath = coreDir + '/' + coreSet
    groups = os.listdir(hmmPath)
    if len(groups) > 0:
        searchPath = '%s/%s/%s/genome_dir' % (outDir, coreSet, queryID)
        # create single fdog job for each core group
        for groupID in groups:
            if os.path.isdir(hmmPath + '/' + groupID):
                groupFa = '%s/%s/%s/%s.fa' % (coreDir, coreSet, groupID, groupID)
                # check refspec
                refspec = checkRefspec(refspecList, groupFa)
                if refspec == '':
                    ignored.append(groupID)
                else:
                    outPath = '%s/%s/%s/fdogOutput/%s' % (outDir, coreSet, queryID, refspec)
                    if not os.path.exists('%s/%s/%s.phyloprofile' % (outPath, groupID, groupID)) or force:
                        fdogJobs.append([groupFa, groupID, refspec, outPath, blastDir, hmmPath, searchPath, force])
                    groupRefspec[groupID] = refspec
    else:
        sys.exit('No core group found at %s' % (coreDir + '/' + coreSet))
    return fdogJobs, ignored, groupRefspec


def prepare_data(args):
    if args.outDir:
        Path(args.outDir).mkdir(parents=True, exist_ok=True)
    else:
        args.outDir = os.getcwd()
    if args.cpus >= mp.cpu_count():
        args.cpus = mp.cpu_count()-1

    args.queryID = parseQueryFa(os.path.abspath(args.querySpecies), str(args.taxid), args.outDir)
    print('Query ID used by fDOG: %s' % args.queryID)
    return args


def runFdog(args):
    (seqFile, seqName, refSpec, outPath, blastPath, hmmPath, searchPath, force) = args
    fdog = 'fdog.run --seqFile %s --seqName %s --refspec %s --outpath %s --blastpath %s --hmmpath %s --searchpath %s ' \
           '--fasoff --reuseCore --checkCoorthologsRef --cpu 1' \
           % (seqFile, seqName, refSpec, outPath, blastPath, hmmPath, searchPath)
    if force:
        fdog = fdog + ' --force'
    fdog = fdog + ' > /dev/null 2>&1'
    try:
        subprocess.run([fdog], shell=True, check=True)
    except:
        print('\033[91mProblem occurred while running fDOG for \'%s\' core group\033[0m\n%s' % (seqName, fdog))


def main():
    version = '0.1'
    parser = argparse.ArgumentParser(description='You are running PathwayTrace version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-c', '--hmmdir', action='store', type=str, required=True,
                          help='Name of the group')
    required.add_argument('-f', '--groupFa', action='store', type=str, required=True,
                          help='Name of the group')
    required.add_argument('-r', '--refspecList', required=True, nargs='+',
                          help='List of reference species')
    required.add_argument('-q', '--querySpecies', action='store', type=str, required=True,
                          help='Path to gene set for species of interest')
    optional.add_argument('-o', '--outDir', action='store', default=None, help='Path to output directory')
    optional.add_argument('-b', '--blastDir', action='store', default=None,
                          help='Path to BLAST directory of all core species')
    optional.add_argument('-i', '--taxid', action='store', default=0, type=int,
                          help='Taxonomy ID of gene set for species of interest')
    optional.add_argument('--cpus', action='store', default=1, type=int,
                          help='Number of CPUs used for annotation. Default = 1')
    optional.add_argument('--force', action='store_true', default=False, help='Force overwrite existing data')
    args = parser.parse_args()

    start = time.time()
    args = prepare_data(args)
    args.refspec = get_refspec(args.refspecList, args.groupFa)
    args.group_id = args.groupFa.split('/')[-1].split('.')[0]
    args.searchPath = '%s/%s/genome_dir' % (args.outDir, args.queryID)
    fdog_args = (args.groupFa, args.group_id, args.refspec, args.outDir, args.blastDir, args.hmmdir,
                 args.searchPath, args.force)
    print('##### Searching for orthologs...')
#    runFdog(fdog_args)
    end = time.time()
    print('Finished in ' + '{:5.3f}s'.format(end-start))


if __name__ == '__main__':
    main()
