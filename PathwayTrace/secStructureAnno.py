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


from Bio import SeqIO
import subprocess
import os
from pathlib import Path
import multiprocessing as mp
from tqdm import tqdm
import argparse


def parse_porter(infile):
    struc = ''
    with open(infile, 'r') as toparse:
        for line in toparse.readlines():
            if not line[0] == '#':
                struc = struc + line.split()[2]
    return struc


def prepare_annojobs(infile, tmppath, porterpath, cpus):
    annojobs = []
    for seq in SeqIO.parse(infile, 'fasta'):
        annojobs.append([seq.id, str(seq.seq), tmppath, porterpath, str(cpus)])
    return annojobs


def run_porter(infile, porter, cpus):
    cmd = 'python ' + porter + ' -i ' + infile + '.fasta --cpu ' + cpus + ' --fast'
    subprocess.run([cmd], shell=True, capture_output=True, check=True)
    ss3 = parse_porter(infile + '.fasta.ss3')
    ss8 = parse_porter(infile + '.fasta.ss8')
    for i in ['.fasta.ss3', '.fasta.ss8', '.fasta.psi', '.fasta.blastpgp', '.hhr']:
        remove_tmp(infile + i)
    return ss3, ss8


def remove_tmp(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)


def make_tmp_fasta(header, seq, path):
    with open(path + '/' + header + '.fasta', 'w') as out:
        out.write('>' + header + '\n' + seq)


def run_anno(inpath, outpath, tmppath, cpus, parallel, porterpath):
    name = ''.join(inpath.split('/')[-1].split('.')[0:-1])
    Path(outpath).mkdir(parents=True, exist_ok=True)
    Path(tmppath + name + '/').mkdir(parents=True, exist_ok=True)
    joblist = prepare_annojobs(inpath, tmppath + name + '/', porterpath, cpus)
    out = []
    pool = mp.Pool(parallel)
    try:
        for _ in tqdm(pool.imap_unordered(run_anno_single, joblist), total=len(joblist), mininterval=1.0):
            out.append(_)
    except subprocess.CalledProcessError:
        raise 'Something went wrong while running Porter.'
    pool.close()
    write_output(out, outpath + '/' + name)
    Path(tmppath + name + '/').rmdir()


def write_output(out, outpath):
    ss3 = open(outpath + '.ss3', 'w')
    ss8 = open(outpath + '.ss8', 'w')
    for seq in out:
        ss3.write('>' + seq[0] + '\n' + seq[1] + '\n')
        ss8.write('>' + seq[0] + '\n' + seq[2] + '\n')
    ss3.close()
    ss8.close()


def run_anno_single(args):
    header, seq, tmppath, porter, cpus = args
    make_tmp_fasta(header, seq, tmppath)
    ss3, ss8 = run_porter(tmppath + header, porter, cpus)
    if os.path.exists(tmppath + header + '.fasta'):
        os.remove(tmppath + header + '.fasta')
    return header, ss3, ss8


def main():
    parser = argparse.ArgumentParser(epilog="This script filters out identical sequences from a multifasta file.")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", default='.', type=str, required=True,
                          help="path to input fasta")
    required.add_argument("-o", "--outPath", default='.', type=str, required=True,
                          help="path to output directory")
    optional.add_argument("-t", "--tmp", default='.', type=str, required=False,
                          help="")
    optional.add_argument("-c", "--cpus", default=mp.cpu_count()-1, type=int, required=False,
                          help="number of cpus for each instance of Porter to use")
    optional.add_argument("--parallel", default=1, type=int, required=False,
                          help="number of parallel runs of porter")
    optional.add_argument("-p", "--porter", default=None, type=str, required=False,
                          help="Path to porter")
    args = parser.parse_args()
    run_anno(args.input, args.outPath, args.tmp, args.cpus, args.parallel, args.porter)


if __name__ == '__main__':
    main()
