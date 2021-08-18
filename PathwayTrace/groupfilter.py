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


import argparse


def get_options():
    parser = argparse.ArgumentParser(epilog="This script filters out identical sequences from a multifasta file.")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-i", "--input", default=None, type=str, required=True,
                          help="path to input fasta")
    required.add_argument("-o", "--output", default=None, type=str, required=True,
                          help="path to output fasta")
    optional.add_argument("-t", "--info", default=None, type=str, required=False,
                          help="file that contains all identical headers")
    args = parser.parse_args()
    main(args.input, args.output, args.info)


def read_fasta(path):
    fasta_dict = {}
    identical = {}
    count = 0
    with open(path, 'r') as infile:
        header = None
        seq = ''
        line = infile.readline()
        while line:
            if line[0] == '>':
                if header:
                    tmp = True
                    for entry in fasta_dict:
                        if fasta_dict[entry] == seq:
                            if entry in identical:
                                identical[entry].append(header)
                            else:
                                identical[entry] = [header]
                            tmp = False
                    if tmp:
                        fasta_dict[header] = seq
                header = line.lstrip('>').rstrip('\n')
                count = count + 1
                seq = ''
            else:
                seq = seq + line.rstrip('\n')
            line = infile.readline()
        fasta_dict[header] = seq
    return fasta_dict, identical, count


def write_fasta(path, fasta_dict):
    with open(path, 'w') as out:
        for entry in fasta_dict:
            out.write('>' + entry + '\n' + fasta_dict[entry] + '\n')


def write_identical(path, identical):
    with open(path, 'w') as out:
        for entry in identical:
            out.write(entry)
            for i in identical[entry]:
                out.write('\t' + i)
            out.write('\n')


def main(input, output, info):
    print('Reading Fasta...')
    fasta_dict, identical, count = read_fasta(input)
    write_fasta(output, fasta_dict)
    if info:
        write_identical(info, identical)
    print('#Original: ' + str(count) + '\n#Filtered: ' + str(len(fasta_dict)))


if __name__ == '__main__':
    get_options()
