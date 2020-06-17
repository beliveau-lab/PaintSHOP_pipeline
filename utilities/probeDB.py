#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:11:26 2019

@author: hershe

A script to generate a database of refseq annotations with a row for each
probe that overlaps the region.

"""

import argparse
import glob
import timeit
import pybedtools
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd

def check_polarity(row):
    """"Checks polarity and flips probe sequence if reference is +."""
    if row[14] == '+':
        return str(Seq(row[3], IUPAC.unambiguous_dna).reverse_complement())
    else:
        return str(row[3])
    
def truncate_refseq(row):
    """"Truncates the Refseq column down to just accession."""
    version = row[12].split('_')[0] + '_' + row[12].split('_')[1]
    accession = version.split('.')[0]
    return accession

# command line interface
def getArgs(strInput=None):
    parser = argparse.ArgumentParser(
        'Concatenate the separate chromosome files of a probe set'
        'for an assembly and perform an intersection against '
        'the annotation file from Refseq. '
        'the output of the intersection is stored in a data frame.')
    parser.add_argument('-p', '--probeFolder', action='store', type=str, 
                        required=True,  help='A folder with probes')
    parser.add_argument('-a', '--annotationFile', action='store', type=str,
                        required=True, help="The name of the annotation file")
    parser.add_argument('-o', '--outputFile', action='store', type=str,
                        required=False, help="The name for output file")
    parser.add_argument('-f', '--file', action='store_true', default=False,
                        required=False, help="Run in file mode")
    
    return parser.parse_args()

def main():

    start_time = timeit.default_timer()
    
    args = getArgs()
    
    # store command line arguments
    folder = args.probeFolder
    annotation = args.annotationFile
    
    # determine output name
    if args.outputFile is None:
        out_name = 'probe_DB_out.csv'
    else:
        out_name = args.outputFile
    
    if args.file:
        probes = []
        with open(args.folder) as probe_file:
            for line in probe_file:
                probes.append(line.strip())
    else:
        # create a list of all files in directory
        files = glob.glob(folder + "/*")
        
        # creates a list of the all of the probes from each
        # probe file, each represented as a single string
        probes = []
        for f in files:
            with open(f) as probe_file:
                for line in probe_file:
                    probes.append(line.strip())
    
    # creates a bedtool object with the entire probe set for the assembly
    probe_bedtool = pybedtools.BedTool(probes)
    
    # creates a bedtool object for the refseq annotations
    # (snakemake input needs to be refseq BED file)
    refseq_bedtool = pybedtools.BedTool(annotation)
    
    # create output name for raw intersect result
    intersect_result = out_name.split('.')[0] + '_raw_intersect.bed'
    
    
    # wo is the same as -wo from BedTools, and writes the original
    # A and B entries plus the number of base pairs overlap of the
    # two features. f is the same as -f in BedTools, which is the
    # minimum overlap percentage accepted. 1 means 100% overlap
    probe_result = probe_bedtool.intersect(refseq_bedtool, wo=True, f=1,
                                           output = intersect_result)
    
    # remove since I parse output file
    del probe_result
    
    # read in the bed file of the results of intersection
    probes = pd.read_csv(intersect_result,
                         sep = '\t',
                         header = None)
    
    # check polarity of annotations, function flips probe sequence if necessary
    probes[3] = probes.apply(check_polarity, axis = 1)
    
    # convert refseq column to just accession
    probes[12] = probes.apply(truncate_refseq, axis = 1)
    
    # drop unnecessary columns
    probes.drop([9, 10, 11, 13, 14, 15], 
                axis = 1,
                inplace = True)
    
    # save data frame to disk
    probes.to_csv(out_name,
                  sep = '\t',
                  index = False,
                  index_label = False,
                  header = False)
    
    print('Program took %f seconds' % (timeit.default_timer() - start_time))
    
    
if __name__ == "__main__":
    main()
