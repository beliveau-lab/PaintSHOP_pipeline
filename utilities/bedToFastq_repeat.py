#!/usr/bin/env python
# --------------------------------------------------------------------------
# OligoMiner
# bedToFastq.py
# 
# (c) 2016 Molecular Systems Lab
# Wyss Institute for Biologically-Inspired Engineering
# Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# --------------------------------------------------------------------------

# Specific script name.
scriptName = 'bedToFastq'

# Specify script version.
Version = '1.7'

def convertBedToFastq(inputFile, outNameVal):
    """Converts a .bed file to a .fastq file."""

    # Determine the stem of the input filename.
    fileName = str(inputFile).split('.')[0]

    # Open input file for reading.
    with open(inputFile, 'r') as f:
        file_read = [line.strip() for line in f]

    # Create list to hold output.
    outList = []

    # A list to hold arbitrary quality scores for each base in the candidate
    # probe.
    quals = ['~' * len(file_read[i].split('\t')[3]) \
             for i in range(len(file_read))]

    # Parse out probe information and write into .fastq format.
    for i in range(0, len(file_read), 1):
        chrom = file_read[i].split('\t')[0]
        start = file_read[i].split('\t')[1]
        stop = file_read[i].split('\t')[2]
        probeSeq = file_read[i].split('\t')[3]
        # include repeat flag (set to 0 for "no repeat")
        outList.append('@%s:%s-%s|0\n%s\n+\n%s' \
                       % (chrom, start, stop, probeSeq, quals[i]))

    # Determine the name of the output file.
    if outNameVal is None:
        outName = fileName
    else:
        outName = outNameVal

    # Create the output file.
    output = open(outName, 'w')

    # Write the output file
    output.write('\n'.join(outList))
    output.close()

# call the function with the snakemake arguments
convertBedToFastq(snakemake.input[0], snakemake.output[0])


