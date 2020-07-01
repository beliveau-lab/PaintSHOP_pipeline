import argparse
import pybedtools
import timeit
import random
import os


# Create a function to caculated shared exons.
def collapseIsos(a, b):
    '''Collapses isoforms to only shared sequence.'''

    # Make random ints to use for file names.
    r1 = random.randint(0, 10000000)
    r2 = random.randint(0, 10000000)

    # Create input .bed files.
    with open('/tmp/%s.bed' % r1, 'w') as f:
        f.write(a)
    with open('/tmp/%s.bed' % r2, 'w') as f:
        f.write(b)

    # Use pybedtools to call coverageBed
    a_bed = pybedtools.BedTool('/tmp/%s.bed' % r1)
    b_bed = pybedtools.BedTool('/tmp/%s.bed' % r2)
    cov = a_bed.coverage(b_bed, d=True)

    # Extract covergeBed ouput to list.
    with open(cov.fn, 'r') as f:
        cov_list = [x.strip() for x in f]

    # Create list to hold output.
    return_list = []

    # Determine is there is an interval shared by >1 isform.
    max_overlaps = max([int(x.split('\t')[7]) for x in cov_list if x != ''])

    # If so, take intervals shared by the greatest number of isoforms.
    if max_overlaps > 1:

        # Process coverageBed output and collapse isoforms to
        # only shared sequence.
        # Initialize position tracking variables.
        start = 0
        counter = 0

        # Loop through coverageBed ouput.
        for i in range(0, len(cov_list)-1, 1):

            # If position covers all isoforms and no interval is being tracked,
            # initialize tracking.
            if int(cov_list[i].split('\t')[7]) == max_overlaps \
            and start == 0:
                start = int(cov_list[i].split('\t')[1]) - 1 \
                + int(cov_list[i].split('\t')[6])

            # If position covers all isoforms and an interval is being tracked,
            # continue counting interval coords.
            elif int(cov_list[i].split('\t')[7]) == max_overlaps and start != 0:

                # Append interval to output if at end of coverageBed output.
                if i == len(cov_list) - 2:
                    return_list.append('%d\t%d\t%d' % \
                    (start, start + counter + 2, max_overlaps))
                else:
                    counter +=1

            # If position no longer covering all isoforms, append previously
            # interval to output and reset position tracking variables.
            else:
                if start != 0:
                    return_list.append('%d\t%d\t%d' % \
                    (start, start + counter + 1, max_overlaps))
                    start = 0
                    counter = 0

    # Else, return the original annotation as is.
    else:
        b_exons = b.split('\n')
        for i in range(0, len(b_exons), 1):
            return_list.append('%s\t%s\t1' % (b_exons[i].split('\t')[1],
            b_exons[i].split('\t')[2]))

    # Clean up temp files.
    pybedtools.helpers.cleanup()
    os.remove('/tmp/%s.bed' % r1)
    os.remove('/tmp/%s.bed' % r2)

    # Return list of flattened isoforms.
    return return_list


def main():
    '''Use pybedtools wrapper to collapse transcript isoforms into only shared
    sequence. Transcripts with 1 isoform written as is. Transcripts without
    shared sequenced among all isoforms have most parsimonious set of intervals
    returned.'''

    # Start timer.
    startTime = timeit.default_timer()

    # Allow user to input parameters on command line.
    userInput = argparse.ArgumentParser()
    userInput.add_argument('-f', '--file', action='store',
                           help='The GTF file contaning transcript annotations')
    userInput.add_argument('-o', '--output', action='store', default=None,
                           type=str,
                           help='Specify the name prefix of the output file')

    # Process argparse input.
    args = userInput.parse_args()

    # Determine the stem of the input filename.
    fileName = str(args.file).split('.')[0]

    # Open and extract contents of annotation file
    with open(args.file, 'r') as f:
        file_read = [line.strip() for line in f \
        if '_' not in line.strip().split('\t')[0]]

    # Make list to hold output.
    outList = []

    # Extract unique gene IDs to a dictionary as keys. Initialize with empty lists
    # for values.
    geneIDs = {x.split('\t')[0] + "|" + \
    x.split('\t')[8].split(' ')[1].strip('"').strip('";'): [] \
    for x in file_read}

    # Iterate through input file mapping isoform information to dictionary.
    for i in range(0, len(file_read), 1):
        ID = file_read[i].split('\t')[0] + "|" + \
        file_read[i].split('\t')[8].split(' ')[1].strip('"').strip('";')
        entry = file_read[i].split('\t')[2]
        if entry == 'exon':
            geneIDs[ID].append('%s\t%s\t%s\t%s\t0\t%s'
            % (file_read[i].split('\t')[0], \
            file_read[i].split('\t')[3], \
            file_read[i].split('\t')[4], \
            file_read[i].split('\t')[8].split(' ')[3].strip('"').strip('";'), \
            file_read[i].split('\t')[6]))

    # Iterate through dictionary and perform isoform flattening.
    for key in geneIDs.keys():

        print (key)

        # Extract info about chromomosome and strand.
        chrom = geneIDs[key][0].split('\t')[0]
        strand = geneIDs[key][0].split('\t')[5]

        # Extract start and stop positions of exons.
        starts = [int(x.split('\t')[1]) for x in geneIDs[key]]
        stops = [int(x.split('\t')[2]) for x in geneIDs[key]]

        # Compute maximum span covered by all isoforms.
        _a = '%s\t%d\t%d\t%s' % (chrom, min(starts), max(stops), \
        '\t'.join(geneIDs[key][0].split('\t')[3:]))

        # Create formatted string of all exon positions.
        _b = '\n'.join(geneIDs[key])

        # Create a list of unique isoforms associated with the gene.
        isos = list({x.split('\t')[3] for x in geneIDs[key]})

        # If there is >1 isoform, call flattening function.
        if len(isos) > 1:
            shared_intervals = collapseIsos(_a, _b)

            # Check to see if there was >1 shared interval
            if int(shared_intervals[0].split('\t')[2]) > 1:

                # Save same shared exon positions for each isoform to output list.
                for i in range(0, len(isos), 1):
                    for j in range(0, len(shared_intervals), 1):
                        outList.append('%s\t%s\t%s\t%s\t0\t%s\t%s\t%s' % \
                        (chrom, shared_intervals[j].split('\t')[0], \
                        shared_intervals[j].split('\t')[1], isos[i], strand, \
                        shared_intervals[j].split('\t')[2], key))

            # Else just re-write original annotation info.
            else:
                for i in range(0, len(geneIDs[key]), 1):
                    outList.append(geneIDs[key][i] + '\t' + '0' + '\t' + key)

        # If there is only 1 isoform, save the info about exons to output list.
        else:
            for i in range(0, len(geneIDs[key]), 1):
                outList.append(geneIDs[key][i] + '\t' + '1' + '\t' + key)

    # Determine the name of the output file.
    if args.output is None:
        outName = fileName + '_iso_flatten'
    else:
        outName = args.output

    # Create the output file.
    with open('%s.bed' % outName, 'w') as f:
        f.write('\n'.join(outList))

    # Print wall-clock runtime to terminal.
    print ('Program took %f seconds' % (timeit.default_timer() - startTime))


if __name__ == '__main__':
    main()
