##### Import modules     
import pybedtools
import glob
from pybedtools import BedTool
import pandas as pd
import csv
import argparse
import logging
import timeit
import datetime

##### Define arguments locallly to test script (comment out when bash submitting)
#bedgraph="GFP__KLF4_sc.sorted.bg"
#threshold=2
#min_length=50
#inter_peak_distance=50
#merge_close_peaks=True
#max_length=10000
#generate_ID=True

##### Define parseArguments to use positional variables in bash submission
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('--bedgraph', help='bedgraph', type=str)

    # Optional arguments with defaults
    parser.add_argument('--threshold', help='threshold', type=int, default=2)
    parser.add_argument('--min_length', help='min_length', type=int, default=50)
    parser.add_argument('--inter_peak_distance', help='inter_peak_distance', type=int, default=50)    
    parser.add_argument('--merge_close_peaks', help='merge_close_peaks', type=bool, default=True)
    parser.add_argument('--max_length', help='max_length', type=int, default=10000)
    parser.add_argument('--generate_ID', help='generate_ID', type=bool, default=True)

    # Print version
    parser.add_argument('--version', action='version', version='%(prog)s - Version 1.0')

    # Parse arguments
    args = parser.parse_args()
    
    return args

def py_peak_calling(bedgraph, threshold, min_length, inter_peak_distance,
                    merge_close_peaks, max_length, generate_ID):
    '''
    - need to install a more up-to-date varsion of bedtools before invoking Jupyter
    type: module load bedtools/2.21.0
    (1) filters bedgraph based on threshold;
    
    (2) merges adjacent basepairs that are over threshold;
    
    (3) retains peaks that satisfy min/max length criteria; 
    
    (4) merges any peaks that are closer than the inter-peak distance cutoff -or-
    alternatively keeps just the highest peak (this is beta functionality)
    
    - max length is typically defaulted to be very large
    - outputs a bed file (default col4 is the sum of the bedgraph scores;
    sorted by chrom;start;stop)
    - generate ID: will auto generate a integer list as a ID number
    (1... number of peaks). This will be reported as column 4 and the bedgraph
    scores will be shifted to column 5 as per standard bed format
    - note the peak score for merged peak is the *just* the sum of the two 
    individual peaks not the total score in the merged region (i.e. there could
    be some sub-threshold scores in the intervening space that won't be included)
    -assumes bedgraph in standard format <chr> <start> <stop> <score>
    '''
    
    #### Set up the log file and timing
    start = timeit.default_timer()
    
    #### generate name for output ####
    bedgraph_name = glob.glob(bedgraph)
    
    threshold_name = 'thr%s' % (threshold)
    filename = bedgraph_name[0].replace('.bg', '.%s.peaks.bed') % (threshold_name)
    
    print('input bedgraph file: ' + bedgraph_name[0])
    print('output filename: ' + filename)
    
    #### import data as BedTool ####
    data = BedTool(bedgraph)
    
    #### retains intervals above threshold ######
    above_thresh = data.filter(lambda b: float(b.name) >= threshold) 
    
    #### merge adjacent above threshold regions and sum bedgraph scores
    #### (assumes bedgraph score in col 4)
    
    #by increasing d value can allow for 
    merge_regions= above_thresh.merge(d=0, c=4, o='sum' )
    
    #### filter based on length criteria
    peaks = BedTool(merge_regions.filter(lambda x: len(x) >= min_length and len(x) <= max_length))
    
#     print('number of regions identified before merging or filtering: ' + str(peaks.count()))
    
    if merge_close_peaks==True:
        #merge the bonafide peaks if they they are shorter than the inter peak distance and sum scores and sort
        print('merging peaks that are closer than: ' + str(inter_peak_distance))
        merge_peaks = peaks.merge(d=inter_peak_distance, c= 4, o='sum').sort()
   
    print('number of peaks found: ' + str(merge_peaks.count()))
        
    if not generate_ID:
        print('saving sorted peak bed file with no ID')
        
        merge_peaks.saveas(filename)
        
    if generate_ID:
        print('saving sorted peak bed file with ID names')
        
        #change to pandas dataframe
        DF_peaks = merge_peaks.to_dataframe()
        
        #insert new column with id: 1.... # of peaks
        DF_peaks.insert(3, 'id', ['id' + str(item) for item in range(1, (len(DF_peaks)+1))])
        
        ['id' + str(item)  for item in range(1, 5)]
        #save output
        DF_peaks.to_csv(filename, sep = '\t', header = False, index = False)
        
    stop = timeit.default_timer()
    print("Run time:", str(datetime.timedelta(seconds=int(stop - start))))
    
    return 'Finished'

if __name__ == '__main__':
    # Parse the arguments
    args = parseArguments()

    # Raw print arguments
    print("You are running the script with arguments: ")
    for a in args.__dict__:
        print(str(a) + ": " + str(args.__dict__[a]))

    # Run function
    py_peak_calling(args.bedgraph, args.threshold, args.min_length,
                    args.inter_peak_distance, args.merge_close_peaks,
                    args.max_length, args.generate_ID)

