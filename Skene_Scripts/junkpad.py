bedgraph='GFP__K7ac.sorted.bg'
threshold=0.01
min_length=20
inter_peak_distance=50
merge_close_peaks=True
keep_highest_close_peak=False
max_length=10000
generate_ID=True
output_name=None
delete_overlap_bed=None

import pybedtools
import glob
from pybedtools import BedTool
import pandas as pd
import csv

if merge_close_peaks==keep_highest_close_peak:
    return 'Exiting... merge_close_peaks and keep_highest_close_peak set the same'

#generate name for output
bedgraph_name = sorted(glob.glob(bedgraph))

threshold_name = "FDR%s" % (int(threshold * 100))
filename = bedgraph_name[0].replace('.bg', '.%s.peaks.bed') % (threshold_name)

print('input bedgraph file: ' + bedgraph_name[0])
print('output filename: ' + filename)

#import data as BedTool
data = BedTool(bedgraph) 

#retains intervals above threshold
above_thresh = data.filter(lambda b: float(b.name) >= threshold) 

#merge adjacent above threshold regions and sum bedgraph scores (assumes bedgraph score in col 4)
#by increasing d value can allow for 
merge_regions= above_thresh.merge(d=0, c=4, o='sum' )

#filter based on length criteria
peaks = BedTool(merge_regions.filter(lambda x: len(x) >= min_length and len(x) <= max_length))

#     print('number of regions identified before merging or filtering: ' + str(peaks.count()))

if merge_close_peaks==True:
    #merge the bonafide peaks if they they are shorter than the inter peak distance and sum scores and sort
    print('merging peaks that are closer than: ' + str(inter_peak_distance))
    merge_peaks = peaks.merge(d=inter_peak_distance, c= 4, o='sum').sort()

print('number of peaks found: ' + str(merge_peaks.count()))

if delete_overlap_bed!=None:
    print('delete_overlap_bed provided: ' + delete_overlap_bed)
    merge_peaks = merge_peaks.intersect(b=delete_overlap_bed, v=True)
    print('number of peaks retained: ' + str(merge_peaks.count()))

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
    
return 'Finished'