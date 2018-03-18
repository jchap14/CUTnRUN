##### Import modules     
import pybedtools 
from pybedtools import BedTool
import argparse
import glob
import os
from subprocess import check_output
import pandas as pd
from pybedtools.contrib.bigwig import bedgraph_to_bigwig
import matplotlib.pyplot as plt
plt.switch_backend('agg')

##### Define arguments locally to test script (comment out when bash submitting)
# bamfile="GFP_K27ac.nmSort.bam"
# spikefile="GFP_K27ac.yst.bam"
# lengths_analysis=True
# lengths_image=True
# size_select_1=True
# size_select_2=True
# bedgraph=True
# chrom_sizes='hg19'
# big_wig=True
# size_min_1=20
# size_max_1=120
# size_min_2=150
# size_max_2=710
# multiplying_factor=10000

##### Define parseArguments to use positional variables in bash submission
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('--bamfile', help='bamfile', type=str)

    # Optional arguments with defaults
    parser.add_argument('--spikefile', help='spikefile', type=str)
    parser.add_argument('--lengths_analysis', help='lengths_analysis', type=bool, default=True)
    parser.add_argument('--lengths_image', help='lengths_image', type=bool, default=True)
    parser.add_argument('--size_select_1', help='size_select_1', type=bool, default=True)
    parser.add_argument('--size_select_2', help='size_select_2', type=bool, default=True)
    parser.add_argument('--bedgraph', help='bedgraph', type=bool, default=True)
    parser.add_argument('--chrom_sizes', help='chrom_sizes', type=str, default='hg19')
    parser.add_argument('--big_wig', help='big_wig', type=bool, default=True)
    parser.add_argument('--size_min_1', help='size_min_1', type=int, default=20)
    parser.add_argument('--size_max_1', help='size_max_1', type=int, default=120)
    parser.add_argument('--size_min_2', help='size_min_2', type=int, default=150)
    parser.add_argument('--size_max_2', help='size_max_2', type=int, default=710)
    parser.add_argument('--multiplying_factor', help='multiplying_factor', type=int, default=10000)

    # Print version
    parser.add_argument('--version', action='version', version='%(prog)s - Version 1.0')

    # Parse arguments
    args = parser.parse_args()
    
    return args

##### Define the main function for bdg conversion

def BAM2BDG(bamfile, spikefile, lengths_analysis,
            lengths_image, size_select_1, size_select_2, bedgraph, chrom_sizes,
            big_wig, size_min_1, size_max_1, size_min_2, size_max_2, multiplying_factor):
    
    from datetime import datetime
    startTime = datetime.now()
        
    #generate list of names for data_files and spike_files from the Folder directory
    
    dataFiles = list()
    dataFiles.append(bamfile)
    
    spikeFiles = list()
    spikeFiles.append(spikefile)
    
    #check for equal numbers and matched data_files and spike_files. If incorrect, exit script
    
    if len(dataFiles)==0:
        return 'No datafiles imported. Exiting...'
    
    if len(dataFiles) != len(spikeFiles):
        return 'Unequal numbers of data files and spike files. Exiting...'
    
    for i in range(len(dataFiles)):
        if dataFiles[i].split('.')[0] != spikeFiles[i].split('.')[0]:
            print('Unmatched pairs and spike files.')
            print(dataFiles[i] + ' not paired with ' + spikeFiles[i])
            return 'Exiting...'
    
    print('Properly matched data_files and spike_files. Continuing...')
    
    ####################################################
    
    ## generate bed file from the bam files (this assumes the bam files are just properly patched pairs)
    
    print('Generating bed files representing whole insert from paired end reads in the data files')
    print('\n')
    
    if size_select_1 or size_select_2:
        print('Generating size selected bed files')
        print('\n')
    
    
    ##generate bed file names from data file names (even if ends/size selection set to false)
    
    bed_names = [f.replace('bam', 'bed') for f in dataFiles]
    
    size_selected_files_1 = [f.replace('bam', str(size_min_1) + '_' + str(size_max_1) + '.bed') for f in dataFiles]
    
    size_selected_files_2 = [f.replace('bam', str(size_min_2) + '_' + str(size_max_2) + '.bed') for f in dataFiles]
    
    all_beds = bed_names + size_selected_files_1 + size_selected_files_2
    
    
    #generate filenames for length analysis as will perform on each datafile on fly rather than reloading
    lengths_names = [f.replace('bam', 'lengths') for f in dataFiles]
    #create empty dataframe to be filled by each length analysis, used to plot lengths distribution 
    lengths_df = pd.DataFrame()

    #####################################################

    #generate bed files using bam_to_bed tool (makes bed12 format)
    for i in range(len(dataFiles)):
        temp_bed = BedTool(dataFiles[i]).bam_to_bed(bedpe=True).to_dataframe()
        
        #need to strip out the start and end position of the whole insert (bed12 is both sequenced reads)
        #note column names actually represent <chrom> <start of insert> <end of insert>
        temp_bed_stripped = temp_bed.iloc[:, [0,1,5]].sort_values(by = ['chrom', 'start', 'strand'])
        
        #calculate insert size and insert as column 4 and save file with bed_name
        #these bed files represent the entire insert
        temp_bed_stripped['length'] = temp_bed_stripped['strand'] - temp_bed_stripped['start']
        
        temp_bed_stripped.to_csv(bed_names[i], sep="\t", header=False, index=False)
        
        #perform analysis on the length of inserts sequenced
        if lengths_analysis:
            temp_lengths = temp_bed_stripped.groupby(by=['length'])['length'].count()
            
            temp_lengths.to_csv(lengths_names[i], sep='\t', header = [bed_names[i]], index=True, index_label='length')
            
            #add the lengths data from this datafile to the lengths_df dataframe, title each series with bed file name
            lengths_df = lengths_df.join(temp_lengths.rename(bed_names[i]), how='outer', lsuffix='_left')
        
        #generate size selected whole insert bed files
        if size_select_1:
            subset_1 = temp_bed_stripped[(temp_bed_stripped.iloc[:,3]>=size_min_1)
                                         & (temp_bed_stripped.iloc[:,3]<=size_max_1)]
            
            subset_1.to_csv(size_selected_files_1[i], sep="\t", header=False, index=False)
        
        if size_select_2:
            subset_2 = temp_bed_stripped[(temp_bed_stripped.iloc[:,3]>=size_min_2)
                                         & (temp_bed_stripped.iloc[:,3]<=size_max_2)]
            
            subset_2.to_csv(size_selected_files_2[i], sep="\t", header=False, index=False)
    
    print('finished generating bed files:')
    print('\n')
    print('whole insert bed files:'+'\n'+'\n'.join(bed_names))
    print('\n')
    if size_select_1:
        print('whole insert bed files with size selection #1:'+'\n'+'\n'.join(size_selected_files_1))
        print('\n')
    
    
    if size_select_2:
        print('whole insert bed files with size selection #2:'+'\n'+'\n'.join(size_selected_files_2))
        print('\n')
    
    
    #####################################################
    #take all the bed files and generate spike normalized bedgraph files
    
    if bedgraph:
        print('Generating spike normalized bedgraphs from all the bed files')
        print('\n')
        
        #generate file names for the bedgraphs
        bdg_names  = [f.replace('bed', 'bdg') for f in bed_names]
        
        
        size_selected_files_1_bdg = [f.replace('bed', 'bdg') for f in size_selected_files_1]
        
        size_selected_files_2_bdg = [f.replace('bed', 'bdg') for f in size_selected_files_2]
        
        all_bdg = [f.replace('bed', 'bdg') for f in all_beds]
        
        #calculate the number of reads in each of the spike sam files
        #spikeCount is a list of reads populated by the 'samtools view -c' shell command for each spikeFile
        spikeCount = []
        
        spike_string = []
        
        for item in spikeFiles:
            spike_string.append('samtools view -c -S ' + item)
        
        for item in spike_string:
            spikeCount.append(int(check_output(item, shell =True))/2)
        
        print('spike counts:')
        print(spikeCount)
        
        #calculating list of scaling factors
        scaling_factor = []
        
        for item in spikeCount:
            scaling_factor.append(float(10000)/item)
        
        #run bedtools genomecov to generate bedgraph files
        for i in range(len(bdg_names)):
            BedTool(bed_names[i]).genome_coverage(bg = True, genome = chrom_sizes, scale = scaling_factor[i]).moveto(bdg_names[i])
        
        if size_select_1:
            for i in range(len(size_selected_files_1_bdg)):
                BedTool(size_selected_files_1[i]).genome_coverage(bg = True, genome = chrom_sizes, scale = scaling_factor[i]).moveto(size_selected_files_1_bdg[i])
        
        if size_select_2:
            for i in range(len(size_selected_files_2_bdg)):
                BedTool(size_selected_files_2[i]).genome_coverage(bg = True, genome = chrom_sizes, scale = scaling_factor[i]).moveto(size_selected_files_2_bdg[i])
        
        print('finished generating bedgraph files:')
        print('\n')
        print('whole insert bedgraph files:'+'\n'+'\n'.join(bdg_names))
        print('\n')
        if size_select_1:
            print('whole insert bedgraph files with size selection #1:'+'\n'+'\n'.join(size_selected_files_1_bdg))
            print('\n')
        if size_select_2:
            print('whole insert bedgraph files with size selection #2:'+'\n'+'\n'.join(size_selected_files_2_bdg))
            print('\n')

    #####################################################
    #make bigwig files from all the bedgraphs generated
    
    if big_wig:
        print('Generating big_wig files from each of the bedgraphs')
        
        if big_wig==True and bedgraph==False:
            return 'WARNING: no bedgraphs to make into big_wig files'
        
        #generate file names for the bigwigs
        bw_names  = [f.replace('bdg', 'bw') for f in bdg_names]
        
        size_selected_files_1_bw = [f.replace('bdg', 'bw') for f in size_selected_files_1_bdg]
        
        size_selected_files_2_bw = [f.replace('bdg', 'bw') for f in size_selected_files_2_bdg]
        
        all_bw = [f.replace('bdg', 'bw') for f in all_beds]
        
        #run bedgraph_to_bigwig tool
        for i in range(len(bdg_names)):
            bedgraph_to_bigwig(BedTool(bdg_names[i]), chrom_sizes, bw_names[i])
        
        if size_select_1:
            for i in range(len(size_selected_files_1_bdg)):
                bedgraph_to_bigwig(BedTool(size_selected_files_1_bdg[i]), chrom_sizes, size_selected_files_1_bw[i])
        
        if size_select_2:
            for i in range(len(size_selected_files_2_bdg)):
                bedgraph_to_bigwig(BedTool(size_selected_files_2_bdg[i]), chrom_sizes, size_selected_files_2_bw[i])
        
        print('finished generating bigwig files:')
        print('\n')
        print('whole insert bigwig files:'+'\n'+'\n'.join(bw_names))
        print('\n')
        if size_select_1:
            print('whole insert bigwig files with size selection #1:'+'\n'+'\n'.join(size_selected_files_1_bw))
            print('\n')
        if size_select_2:
            print('whole insert bigwig files with size selection #2:'+'\n'+'\n'.join(size_selected_files_2_bw))
            print('\n')
    
    if lengths_image:
        if lengths_analysis==False:
            print('lengths analysis set to false, so no image to display')
        
        else:
            print('saving combined lengths distribution to file: ' + os.getcwd().split('/')[-1] + str('.lengths'))
            
            temp_name = os.getcwd().split('/')[-1] + str('.lengths')
            lengths_df.to_csv(temp_name, sep='\t', header=True, index=True, index_label='bp')
            
            print('generating image of lengths distribution')
            
            temp_plot_name = temp_name + str('_plot.png')
            
            fig = plt.figure(figsize=(12,6))
            ax  = fig.add_subplot(111)
            ax.set_position([0.1,0.1,0.5,0.8])
            ax.plot(lengths_df)
            leg = ax.legend(lengths_df.columns.values.tolist(), loc = 'center left', bbox_to_anchor = (1.0, 0.5))
            plt.title('Length Distribution')
            plt.xlabel('Insert Lengths (bp)')
            plt.ylabel('Count')
            fig.savefig(temp_plot_name)
    
    print('Runtime (hh:mm:ss): ' + str(datetime.now() - startTime))
    return 'Finished'

if __name__ == '__main__':
    # Parse the arguments
    args = parseArguments()

    # Raw print arguments
    print("You are running the script with arguments: ")
    for a in args.__dict__:
        print(str(a) + ": " + str(args.__dict__[a]))

    # Run function
    BAM2BDG(args.bamfile, args.spikefile,
            args.lengths_analysis, args.lengths_image, args.size_select_1,
            args.size_select_2, args.bedgraph, args.chrom_sizes, args.big_wig, args.size_min_1,
            args.size_max_1, args.size_min_2, args.size_max_2, args.multiplying_factor)
