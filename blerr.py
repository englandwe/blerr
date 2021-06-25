#!/usr/bin/env python

#Takes two bedfiles and intersects each with a third bedfile
#But quickly, because it will break the third bedfile into chunks
#and parallelizes the search

import sys
import os
import math
import operator
import subprocess
import statistics
import argparse
import multiprocessing as mp
import random
from collections import Counter
from datetime import datetime

class BLERRception(Exception):
    pass

def chunkParse(rawval):
    if rawval.isdigit():
        #bytes
        parsed_val=int(rawval)
    elif rawval[0:-1].isdigit():
        if rawval[-1] in ('m','M'):
            #megabytes
            parsed_val=int(rawval[0:-1])*1000*1000
        elif rawval[-1] in ('g','G'):
            #gigabytes
            parsed_val=int(rawval[0:-1])*1000*1000*1000
        else:
            sys.stderr.write('Incorrect chunk size suffix detected!  You entered %s.  Please use M for megabytes, G for gigabytes, or no suffix for bytes.\n' % rawval)
            sys.stderr.write('Stopping...\n')
            sys.exit(1)
    else:
        sys.stderr.write('Incorrect chunk size format detected!  You entered %s.  Value should be an integer, optionally followed by M for megabytes or G for gigabytes.\n' % rawval)
        sys.stderr.write('Stopping...\n')
        sys.exit(1)

    return parsed_val

def chunkify(filename,size):
    file_end = os.path.getsize(filename)
    with open(filename,'rb') as f:
        chunk_end = f.tell()
        while True:
            chunk_start = chunk_end
            f.seek(size,1)
            f.readline()
            chunk_end = f.tell()
            yield chunk_start, chunk_end - chunk_start
            if chunk_end > file_end:
                break

def poolDump(poolname,afile_sorted,featurefile,target_chunksize):
    outjobs=[]
    outjobs=[poolname.apply_async(intFeatureBed, (featurefile,chunk_start,chunk_size,afile_sorted)) for chunk_start,chunk_size in chunkify(featurefile,target_chunksize)]
    return outjobs

def poolFetch(joblist,outfile):
    overlaps=[]
    big_featset=set()
    with open(outfile,'w') as fout:
        for job in joblist:
            bedout,featset=job.get()
            debedout=bedout.decode()
            fout.write(debedout)
            overlaps+=[x.split('\t') for x in debedout.split('\n')]
            big_featset=big_featset | featset
    return overlaps,big_featset

def sortInput(afile,namebreaker):
    #also counts total input ranges
    with open(afile) as infile:
        afile_lines=infile.read().splitlines()
        afile_split=[x.split('\t') for x in afile_lines]
        badlines=[[idx,line] for idx,line in enumerate(afile_split) if len(line) < 3]
        if len(badlines) > 0:
            for badline in badlines:
                sys.stderr.write('ERROR: Nonconforming line found in %s.  Offending line is %s:%s\n' % (afile,badline[0]+1,badline[1]))
            sys.stderr.write('Stopping...\n')
            sys.exit(1)
            #return 1
        else:   
            afile_sort = sorted(afile_split, key = lambda x: (x[0], int(x[1])))
            afile_out='\n'.join(['\t'.join(x) for x in afile_sort])
            outname=afile+'_blerr_sorted_'+namebreaker
            with open(outname,'w') as outfile:
                outfile.write(afile_out+'\n')
            return outname, len(afile_lines)

def intFeatureBed(featurefile,chunkStart,chunkSize,afile):
    featset=set()
    with open(featurefile,'rb') as infile:
        infile.seek(chunkStart)
        bedChunk=infile.read(chunkSize).splitlines()
        chunklines=[x.split(b'\t') for x in bedChunk]
        badlines=[[idx,line] for idx,line in enumerate(chunklines) if len(line) < 4]
        if len(badlines) > 0:
            for badline in badlines:
                sys.stderr.write('ERROR: Nonconforming line found in %s.  Offending line is %s:%s\n' % (featurefile,badline[0]+1,badline[1]))
            sys.stderr.write('Stopping...\n')
            sys.stderr.flush()
            raise BLERRception('ERROR: Nonconforming lines found in feature file.  See above for details.')
        else:
            sortchunk = sorted(chunklines, key = lambda x: (x[0], int(x[1])))
            for line in sortchunk:
                featset.add(line[3].decode())
            outchunk=b'\n'.join([b'\t'.join(x) for x in sortchunk])
            #intersectBed -wa -wb -sorted -a background.bed3 -b $FEATURE_FILENAME > background_overlap.temp
            bedproc = subprocess.Popen(['intersectBed', '-wa', '-wb', '-sorted', '-a', afile, '-b', 'stdin'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            bedout,bederr=bedproc.communicate(outchunk)
            sys.stderr.write('INFO: another chunk bites the dust @ %s' % (str(datetime.now())) + '\n')
            if len(bederr) > 0:
                sys.stderr.write('WARNING: bedtools produced some error output: %s\n' % bederr.decode())
            return bedout,featset

def getCounts(overlaps):
    overlaps_uniq=list(set([(x[0],x[1],x[2],x[7]) for x in overlaps if len(x) > 1]))
    overlap_counts=dict(Counter([x[3] for x in overlaps_uniq]).most_common())
    return overlap_counts

def featureEnrich(back_overlaps,back_total,sub_overlaps,sub_total,feat_list,cutoff):
    enrichscore_dict={}
    for feat in feat_list:
        try:
            back_count=back_overlaps[feat]
        except KeyError:
            back_count=0
        try:
            sub_count=sub_overlaps[feat]
        except KeyError:
            sub_count=0
        if back_count > 0 and sub_count/sub_total >= cutoff:
            enrichscore=math.log2((sub_count/sub_total)/(back_count/back_total))
        else:
            enrichscore=0
        enrichscore_dict[feat]=[back_count,back_total,sub_count,sub_total,enrichscore]

    return enrichscore_dict

def calcZ(enrichscore_dict):
    zscore_dict={}
    allscores=[x[-1] for x in enrichscore_dict.values()]
    scoremean=sum(allscores)/len(allscores)
    scorestdev=statistics.pstdev(allscores)
    for feat,values in enrichscore_dict.items():
        try:
            score=values[-1]
            zscore = (score - scoremean) / scorestdev
        except ZeroDivisionError:
            zscore='NA'
            sys.stderr.write('WARNING: %s had a standard deviation of 0; zscore cannot be calculated.\n' % feat)
        zscore_dict[feat] = values+[zscore]

    return zscore_dict

def writeOutput(outfile_name,zscores):
    with open(outfile_name, 'w') as outfile:
        for key,vals in zscores.items():
            outline='\t'.join([str(x) for x in vals])
            outfile.write('%s\t%s\n' % (key,outline))

if __name__ == '__main__':  

    ap=argparse.ArgumentParser()
    ap.add_argument('-b', '--background', type=str, help='BED file of background ranges, or previously calculated overlap file (see -rb)',required=True)
    ap.add_argument('-s', '--subset', type=str, help='BED file of subset ranges',required=True)
    ap.add_argument('-f', '--feature', type=str, help='BED file of features',required=True)
    ap.add_argument('-rb', '--reuse_background', type=int, default=0, help='Reuse previously calculated background overlap file. Value should be number of entries in original background file. -b should point to the overlap file.')
    ap.add_argument('-rs', '--reuse_subset', type=int, default=0, help='NOT CURRENTLY IMPLEMENTED. Reuse previously calculated subset overlap file. Value shold be number of entries in original subset file. -s should point to the overlap file.')
    ap.add_argument('-o', '--output', type=str, default='blerr_output.tsv',help='Output file name (default=%(default)s)')
    ap.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use (default=%(default)s)')
    ap.add_argument('-c', '--cutoff', type=float, default=0.05, help='Overlap cutoff; minimum ratio of subset feature overlaps to total subset ranges (default=%(default)s)')
    ap.add_argument('-cs', '--chunksize', type=str, default='100M', help='Break feature file into chunks of this size (default=%(default)s). Understands bytes (no suffix), megabytes (suffix M) and gigabytes (suffix G)')
    args = ap.parse_args()

    backfile=args.background
    subfile=args.subset
    featurefile=args.feature
    cores=args.threads
    overlap_cutoff=args.cutoff
    outfile_name=args.output
    raw_chunksize=args.chunksize
    target_chunksize=chunkParse(raw_chunksize)
    rb=args.reuse_background
    rs=args.reuse_subset
    namebreaker=''.join(random.choices('0123456789',k=6))
    
    if rb == 0:
        backfile_sorted,back_total = sortInput(backfile,namebreaker)
    else:
        back_total=args.reuse_background

    if rs == 0:
        subfile_sorted,sub_total = sortInput(subfile,namebreaker)
    else:
        sub_total_total=args.reuse_subset


    bedpool = mp.Pool(cores)

    if rb == 0:
        backjobs = poolDump(bedpool,backfile_sorted,featurefile,target_chunksize)
        back_overlaps,featset = poolFetch(backjobs,'background_blerr_'+namebreaker+'.bed')
    if rs == 0:
        subjobs = poolDump(bedpool,subfile_sorted,featurefile,target_chunksize)
        sub_overlaps,featset2 = poolFetch(subjobs,'subset_blerr_'+namebreaker+'.bed')

    bedpool.close()
    bedpool.join()

    if rb != 0:
        back_overlaps=[]
        with open(backfile) as infile:
            for line in infile:
                back_overlaps.append(line.strip().split('\t'))
        featset=featset2
    
    back_counts=getCounts(back_overlaps)
    sub_counts=getCounts(sub_overlaps)

    enrichscores=featureEnrich(back_counts,back_total,sub_counts,sub_total,featset,overlap_cutoff)
    zscores=calcZ(enrichscores)

    writeOutput(outfile_name,zscores)
   

