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
from collections import Counter
from multiprocessing import Pool
from datetime import datetime

#Break file into chunks to distribute to the pool
def chunkify(fname,size):
    fileEnd = os.path.getsize(fname)
    with open(fname,'rb') as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size,1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break

def sortInput(afile):
    #also counts total input ranges
    with open(afile) as infile:
        afile_lines=infile.read().splitlines()
        afile_split=[x.split('\t') for x in afile_lines]
        afile_sort = sorted(afile_split, key = lambda x: (x[0], int(x[1])))
        afile_out='\n'.join(['\t'.join(x) for x in afile_sort])
        outname=afile+'_blerr_sorted'
        with open(outname,'w') as outfile:
            outfile.write(afile_out+'\n')
    return outname, len(afile_lines)

def intFeatureBed(featurefile,chunkStart,chunkSize,afile):
    featset=set()
    with open(featurefile,'rb') as infile:
        infile.seek(chunkStart)
        bedChunk=infile.read(chunkSize).splitlines()
        chunklines=[x.split(b'\t') for x in bedChunk]
        #sortchunk = sorted(chunklines, key = operator.itemgetter(0, 1))
        sortchunk = sorted(chunklines, key = lambda x: (x[0], int(x[1])))
        #[3] is feature name
        for line in sortchunk:
            featset.add(line[3].decode())
        outchunk=b'\n'.join([b'\t'.join(x) for x in sortchunk])
        #intersectBed -wa -wb -sorted -a background.bed3 -b $FEATURE_FILENAME > background_overlap.temp
        bedproc = subprocess.Popen(['intersectBed', '-wa', '-wb', '-sorted', '-a', afile, '-b', 'stdin'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bedout,bederr=bedproc.communicate(outchunk)
        sys.stderr.write('bederr: %s\n' % bederr.decode())
    sys.stderr.write('BLERRP: another chunk bites the dust @ %s' % (str(datetime.now())) + '\n')
    return bedout,featset


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
    print(allscores)
    scoremean=sum(allscores)/len(allscores)
    scorestdev=statistics.pstdev(allscores)
    for feat,values in enrichscore_dict.items():
        score=values[-1]
        zscore = (score - scoremean) / scorestdev
        zscore_dict[feat] = values+[zscore]
    return zscore_dict
  
backfile=sys.argv[1]
subfile=sys.argv[2]
featurefile=sys.argv[3]
cores=int(sys.argv[4])
#overlap_cutoff=float(sys.argv[5])
overlap_cutoff=0.001

#sort shit
backfile_sorted,back_total = sortInput(backfile)
subfile_sorted,sub_total = sortInput(subfile)

#run the background intersections
bedpool = Pool(cores)

backjobs = []
backjobs=[bedpool.apply_async(intFeatureBed, (featurefile,chunkStart,chunkSize,backfile_sorted)) for chunkStart,chunkSize in chunkify(featurefile,1024*1024*1024)]

subjobs = []
subjobs=[bedpool.apply_async(intFeatureBed, (featurefile,chunkStart,chunkSize,subfile_sorted)) for chunkStart,chunkSize in chunkify(featurefile,1024*1024*1024)]

#assuming everything in the background and subset files is one class of thing
back_overlaps=[]
big_featset=set()
with open('background_blerr.bed','w') as fout:
    for backjob in backjobs:
        bedout,featset=backjob.get()
        debedout=bedout.decode()
        fout.write(debedout)
        back_overlaps+=[x.split('\t') for x in debedout.split('\n')]
        big_featset=big_featset | featset

sub_overlaps=[]
with open('subset_blerr.bed','w') as fout:
    for subjob in subjobs:
        #featlist is redundant, should probably fix this
        bedout,featlist2=subjob.get()
        debedout=bedout.decode()
        fout.write(debedout)
        sub_overlaps+=[x.split('\t') for x in debedout.split('\n')]

bedpool.close()
bedpool.join()

#print(back_overlaps)
#count occurrences of feature overlaps
back_uniq=list(set([(x[0],x[1],x[2],x[7]) for x in back_overlaps if len(x) > 1]))
back_counts=dict(Counter([x[3] for x in back_uniq]).most_common())
sub_uniq=list(set([(x[0],x[1],x[2],x[7]) for x in sub_overlaps if len(x) > 1]))
sub_counts=dict(Counter([x[3] for x in sub_uniq]).most_common())
#back_counts=dict(Counter([x[7] for x in back_overlaps if len(x) > 1]).most_common())
#sub_counts=dict(Counter([x[7] for x in sub_overlaps if len(x) > 1]).most_common())
enrichscores=featureEnrich(back_counts,back_total,sub_counts,sub_total,big_featset,overlap_cutoff)
zscores=calcZ(enrichscores)
with open('blerr_testout', 'w') as outfile:
    for key,vals in zscores.items():
        outline='\t'.join([str(x) for x in vals])
        outfile.write('%s\t%s\n' % (key,outline))
