#!/usr/bin/env python

#Takes two bedfiles and intersects each with a third bedfile
#But quickly, because it will break the third bedfile into chunks
#and parallelizes the search

import sys
import os
import math
import operator
import subprocess
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

def calcEnrich(back_overlap,back_total,sub_overlap,sub_total):
    enrichscore=math.log2((sub_overlap/sub_total)/(back_overlap/back_total))
    return enrichscore

'''
def featureEnrich(back_overlaps,sub_overlaps,feat_list):
    for feat in feat_list:
        #make dict instead?
        feat_back=[x for x in back_overlaps
'''    


backfile=sys.argv[1]
subfile=sys.argv[2]
featurefile=sys.argv[3]
cores=int(sys.argv[4])

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
with open('background_blerr.bed','w') as fout:
    for backjob in backjobs:
        bedout,featlist=backjob.get()
        debedout=bedout.decode()
        fout.write(debedout)
        back_overlaps+=[x.split('\t') for x in debedout.split('\n')]
#print(back_overlaps)
#count occurrences of feature overlaps
#should we be counting unique overlaps instead?  is it even possible for them to be be non-unique?
back_counts=Counter([x[7] for x in back_overlaps if len(x) > 1]).most_common()
print(back_counts)        

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


