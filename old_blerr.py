#!/usr/bin/env python

#Takes two bedfiles and intersects each with a third bedfile
#But quickly, because it will break the third bedfile into chunks
#and parallelizes the search

import sys
import os
import subprocess
from multiprocessing import Pool
from datetime import datetime

#Multiprocessing
def chunkify(fname,size):
    fileEnd = os.path.getsize(fname)
    with open(fname,'r') as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size,1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break

def intFeatureBed(featurefile,chunkStart,chunkSize,afile):
    with open(featurefile) as infile:
        infile.seek(chunkStart)
        #lines = infile.read(chunkSize).splitlines()
        #print(lines)
        bedChunk=infile.read(chunkSize)
        #intersectBed -wa -wb -sorted -a background.bed3 -b $FEATURE_FILENAME > background_overlap.temp
        bedproc = subprocess.Popen(['intersectBed', '-wa', '-wb', '-sorted', '-a', afile, '-b', 'stdin'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bedout,bederr=bedproc.communicate(bedChunk)
        sys.stderr.write('bederr: %s\n' % bederr)
    sys.stderr.write('BLERRP: another chunk bites the dust @ %s' % (str(datetime.now())) + '\n')
    return bedout

#read big bedfile in in chunks

#write the chunk to stdin, feed to bedtools
#bedproc = subprocess.Popen(["grep", "-c", "test"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)


backfile=sys.argv[1]
subfile=sys.argv[2]
featurefile=sys.argv[3]
cores=int(sys.argv[4])

bedpool = Pool(cores)
backjobs = []
backjobs=[bedpool.apply_async(intFeatureBed, (featurefile,chunkStart,chunkSize,backfile)) for chunkStart,chunkSize in chunkify(featurefile,1024*1024*1024)]
#bedoutput=[bedjob.get() for bedjob in bedjobs]

fout=open('background_blerr.bed','w')
for backjob in backjobs:
    fout.write(backjob.get())
fout.close()

subjobs = []
subjobs=[bedpool.apply_async(intFeatureBed, (featurefile,chunkStart,chunkSize,subfile)) for chunkStart,chunkSize in chunkify(featurefile,1024*1024*1024)]

fout2=open('subset_blerr.bed','w')
for subjob in subjobs:
    fout2.write(subjob.get())
fout2.close()

bedpool.close()
bedpool.join()


