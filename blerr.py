#!/usr/bin/env python

#Takes two bedfiles and intersects each with a third bedfile
#But quickly, because it will break the third bedfile into chunks
#and parallelizes the search

import sys
import os
import operator
import subprocess
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
    with open(afile) as infile:
        afile_lines=infile.read().splitlines()
        afile_split=[x.split('\t') for x in afile_lines]
        afile_sort = sorted(afile_split, key = lambda x: (x[0], int(x[1])))
        afile_out='\n'.join(['\t'.join(x) for x in afile_sort])
        outname=afile+'_blerr_sorted'
        with open(outname,'w') as outfile:
            outfile.write(afile_out+'\n')
    return outname


def intFeatureBed(featurefile,chunkStart,chunkSize,afile):
    with open(featurefile,'rb') as infile:
        infile.seek(chunkStart)
        #lines = infile.read(chunkSize).splitlines()
        #print(lines)
        bedChunk=infile.read(chunkSize).splitlines()
        chunklines=[x.split(b'\t') for x in bedChunk]
        #sortchunk = sorted(chunklines, key = operator.itemgetter(0, 1))
        sortchunk = sorted(chunklines, key = lambda x: (x[0], int(x[1])))
        outchunk=b'\n'.join([b'\t'.join(x) for x in sortchunk])
        #bytechunk=bytes(outchunk, 'utf-8')
        #intersectBed -wa -wb -sorted -a background.bed3 -b $FEATURE_FILENAME > background_overlap.temp
        bedproc = subprocess.Popen(['intersectBed', '-wa', '-wb', '-sorted', '-a', afile, '-b', 'stdin'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #bedout,bederr=bedproc.communicate(bytechunk)
        bedout,bederr=bedproc.communicate(outchunk)
        sys.stderr.write('bederr: %s\n' % bederr.decode())
    sys.stderr.write('BLERRP: another chunk bites the dust @ %s' % (str(datetime.now())) + '\n')
    return bedout


#read big bedfile in in chunks

#write the chunk to stdin, feed to bedtools
#bedproc = subprocess.Popen(["grep", "-c", "test"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)


backfile=sys.argv[1]
subfile=sys.argv[2]
featurefile=sys.argv[3]
cores=int(sys.argv[4])

#sort shit
backfile_sorted = sortInput(backfile)
subfile_sorted = sortInput(subfile)


bedpool = Pool(cores)


backjobs = []
backjobs=[bedpool.apply_async(intFeatureBed, (featurefile,chunkStart,chunkSize,backfile_sorted)) for chunkStart,chunkSize in chunkify(featurefile,1024*1024*1024)]
#bedoutput=[bedjob.get() for bedjob in bedjobs]

fout=open('background_blerr.bed','w')
for backjob in backjobs:
    fout.write(backjob.get().decode())
fout.close()

subjobs = []
subjobs=[bedpool.apply_async(intFeatureBed, (featurefile,chunkStart,chunkSize,subfile_sorted)) for chunkStart,chunkSize in chunkify(featurefile,1024*1024*1024)]

fout2=open('subset_blerr.bed','w')
for subjob in subjobs:
    fout2.write(subjob.get().decode())
fout2.close()

''''
subjobs = []
subjobs=[bedpool.apply_async(testBed, (featurefile,chunkStart,chunkSize,subfile)) for chunkStart,chunkSize in chunkify(featurefile,1024*1024*1024)]

fout2=open('subset_blerr.bed','w')
for subjob in subjobs:
   fout2.write(subjob.get().decode())
fout2.close()
'''

bedpool.close()
bedpool.join()


