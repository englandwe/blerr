# blerr
Identify enriched overlaps between genome features. Fast.

blerr is designed to identify overlaps between genomic regions of interest and a set of genomic features - for example, identifying potential transcription factor binding sites upstream of genes of interest - and find those which occur more or less frequently than expected by chance.

It takes three inputs in BED format: a subset of genomic ranges of interest (`--subset`), the background ranges in your dataset (i.e. the ones not of interest; `--background`), and the features you want to overlap with these two sets of ranges (`--features`).  blerr is designed to efficiently handle very large inputs, such as the entire set of human transcription factor binding motifs.  To enable this, the features file is broken into smaller chunks, and the features search is parallelized.  Previously-computed background overlaps can also be reused to save time.

After identifying overlaps, blerr calculates enrichment and Z-scores to help identify features with significantly more or fewer overlaps with your dataset than expected by chance.

blerr requires only python 3.x.

## Usage

blerr.py [-h] -b BACKGROUND -s SUBSET -f FEATURE [-rb REUSE_BACKGROUND] [-o OUTPUT] [-t THREADS] [-c CUTOFF] [-cs CHUNKSIZE]

  -h, --help
  
  show this help message and exit
  
  -b BACKGROUND, --background BACKGROUND
   
  BED file of background ranges, or previously calculated overlap file (see -rb)
                        
  -s SUBSET, --subset SUBSET

  BED file of subset ranges
                        
  -f FEATURE, --feature FEATURE
  
   BED file of features
                        
  -rb REUSE_BACKGROUND, --reuse_background REUSE_BACKGROUND
  
  Reuse previously calculated background overlap file. Value should be number of entries in original background file. -b should point to the overlap file.
                        
  -o OUTPUT, --output OUTPUT
  
   Output file name (default=blerr_output.tsv)
                        
  -t THREADS, --threads THREADS
  
   Number of threads to use (default=1)
                        
  -c CUTOFF, --cutoff CUTOFF
  
  Overlap cutoff; minimum ratio of subset feature overlaps to total subset ranges (default=0.05)
                        
  -cs CHUNKSIZE, --chunksize CHUNKSIZE
  
  Break feature file into chunks of this size (default=100M). Understands bytes (no suffix), megabytes (suffix M) and gigabytes (suffix G)
