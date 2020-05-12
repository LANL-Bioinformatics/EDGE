#!/bin/bash
# Usage: deinterleave_fastq.sh < interleaved.fastq f.fastq r.fastq [compress]
# 
# Deinterleaves a FASTQ file of paired reads into two FASTQ
# files specified on the command line. Optionally GZip compresses the output
# FASTQ files using pigz if the 3rd command line argument is the word "compress"
# 
# Can deinterleave 100 million paired reads (200 million total
# reads; a 43Gbyte file), in memory (/dev/shm), in 4m15s (255s)
# 
# Latest code: https://gist.github.com/3521724
# Also see my interleaving script: https://gist.github.com/4544979
# 
# Inspired by Torsten Seemann's blog post:
# http://thegenomefactory.blogspot.com.au/2012/05/cool-use-of-unix-paste-with-ngs.html

# Set up some defaults
GZIP_OUTPUT=0

# If the third argument is the word "compress" then we'll compress the output using pigz
if [[ $3 == "compress" ]]; then
  GZIP_OUTPUT=1
fi

if [[ ${GZIP_OUTPUT} == 0 ]]; then
  paste - - - - - - - -  | tee >(cut -f 1-4 | perl -ne 's/(\.\d+)\.[12] /$1 /g; print;' | tr "\t" "\n" | egrep -v '^$' > $1) | cut -f 5-8 |  perl -ne 's/(\.\d+)\.[12] /$1 /g; print;' | tr "\t" "\n" | egrep -v '^$' > $2
else
  paste - - - - - - - -  | tee >(cut -f 1-4 |  perl -ne 's/(\.\d+)\.[12] /$1 /g; print;' |  tr "\t" "\n" | egrep -v '^$' | gzip --best > $1) | cut -f 5-8 |  perl -ne 's/(\.\d+)\.[12] /$1 /g; print;' | tr "\t" "\n" | egrep -v '^$' | gzip --best > $2
fi
