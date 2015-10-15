#!/bin/csh

if ($#argv < 3) then
  echo "Usage:  g3-from-training.csh  <genome> <coords> <tag>"
  echo "           [step<i>  [only]]"
  echo ""
  echo "Run Glimmer3 on the sequence in file <genome> using the genes"
  echo "in file <coords> to extract a training set.  Use <tag> to prefix"
  echo "output files, which are:"
  echo "<tag>.train is the multifasta file of training sequences"
  echo "<tag>.icm is the model"
  echo "<tag>.upstream are the regions before the starts in <coords>"
  echo "<tag>.motif is a PWM of the upstream regions"
  echo "<tag>.detail is Glimmer3 output detail"
  echo "<tag>.predict is Glimmer3 predictions"
  echo ""
  echo "If the 6th argument is filled in, then jump to that step"
  echo "If the 7th argument is also set to 'only', then do only that step"

  exit -1;
endif


set genome = $1
set coords = $2
set tag = $3
set step = $4
set onestep = $5

set awkpath = /home/edge/edge/scripts
set glimmerpath = /home/edge/edge/bin
set elphbin = /nfshomes/adelcher/bin/elph

# add/change glimmer options here
set glimmeropts = "-o50 -g110 -t30"

set numsteps = 5

if  ($step != "")  goto $step
    

step1:
# Extract the training sequences from the genome file
echo "Step 1 of ${numsteps}:  Extracting training sequences"
$glimmerpath/extract -t $genome $coords > $tag.train
if  ($status != 0)  then
  echo "Failed to extract training sequences"
  exit
endif
if  ($onestep == "only")  exit


step2:
# Build the icm from the training sequences
echo "Step 2 of ${numsteps}:  Building ICM"
$glimmerpath/build-icm -r $tag.icm < $tag.train
if  ($status != 0)  then
  echo "Failed to build ICM"
  exit
endif
if  ($onestep == "only")  exit


step3:
# Create a position weight matrix (PWM) from the regions
# upstream of the start locations in $coords
echo "Step 3 of ${numsteps}:  Making PWM from upstream regions"
$awkpath/upstream-coords.awk 25 0 $coords \
   | $glimmerpath/extract $genome - > $tag.upstream
$elphbin $tag.upstream LEN=6 | $awkpath/get-motif-counts.awk > $tag.motif
if  ($status != 0)  then
  echo "Failed to create PWM"
  exit
endif
if  ($onestep == "only")  exit


step4:
# Determine the distribution of start-codon usage in $coords
echo "Step 4 of ${numsteps}:  Getting start-codon usage"
set startuse = `$glimmerpath/start-codon-distrib -3 $genome $coords`
if  ($onestep == "only")  exit


step5:
# Run Glimmer
echo "Step 5 of ${numsteps}:  Running Glimmer3"
$glimmerpath/glimmer3 $glimmeropts -b $tag.motif -P $startuse $genome $tag.icm $tag
if  ($status != 0)  then
  echo "Failed to run Glimmer3"
  exit
endif
if  ($onestep == "only")  exit


