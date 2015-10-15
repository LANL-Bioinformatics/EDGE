#!/bin/csh

if ($#argv < 2) then
  echo "Usage:  g3-iterated.csh  <genome> <tag>"
  echo "           [step<i>  [only]]"
  echo ""
  echo "Run Glimmer3 on the sequence in file <genome> using program"
  echo " long-orfs  to find a training set.  Use those predictions as a"
  echo "training set for a second run of Glimmer3.  Use <tag> to prefix"
  echo "output files, which are:"
  echo "<tag>.longorfs is coordinate file of training sequences"
  echo "<tag>.train is the multifasta file of training sequences"
  echo "<tag>.icm is the model"
  echo "<tag>.run1.detail is the first Glimmer3 output detail"
  echo "<tag>.run1.predict is the first Glimmer3 predictions"
  echo "<tag>.coords are the training coordinates for starts and the PWM"
  echo "<tag>.upstream are the regions before the starts in <tag>.coords"
  echo "<tag>.motif is a PWM of the upstream regions"
  echo "<tag>.detail is the final Glimmer3 output detail"
  echo "<tag>.predict is the final Glimmer3 predictions"
  echo ""
  echo "If the 6th argument is filled in, then jump to that step"
  echo "If the 7th argument is also set to 'only', then do only that step"

  exit -1;
endif


set genome = $1
set tag = $2
set step = $3
set onestep = $4

set awkpath = /home/edge/edge/scripts
set glimmerpath = /home/edge/edge/bin
set elphbin = /nfshomes/adelcher/bin/elph

# add/change glimmer options here
set glimmeropts = "-o50 -g110 -t30"

set numsteps = 8

if  ($step != "")  goto $step
    

step1:
# Find long, non-overlapping orfs to use as a training set
echo "Step 1 of ${numsteps}:  Finding long orfs for training"
$glimmerpath/long-orfs -n -t 1.15 $genome $tag.longorfs
if  ($status != 0)  then
  echo "Failed to find long-orf training set"
  exit
endif
if  ($onestep == "only")  exit


step2:
# Extract the training sequences from the genome file
echo "Step 2 of ${numsteps}:  Extracting training sequences"
$glimmerpath/extract -t $genome $tag.longorfs > $tag.train
if  ($status != 0)  then
  echo "Failed to extract training sequences"
  exit
endif
if  ($onestep == "only")  exit


step3:
# Build the icm from the training sequences
echo "Step 3 of ${numsteps}:  Building ICM"
$glimmerpath/build-icm -r $tag.icm < $tag.train
if  ($status != 0)  then
  echo "Failed to build ICM"
  exit
endif
if  ($onestep == "only")  exit


step4:
# Run first Glimmer
echo "Step 4 of ${numsteps}:  Running first Glimmer3"
$glimmerpath/glimmer3 $glimmeropts $genome $tag.icm $tag.run1
if  ($status != 0)  then
  echo "Failed to run Glimmer3"
  exit
endif
if  ($onestep == "only")  exit


step5:
# Get training coordinates from first predictions
echo "Step 5 of ${numsteps}:  Getting training coordinates"
tail -n +2 $tag.run1.predict > $tag.coords
if  ($status != 0)  then
  echo "Failed to get training coordinates"
  exit
endif
if  ($onestep == "only")  exit


step6:
# Create a position weight matrix (PWM) from the regions
# upstream of the start locations in $tag.coords
echo "Step 6 of ${numsteps}:  Making PWM from upstream regions"
$awkpath/upstream-coords.awk 25 0 $tag.coords \
   | $glimmerpath/extract $genome - > $tag.upstream
$elphbin $tag.upstream LEN=6 | $awkpath/get-motif-counts.awk > $tag.motif
if  ($status != 0)  then
  echo "Failed to create PWM"
  exit
endif
if  ($onestep == "only")  exit


step7:
# Determine the distribution of start-codon usage in $tag.coords
echo "Step 7 of ${numsteps}:  Getting start-codon usage"
set startuse = `$glimmerpath/start-codon-distrib -3 $genome $tag.coords`
if  ($onestep == "only")  exit


step8:
# Run second Glimmer
echo "Step 8 of ${numsteps}:  Running second Glimmer3"
$glimmerpath/glimmer3 $glimmeropts -b $tag.motif -P $startuse $genome $tag.icm $tag
if  ($status != 0)  then
  echo "Failed to run Glimmer3"
  exit
endif
if  ($onestep == "only")  exit


