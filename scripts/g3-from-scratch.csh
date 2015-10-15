#!/bin/csh

if ($#argv < 2) then
  echo "Usage:  g3-from-scratch.csh  <genome> <tag>"
  echo "           [step<i>  [only]]"
  echo ""
  echo "Run Glimmer3 on the sequence in file <genome> using program"
  echo " long-orfs  to find a training set.  Use <tag> to prefix"
  echo "output files, which are:"
  echo "<tag>.longorfs is coordinate file of training sequences"
  echo "<tag>.train is the multifasta file of training sequences"
  echo "<tag>.icm is the model"
  echo "<tag>.detail is Glimmer3 output detail"
  echo "<tag>.predict is Glimmer3 predictions"
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

# add/change glimmer options here
set glimmeropts = "-o50 -g110 -t30"

set numsteps = 4

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
# Run Glimmer
echo "Step 4 of ${numsteps}:  Running Glimmer3"
$glimmerpath/glimmer3 $glimmeropts $genome $tag.icm $tag
if  ($status != 0)  then
  echo "Failed to run Glimmer3"
  exit
endif
if  ($onestep == "only")  exit


