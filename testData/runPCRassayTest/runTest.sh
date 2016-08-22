#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi

test_result(){
	MainErrLog=$rootdir/TestOutput/error.log
	Test=$rootdir/TestOutput/AssayCheck/pcrContigValidation.txt
	Expect=$rootdir/pcrContigValidation.txt
	Test2=$rootdir/TestOutput/AssayCheck/PCR.design.primers.gff3
	testName="EDGE PCR Assay Analysis test";
	Test2Log=$rootdir/TestOutput/AssayCheck/PCR.design.log
	if cmp -s "$Test" "$Expect"
	then
		if [ -s "$Test2" ]
		then
			echo "$testName passed!"
			touch "$rootdir/TestOutput/test.success"
		else
			echo "PCA design failed!"
			if [ -f "$Test2Log" ] 
			then
				cat $Test2Log >> $MainErrLog
			fi
			touch "$rootdir/TestOutput/test.fail"
		fi
	else
		echo "PCR validation failed!"
		if [ -f "$Test" ] 
		then
			cat $Test >> $MainErrLog
		fi
		touch "$rootdir/TestOutput/test.fail"
	fi
}

cd $rootdir
echo "Working Dir: $rootdir";
echo "EDGE HOME Dir: $EDGE_HOME";

if [ ! -f "$rootdir/TestOutput/test.success" ]
then
	rm -rf $rootdir/TestOutput
fi

perl $EDGE_HOME/runPipeline -c $rootdir/config.txt -o $rootdir/TestOutput -cpu 4 -noColorLog  -p $rootdir/../Ecoli_10x.1.fastq $rootdir/../Ecoli_10x.2.fastq || true

rm -rf $rootdir/TestOutput/QcReads

test_result;
