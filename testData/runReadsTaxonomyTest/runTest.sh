#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi

test_result(){
	MainErrLog=$rootdir/TestOutput/error.log
	TestLog=$rootdir/TestOutput/ReadsBasedAnalysis/Taxonomy/error.log
	Test=$rootdir/TestOutput/ReadsBasedAnalysis/Taxonomy/report/summary.txt
	Expect=$rootdir/summary.txt
	Expect2=$rootdir/summary2.txt
	Expect3=$rootdir/summary3.txt
	Expect4=$rootdir/summary4.txt
	Expect5=$rootdir/summary5.txt
	Expect6=$rootdir/summary6.txt
	testName="EDGE Reads Taxonomy test";
	if cmp -s "$Test" "$Expect" || cmp -s "$Test" "$Expect2" || cmp -s "$Test" "$Expect3" || cmp -s "$Test" "$Expect4" || cmp -s "$Test" "$Expect5" || cmp -s "$Test" "$Expect6"
	then
		echo "$testName passed!"
		touch "$rootdir/TestOutput/test.success"
	else
		echo "$testName failed!"
		if [ -f "$TestLog" ]
		then
			cat $TestLog >> $MainErrLog
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
