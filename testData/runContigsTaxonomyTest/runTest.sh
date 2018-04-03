#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi

test_result(){
	MainErrLog=$rootdir/TestOutput/error.log
	TestLog=$rootdir/TestOutput/AssemblyBasedAnalysis/Taxonomy/log.txt
	Test=$rootdir/TestOutput/AssemblyBasedAnalysis/Taxonomy/summary_by_topHitCount.txt
	Expect=$rootdir/summary_by_topHitCount.txt
	Expect2=$rootdir/summary_by_topHitCount2.txt
	testName="EDGE Contigs Taxonomy test";
	if cmp -s "$Test" "$Expect" || cmp -s "$Test" "$Expect2"
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
