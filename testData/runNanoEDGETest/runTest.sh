#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

#if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
#fi

test_result(){
	MainErrLog=$rootdir/TestOutput/error.log
	Test=$rootdir/TestOutput/ReadsBasedAnalysis/readsMappingToRef/readsToRef.alnstats.txt
	Test2=$rootdir/TestOutput/AssemblyBasedAnalysis/readsMappingToContig/readsToContigs.alnstats.txt
	Test3=$rootdir/TestOutput/QcReads/NanoStats.txt
	Expect=$rootdir/readsToRef.alnstats.txt
	Expect2=$rootdir/readsToContigs.alnstats.txt
	Expect3=$rootdir/NanoStats.txt
	Expect4=$rootdir/readsToRef.alnstats2.txt
	testName="EDGE Nanopore data analysis test";
	if cmp -s "$Test" "$Expect" || cmp -s "$Test" "$Expect4"
	then
		if cmp -s "$Test2" "$Expect2" && cmp -s "$Test3" "$Expect3"
		then
			echo "$testName passed!"
			touch "$rootdir/TestOutput/test.success"
		else
			echo "$testName failed!"
			touch "$rootdir/TestOutput/test.fail"
		fi
	else
		echo "$testName failed!"
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

perl $EDGE_HOME/runPipeline -c $rootdir/config.txt -o $rootdir/TestOutput -cpu 4 -noColorLog  -u lambda.fastq.gz || true

#rm -rf $rootdir/TestOutput/QcReads
test_result;
