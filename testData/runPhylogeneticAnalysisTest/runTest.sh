#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi

test_result(){
	MainErrLog=$rootdir/TestOutput/error.log
	TestLog=$rootdir/TestOutput/SNP_Phylogeny/log.txt
	Test=$rootdir/TestOutput/SNP_Phylogeny/testPhylogeneticAnalysis_summaryStatistics.txt
	Expect=$rootdir/testPhylogeneticAnalysis_summaryStatistics.txt
	Expect2=$rootdir/testPhylogeneticAnalysis_summaryStatistics2.txt
	testName="EDGE Phylogenetic Analysis test";
	if cmp -s "$Test" "$Expect"
	then
		echo "$testName passed!"
		touch "$rootdir/TestOutput/test.success"
	else
		if cmp -s "$Test" "$Expect2"
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
	fi
}

cd $rootdir
echo "Working Dir: $rootdir";
echo "EDGE HOME Dir: $EDGE_HOME";

if [ ! -f "$rootdir/TestOutput/test.success" ]
then
	rm -rf $rootdir/TestOutput
fi

perl $EDGE_HOME/runPipeline -c $rootdir/config.txt -o $rootdir/TestOutput -cpu 4 -noColorLog  -p $rootdir/ebola_R1.fastq $rootdir/ebola_R2.fastq || true
rm -rf $rootdir/TestOutput/QcReads

test_result;
