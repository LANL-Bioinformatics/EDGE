#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi

test_result(){
	MainErrLog=$rootdir/TestOutput/error.log
	TestLog=$rootdir/TestOutput/AssemblyBasedAnalysis/assembly.log
	Test=$rootdir/TestOutput/AssemblyBasedAnalysis/testAssembly_contigs.fa
	testName="EDGE Reads Assembly test";
	if [[ $(find $Test -type f -size +4500000c 2>/dev/null) ]]
	then
		grep -c ">" $Test | awk '{print "Contigs number: " $1}'
		echo "$testName finished"
		touch "$rootdir/TestOutput/test.success"
	else
		echo "$testName failed!"
		if [ -f "$TestLog" ]
		then
			cat $TestLog >>  $MainErrLog
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
