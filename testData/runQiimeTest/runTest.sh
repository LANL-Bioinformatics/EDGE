#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi
test_result(){
	Test=$rootdir/TestOutput/QiimeAnalysis/analysis/biom_table_summary.txt
	Expect=$rootdir/biom_table_summary.txt
	testName="EDGE Qiime Analysis pipeline test";
	if cmp -s "$Test" "$Expect"
	then
   		echo "$testName passed!"
   		touch "$rootdir/TestOutput/test.success"
	else
   		echo "$testName failed!"
   		touch "$rootdir/TestOutput/test.fail"
	fi
}

if [ ! -f "$rootdir/TestOutput/test.success" ]
then
        rm -rf $rootdir/TestOutput
fi

cd $rootdir
echo "Working Dir: $rootdir";
echo "EDGE HOME Dir: $EDGE_HOME";

perl $EDGE_HOME/runPipeline -c $rootdir/config.txt -o $rootdir/TestOutput -cpu 4 -noColorLog  -u $rootdir/forward_reads.fastq.gz || true

rm -rf $rootdir/TestOutput/QcReads

test_result;
