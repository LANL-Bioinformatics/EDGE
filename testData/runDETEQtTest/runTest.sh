#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi

test_result(){
	Test=$rootdir/TestOutput/DETEQT/stats/testDETEQT.report.txt 
	Expect=$rootdir/report.txt
	testName="EDGE DETEQT test";
	if cmp -s "$Test" "$Expect"
	then
		echo "$testName passed!"
		touch "$rootdir/TestOutput/test.success"
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

sed -e "s;=\.;=$PWD/;" -e "s;=sample_test.txt;=$PWD/sample_test.txt;" -e "s;=targeted_reference.fa;=$PWD/targeted_reference.fa;" $rootdir/config.txt > $rootdir/config_run.txt
perl $EDGE_HOME/runPipeline -c $rootdir/config_run.txt -o $rootdir/TestOutput -cpu 4 -noColorLog || true


test_result;
