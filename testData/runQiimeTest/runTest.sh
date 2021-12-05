#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi
test_result(){
	MainErrLog=$rootdir/TestOutput/error.log
	TestLog=$rootdir/TestOutput/QiimeAnalysis/errorLog.txt
	Test=$rootdir/TestOutput/QiimeAnalysis/TaxonomyAnalysis/Table/feature-table-taxanomy.tsv
	Test2=$rootdir/TestOutput/QiimeAnalysis/TaxonomyAnalysis/Table/metadata.tsv
	Test3=$rootdir/TestOutput/QiimeAnalysis/DiversityAnalysis/table.html
	Expect=$rootdir/feature-table-taxanomy.tsv
	testName="EDGE Qiime2 Analysis pipeline test";
	if [ -s $Test3 ] 
	then
		if [ -s $Test2 -a -s $Test ]
		then
   			echo "$testName passed!"
   			touch "$rootdir/TestOutput/test.success"
		else
   			echo "$testName failed on Taxonomy Anaylysis!"
			if [ -f "$TestLog" ]
			then
				cat $TestLog >> $MainErrLog
			fi
   			touch "$rootdir/TestOutput/test.fail"
		fi
	else
		echo "$testName failed! on Diversity Analysis!"
		if [ -f "$TestLog" ]
		then
			cat $TestLog >> $MainErrLog
		fi
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

sed -e "s;=barcodes.fastq.gz;=$PWD/barcodes.fastq.gz;" -e "s;=map.txt;=$PWD/map.txt;" $rootdir/config.txt > $rootdir/config_run.txt
perl $EDGE_HOME/runPipeline -c $rootdir/config_run.txt -o $rootdir/TestOutput -cpu 4 -noColorLog  -u $rootdir/forward_reads.fastq.gz || true

rm -rf $rootdir/TestOutput/QcReads

test_result;
