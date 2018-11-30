#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi
test_result(){
	MainErrLog=$rootdir/TestOutput/error.log
	TestLog=$rootdir/TestOutput/ReadsBasedAnalysis/SpecialtyGenes/log.txt
	Test=$rootdir/TestOutput/AssemblyBasedAnalysis/SpecialtyGenes/testSpecialtyGenesProfiling_AR_genes_rgi.txt
	Expect=$rootdir/Profiling_AR_genes_rgi.txt
	Expect3=$rootdir/Profiling_AR_genes_rgi2.txt
	Expect4=$rootdir/Profiling_AR_genes_rgi3.txt
	Test2=$rootdir/TestOutput/ReadsBasedAnalysis/SpecialtyGenes/testSpecialtyGenesProfiling_VF_genes_ShortBRED_table.txt
	Expect2=$rootdir/VF_genes_ShortBRED_table.txt
	Expect5=$rootdir/VF_genes_ShortBRED_table2.txt
	testName="EDGE Specialty Genes Profiling test";
	if cmp -s "$Test" "$Expect" || cmp -s "$Test" "$Expect3" || cmp -s "$Test" "$Expect4"
	then
		if cmp -s "$Test2" "$Expect2" || cmp -s "$Test2" "$Expect5"
		then
   			echo "$testName passed!"
 			touch "$rootdir/TestOutput/test.success"
		else
			echo "Virluence Gene profiling test failed"
			if [ -f "$TestLog" ]
			then
				cat $TestLog >> $MainErrLog
			fi
   			touch "$rootdir/TestOutput/test.fail"
		fi
	else
   		echo "$testName failed!"
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

perl $EDGE_HOME/runPipeline -c $rootdir/config.txt -o $rootdir/TestOutput -cpu 4 -noColorLog -p $rootdir/../Ecoli_10x.1.fastq $rootdir/../Ecoli_10x.2.fastq || true

rm -rf $rootdir/TestOutput/QcReads
rm -f $EDGE_HOME/scripts/specialty_genes/mysql.ini

test_result;
