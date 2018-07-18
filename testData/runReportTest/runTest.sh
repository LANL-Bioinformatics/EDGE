#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
	EDGE_HOME="$rootdir/../../"
fi

test_result(){
	failTestNum=0
	TotalTestFileNum=6
  	testName="EDGE Report test";
	final_PDF=$rootdir/TestOutput/final_report.pdf
	SVGFile=$rootdir/TestOutput/ReadsBasedAnalysis/Taxonomy/report/1_allReads/metaphlan2/allReads-metaphlan2.tree.svg
	PDFFile=$rootdir/TestOutput/ReadsBasedAnalysis/Taxonomy/report/1_allReads/metaphlan2/allReads-metaphlan2.tree.pdf
	PNGFile=$rootdir/TestOutput/HTML_Report/images/allReads-metaphlan2.tree.png
	HEATMAPPDFFile=$rootdir/TestOutput/ReadsBasedAnalysis/Taxonomy/report/heatmap_TOOL-metaphlan2.species.pdf
	HTMLFile=$rootdir/TestOutput/HTML_Report/report.html
	HTMLlog=$rootdir/TestOutput/HTML_Report/log.txt
	
	if [ -s "$final_PDF" ]
	then
		echo "Final PDF file generated successfully";
	else
		echo "Final PDF file generated failed";
		failTestNum=$((failTestNum+1))
	fi

	# $EDGE_HOME/scripts/microbial_profiling/script/phylo_dot_plot.pl
	if [ -s "$SVGFile" ]
	then
		echo "Tree SVG file generated successfully";
	else
		echo "Tree SVG file generated failed";
		failTestNum=$((failTestNum+1))
	fi
	# $EDGE_HOME/scripts/svg2pdf  (dependency: inkscape)
	if [ -s "$PDFFile" ]
	then
		echo "SVG to PDF file conversion successfully";
	else
		echo "SVG to PDF file conversion failed";
		failTestNum=$((failTestNum+1))
	fi
	# Dependency: ImageMagic (convert)
	if [ -s "$PNGFile" ]
	then
		echo "PDF to PNG file conversion successfully";
	else
		echo "PDF to PNG file conversion failed";
		failTestNum=$((failTestNum+1))
	fi
	# $EDGE_HOME/scripts/microbial_profiling/script/heatmap_distinctZ_noClust_zeroRowAllow.py
	if [ -s "$HEATMAPPDFFile" ]
	then
		echo "HEATMAP PDF file generated successfully";
	else
		echo "HEATMAP PDF file generated failed";
		failTestNum=$((failTestNum+1))
	fi
	# $EDGE_HOME/scripts/outputMunger_w_temp.pl
	if [ -s "$HTMLFile" ] && [ ! -s "$HTMLlog" ]
	then
		echo "HTML report Success"
	else
		echo "HTML report failed"
		failTestNum=$((failTestNum+1))
	fi

	if [ $failTestNum -gt 0 ] 
	then
                echo -e "\n$testName failed!"
                touch "$rootdir/TestOutput/test.fail"
	else
                echo -e "\n$testName passed!"
                touch "$rootdir/TestOutput/test.success"
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
