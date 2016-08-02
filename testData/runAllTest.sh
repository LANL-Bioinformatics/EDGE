#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )

if [ -z ${EDGE_HOME+x} ]; then
        EDGE_HOME="$rootdir/.."
fi

cd $rootdir
echo "Working Dir: $rootdir";
echo "EDGE HOME Dir: $EDGE_HOME";
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

TotalTestNUM=0;
for d in $(ls -d */);
do
	if [[ $d =~ run([A-Za-z0-9]+)Test ]]
	then
		TotalTestNUM=$((TotalTestNUM+1))
	fi
done

TestNUM=0;
FailTestNum=0;
SuccessTestNum=0;
for d in $(ls -d */);
do
	if [[ $d =~ run([A-Za-z0-9]+)Test ]]
	then
		TestNUM=$((TestNUM+1))
		TestName=${BASH_REMATCH[1]};
		# Use the bash built in variable SECONDS for the elapsed time since the script invocation.
		printf "[%02d:%02d:%02d]  " "$((SECONDS/3600%24))" "$((SECONDS/60%60))" "$((SECONDS%60))"
		echo -en "[$((TestNUM*100/TotalTestNUM)) %]\t"
		echo -en "Test $TestName ......\t";
		$d/runTest.sh >/dev/null || true
		if [ -f "$d/TestOutput/test.success" ]
		then 
			echo -e "$GREEN[OK]$NC"
			SuccessTestNum=$((SuccessTestNum+1))
		else
			echo -e "$RED[Failed]$NC"
			FailTestNum=$((FailTestNum+1))
		fi
		
	fi
done

echo -e "\n $SuccessTestNum/$TotalTestNUM test(s) passed";

printf "\n Total Running Time: %02d:%02d:%02d\n\n" "$((SECONDS/3600%24))" "$((SECONDS/60%60))" "$((SECONDS%60))"


