#!/bin/bash

EDGE_HOME=$( cd $(dirname $0)/.. ; pwd -P )
BOKEH_URL="http://localhost:5100/pangia-vis?r=pangia-vis/data/test.tsv"
BOKEH_LOG=$EDGE_HOME/edge_ui/bokeh.log
minimumsize=10240

# check if bokeh is still healthy
curl $BOKEH_URL &> /dev/null
EXITCODE=$?

if [[ $(cat $BOKEH_LOG | grep -E 'NoneType|No active exception') ]]; then
	DATETIME=`date "+%x %H:%M:%S"`
	echo "[$DATETIME] [WARNING] Stale Bokeh server detected. Restarting..."

	# stop current bokeh
	if [[ $(ps auxww | grep 'serve pangia-vis' | awk '$11!~/^grep/ {print $2}') ]]; then
 		$SUDO kill -9 $(ps aux | grep 'serve pangia-vis' | awk '$11!~/^grep/ {print $2}');
	fi
        
	# start bokeh server
	cd ${EDGE_HOME}/thirdParty/pangia
	nohup ${EDGE_HOME}/thirdParty/Anaconda3/bin/bokeh serve pangia-vis --port 5100 --allow-websocket-origin '*' &> ${EDGE_HOME}/edge_ui/bokeh.log &
	#nohup ${EDGE_HOME}/thirdParty/Anaconda3/bin/bokeh serve pangia-vis --port 5100 --keep-alive 10000 --allow-websocket-origin '*' --use-xheaders &> ${EDGE_HOME}/edge_ui/bokeh.log &
else
	DATETIME=`date "+%x %H:%M:%S"`
	if [[ $EXITCODE -ne 0 ]]; then
		echo "[$DATETIME] [INFO] Bokeh server is OFFLINE."
		# start bokeh server	
		cd ${EDGE_HOME}/thirdParty/pangia
		nohup ${EDGE_HOME}/thirdParty/Anaconda3/bin/bokeh serve pangia-vis --port 5100 --allow-websocket-origin '*' &> ${EDGE_HOME}/edge_ui/bokeh.log &
		# --use-xheaders for https
		#nohup ${EDGE_HOME}/thirdParty/Anaconda3/bin/bokeh serve pangia-vis --port 5100 --keep-alive 10000 --allow-websocket-origin '*' --use-xheaders &> ${EDGE_HOME}/edge_ui/bokeh.log &
	else
		echo "[$DATETIME] [INFO] Bokeh server seems ONLINE."
	fi
fi

if [ -f $BOKEH_LOG ];then
	actualsize=$(du -k "$BOKEH_LOG" | cut -f 1)
	if [ $actualsize -ge $minimumsize ]; then
		echo > $BOKEH_LOG
	fi
fi

