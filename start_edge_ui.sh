#!/usr/bin/env bash

sleep 1
clear

EDGE_HOME=$( cd $(dirname $0) ; pwd -P )

function checkSystemInstallation
{
    IFS=:
    for d in $PATH; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}

function open_index
{
     if [[ "$OSTYPE" == "darwin"* ]]
     then
        open http://127.0.0.1$1/index.html
     elif ( checkSystemInstallation x-www-browser )
     then
         x-www-browser http://127.0.0.1$1/index.html
     elif ( checkSystemInstallation firefox )
     then
         firefox http://127.0.0.1$1/index.html
     elif ( checkSystemInstallation google-chrome )
     then
         google-chrome http://127.0.0.1$1/index.html
     else
         echo ""
     fi
}

function check {
     if [ $? -ne 0 ] ; then
         curl -s -f -o "/dev/null" --noproxy 127.0.0.1 http://127.0.0.1:8080/index.html
         if [ $? -ne 0 ] ; then
             echo "Turning on localhost"
             python $EDGE_HOME/scripts/httpserver.py --webdir $EDGE_HOME/edge_ui --port 8080 &
             sleep 3
             open_index :8080
             $SHELL
         else
             open_index :8080
             $SHELL
         fi
     else
       open_index
       $SHELL
     fi

}

curl -s -f -o "/dev/null" --noproxy 127.0.0.1 http://127.0.0.1/index.html
check;
