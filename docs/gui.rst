Graphic User Interface
######################

EDGE can be run by a Graphic User Interface (GUI).

To run gui, type::

    $EDGE_HOME/start_edge_ui.sh

will start a localhost and the GUI html page will be opened by the default Browser. 

.. note:: If desktop environment is available, after installation, a "Start EDGE UI" icon should be on the desktop. Click on the green icon" and choose "Run in Terminal" should be the same as above method to start the GUI.

.. image:: img/edge_desktop_icon.png
.. image:: img/start_ui_in_terminal.png
 
The URL address is 127.0.0.1:8080/index.html. It may be not powerful as it is hosted by Apache HTTP Server but it works. With system administrator help, the Apache HTTP Server is suggested to host the gui interface. 
 
.. note:: You may need to configure the edge_wwwroot and input and output in the edge_ui/edge_config.tmpl file while configuring the Apache HTTP Server and link to external drive or network drive if needed.

A Terminal window will display messages and errors as you run EDGE. Under normal operating conditions you can minimize this window. Should an error/problem arise, you may maximize this window to view the error. 

.. image:: img/Terminal_log.png

.. Warning:: IMPORTANT: Do not close this window!

The Browser window is the window in which you will interact with EDGE.