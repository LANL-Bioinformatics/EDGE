Graphic User Interface
######################

The User Interface was mainly implemented in `JQuery Mobile <http://jquerymobile.com>`_ CSS, javascript and perl CGI. It is a HTML5-based user interface system designed to make responsive web sites and apps that are accessible on all smartphone, tablet and desktop devices.

Start GUI
=========

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

The UI composed with Home page, Run EDGE, and Project list in the Left navigation widget, Input and Analyses modules in the Run EDGE page, Job progress widget (right navigation), Action Widget, and Report page.



Initiating an analysis job
==========================


Output path
-----------

Number of CPUs
--------------


Config file
-----------

Batch project submit
--------------------


Choose processes/analyses
=========================

Quality trim and filter
-----------------------

Host removal
------------

Assembly
--------

Community profiling
-------------------


Reference-based Analysis
------------------------

SNP Phylogeny
-------------

Primer validation
-----------------


Primer design
-------------

Blast Contigs
-------------

Job submission
==============

Checking the status of an analysis job
======================================

Monitoring the Resource Usage
=============================

Job Management
==============