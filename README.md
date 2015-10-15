# EDGE Bioinformatics

This is docker image documentation for version 1.1 of EDGE Bioinformatics, a product of collaboration between Los Alamos National Laboratory and the Naval Medical Research Center sponsored by the Defense Threat Reduction Agency.

EDGE is a highly adaptable bioinformatics platform that allows laboratories to quickly analyze and interpret genomic sequence data. The bioinformatics platform allows users to address a wide range of use cases including assay validation and the characterization of novel biological threats, clinical samples, and complex environmental samples.

The EDGE docker image is available at https://hub.docker.com/r/chienchilo/bioedge/ 

# How to use this image

## Install Docker

See Docker at https://www.docker.com/

## Obtain the docker image

    $ docker pull chienchilo/bioedge

## Obtain inital mysql database

    $ git clone -b docker https://github.com/LANL-Bioinformatics/EDGE.git

## Download EDGE database from ftp server

    $ ftp://ftp.lanl.gov/public/genome/EDGE/1.1/

## Start EDGE bioinformatics instance

    $ docker run -d --cap-add SYS_PTRACE -v /path/to/mysql:/var/lib/mysql -v /path/to/database:/opt/apps/edge/database -v /path/to/EDGE_output:/opt/apps/edge/edge_ui/EDGE_output -v /path/to/EDGE_input:/opt/apps/edge/edge_ui/EDGE_input -p 80:80 -p 8080:8080 --name edge chienchilo/bioedge
    
Wait for few seconds for the docker image to start EDGE service and Open http://localhost/ on the browser to start experience EDGE.

* The -v /path/to/mysql:/var/lib/mysql part of the command mounts the /my/own/mysql (obtain from the git clone above) directory from the underlying host system as /var/lib/mysql inside the container, where MySQL by default will write its data files. Using this to persist the database data in the host.
* The -v /path/to/database:/opt/apps/edge/database mounts the databse obtained from the above download step. 
* The -v /path/to/EDGE_input://opt/apps/edge/edge_ui/EDGE_input mounts the EDGE input directory structure (obtain from the git clone above) to persist the input/upload files/user projects in the host. 
* The -v /path/to/EDGE_output://opt/apps/edge/edge_ui/EDGE_output mounts the EDGE output directory to persist the output files in the host. 
* The -p host:container bind the host port 80 and 8080 to container port 80 and 8080 inside the container. You can change the 80 and 8080 to fit your host system requirements.

## Note

* This image is built on top of offical Ubuntu 14.04.3 LTS and is officially supported on Docker version 1.8.2.
* A local web service (e.g. Apache) is required for viewing http://localhost/
* There is an issue to mount host volume using MAC OSX. https://github.com/boot2docker/boot2docker/issues/581.


## Contact Info
Chien-Chi Lo: <chienchi@lanl.gov>  
Paul Li: <po-e@lanl.gov>  
Anderson, Joseph J. CIV: <Joseph.Anderson@dtra.mil>

