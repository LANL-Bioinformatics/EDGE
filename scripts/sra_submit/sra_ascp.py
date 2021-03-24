#!/usr/bin/env python3

__email__ = "conrad.shyu@nih.gov"
__author__ = "Conrad Shyu"
__version__ = "2.1"
__branch__ = "Bioinformatics and Computational Biosciences Branch"
__company__ = "National Institute of Allergy and Infectious Diseases"
__address__ = "5601 Fishers Lane, Rockville, MD 20852"
__update__ = "10/4/2019"
__project__ = "METAGENOTE"

import os
import sys
import argparse

def main(argv):
    """
    send all files to NCBI
    """
    dst_dir = "%s@%s:%s" % (argv.ncbi_user, "upload.ncbi.nlm.nih.gov", argv.ncbi_sra_dir)

    # send all the files in the directory to SRA
    cmd = "ascp -i %s -v -T -r %s %s " % (argv.private_key, argv.input_dir, dst_dir)
    if os.system(cmd) != 0:
        print("Failed to upload %s to NCBI" % argv.input_dir)
        return(False)

    # write an empty ready file
    ready_file = os.path.join(argv.input_dir, "submit.ready")
    with open(ready_file, "w") as f:
        f.write("    ")

    # send the submit.ready file to SRA
    input_dir_name = os.path.basename(os.path.normpath(argv.input_dir))
    dst_input_dir = "%s%s" % (dst_dir, input_dir_name)
    cmd = "ascp -i %s -v -T %s %s " % (argv.private_key, ready_file, dst_input_dir)
    if os.system(cmd) != 0:
        print("Failed to upload %s to NCBI" % ready_file)
        return(False)

    return(True)

if __name__ == '__main__':
    argv = argparse.ArgumentParser(description = "Send all files to NCBI SRA using ascp command.",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter, epilog = """
        This program is free software: you can redistribute it and/or modify it under the terms of
        the GNU General Public License as published by the Free Software Foundation, either version
        3 of the License, or (at your option) any later version. This program is distributed in the
        hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
        more details. You should have received a copy of the GNU General Public License along with
        this program. If not, see <http://www.gnu.org/licenses/>.""")
    argv.add_argument('-f', '--input-dir', dest = 'input_dir', required = True,
        help = 'Diretory where all of your FASTQ files are located.')
    argv.add_argument('-d', '--ncbi-sra-dir', dest = 'ncbi_sra_dir', required = True,
        help = 'Specify the path to the destination BioProject SRA submission folder')
    argv.add_argument('-i', '--ncbi-private-key', type = str, dest = 'private_key', required = True,
        help = 'Specify the path to your private key file.')
    argv.add_argument('-u', '--ncbi-username', dest = 'ncbi_user', required = True,
        help = 'Username for uploading files to SRA.')

    sys.exit(not main(argv.parse_args()))
