#!/usr/bin/env python3
import sys
import requests
import argparse as ap

#example curl cmd
#curl -X POST -F 'db=wuhCor1' -F 'sarsCoV2File=@EPI_ISL_462810_gappy.fa' https://hgwdev-demo-angie.gi.ucsc.edu/cgi-bin/hgPhyloPlace 
def parse_params():
    p = ap.ArgumentParser(prog='run_remote_UShER.py',
            description="""hgPhyloPlace: Ultrafast Sample placement on Existing SARS-CoV-2 Tree """)
    p.add_argument('-f', '--fasta',
                        metavar='[FILE]', type=ap.FileType(), required=True,
                        help="sequence file in FASTA format")
    p.add_argument('-o', '--out',
                        metavar='[HTML]', type=str, default='out.html',
                        help="output HTML file [out.html]")
    args_parsed = p.parse_args()
   # if not args_parsed.out:
   #     args_parsed.out = "out.html"
    return args_parsed

def main():
    argvs = parse_params()

    # api-endpoint
    URL = "https://hgwdev-demo-angie.gi.ucsc.edu/cgi-bin/hgPhyloPlace"

    sendData = {'db':'wuhCor1'}
    sendFiles = {'sarsCoV2File':argvs.fasta}

    try:
        # sending get request and saving the response as out file
        r = requests.post(url = URL, data = sendData, files= sendFiles)
        with open(argvs.out, 'w') as outfile:
            outfile.write(r.text.replace("../js/", "js/").replace("'../","'https://hgwdev-demo-angie.gi.ucsc.edu/").replace('"../','"https://hgwdev-demo-angie.gi.ucsc.edu/'))
            #outfile.write(r.text.replace("../style/", "style/").replace("../js/", "js/").replace("../index.html", "index.html"))
    finally:
        pass


if __name__ == "__main__":
    main()


