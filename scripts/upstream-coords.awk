#!/bin/awk -f
# Usage:  upstream-coords.awk  <len>  <separation>
#   Read gene prediction coordinates from standard input
#   and output the coordinates of the region of length
#    <len>  that is  <sep>  bases before the 5' start
#   of the gene.  Input format is:
#     <tag>  <start>  <stop>
#   Output format is the same.
#   If the length of the gene is longer than  MAX_GENE_LEN ,
#   then the gene is assumed to wrap around a circular genome
#   Note that output coordinates can be negative or longer
#   than the genome length (which is unknown).


BEGIN   {
         if  (ARGC < 3)
             Usage_Exit();

         if  (MAX_GENE_LEN == 0)
             MAX_GENE_LEN = 100000;

         len = ARGV [1];
         delete ARGV [1];

         sep = ARGV [2];
         delete ARGV [2];
        }


        {
         if  (1 * $2 < $3)
             {
              gene_len = 1 + $3 - $2;
              dir = 1;
             }
           else
             {
              gene_len = 1 + $2 - $3;
              dir = -1;
             }
         if  (gene_len > MAX_GENE_LEN)
             dir *= -1;

         printf "%s %8d %8d\n", $1, $2 - dir * (sep + len),
              $2 - dir * (sep + 1);
        }



function  Usage_Exit  ()
  {
   print "# Usage:  upstream-coords.awk  <len>  <separation>";
   print "#   Read gene prediction coordinates from standard input";
   print "#   and output the coordinates of the region of length";
   print "#    <len>  that is  <sep>  bases before the 5' start";
   print "#   of the gene.  Input format is:";
   print "#     <tag>  <start>  <stop>";
   print "#   Output format is the same.";
   print "#   If the length of the gene is longer than  MAX_GENE_LEN ,";
   print "#   then the gene is assumed to wrap around a circular genome";
   print "#   Note that output coordinates can be negative or longer";
   print "#   than the genome length (which is unknown).";

   exit;
  }
