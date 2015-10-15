#!/bin/awk -f
# Usage:  match-list-col.awk  <list-file>  <col>
#   Print lines from stdin whose entry in column <col> is one of the entries
#   occurring in <list-file>.

BEGIN   {
         if  (ARGC < 3)
             Usage();

         fp = ARGV [1];
         delete ARGV [1];
         while  ((getline < fp) > 0)
             {
              list [$1] = 1;
             }

         col = ARGV [2];
         delete ARGV [2];
         match (col, /[0-9]*/);
         if  (RSTART != 1 || RLENGTH != length (col))
             {
              printf "ERROR:  Bad column value = %s\n", col;
              Usage();
             }
        }

        {
         if  (list [$(col)] == 1)
             print;
        }

function  Usage  ()
  {
   print "# Usage:  match-list-col.awk  <list-file> <col>";
   print "#   Print lines from stdin whose entry in column <col> is one of the";
   print "#   entries occurring in <list-file>.";
   exit;
  }
