#!/bin/awk -f
# Usage:  get-motif-counts.awk
#   Read output of  elph  program and extract the motif
#   count information.  Output it in a format suitable for
#   reading by  glimmer3 .


/^Motif counts:/  {
         state = 1;
         next;
        }

        {
         if  (state && match ($0, /^[acgt]:/))
             {
              if  (width == 0)
                  {
                   width = NF - 1;
                   print width;
                  }
              printf "%s", substr ($1, 1, 1);
              for  (i = 2;  i <= NF;  i ++)
                printf " %7d", $(i);
              printf "\n";
             }
        }


