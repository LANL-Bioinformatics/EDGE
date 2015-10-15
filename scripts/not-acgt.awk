#!/bin/awk -f
# Usage:  not-acgt.awk
#   Read a fasta input file and find regions consisting of MIN_RUN
#   or more consecutive non-acgt characters in the first string.
#   If there is more than one string, all strings after the first
#   are ignore.  Output is one line
#   per region, with start position and end position on each line.
#   Positions are inclusive, counting from 1 so that the first 10
#   positions of the file are indicated as "1 10".  The value of
#   MIN_RUN can be set below.


BEGIN {
  MIN_RUN = 5;
  ct = pos = start = 0;
}


/^>/ {
  line_ct ++;
  if (line_ct == 1)
    next;
  else
    exit;
}


{
  n = length ($1);
  for (i = 1; i <= n; i ++)
    {
      if (match (substr ($1, i, 1), /[acgtACGT]/))
        Pr();
      else
        {
          if (ct == 0)
            start = pos + 1;
          ct ++;
        }
      pos ++;
    }
}


END {
  Pr();
}


function  Pr  ()
{
  if (ct >= MIN_RUN)
    printf "%8d %8d\n", start, pos;
  ct = 0;
}
