#!/bin/awk -f
# Usage:  glim-diff.awk  <a-pred> <b-pred>
#   Read gene predictions in <a-pred> and <b-pred>
#   and output them side by side.  Both must be
#   in sorted order by stop codon and the format for
#   each must be:
#     <id>  <start>  <stop>  [additional columns irrelevant]
#   Also print summary info at end.


BEGIN   {
         if  (ARGC < 3)
             Usage_Exit();

         afp = ARGV [1];
         delete ARGV [1];
         bfp = ARGV [2];
         delete ARGV [2];

         Read_A();
         Read_B();

         while  (! (adone || bdone))
           {
            if  (1 * aend < 1 * bend)
                {
                 printf "%-8s %7d %7d  <\n", aid, astart, aend;
                 aonly ++;
                 Read_A();
                }
            else if  (1 * bend < 1 * aend)
                {
                 printf "%24s  >  %-8s %7d %7d\n", "", bid, bstart, bend;
                 bonly ++;
                 Read_B();
                }
              else
                {
                 if  (1 * astart < 1 * aend)
                     diff = bstart - astart;
                   else
                     diff = astart - bstart;
                 if  (diff == 0)
                     {
                      ch = "=";
                      exact_ct ++;
                     }
                   else
                     ch = "|";
                 printf "%-8s %7d %7d  %s  %-8s %7d %7d\n",
                      aid, astart, aend, ch, bid, bstart, bend;
                 match_ct ++;
                 diff_sum += diff;
                 Read_A();
                 Read_B();
                }
           }

         while  (! adone)
           {
            printf "%-8s %7d %7d  <\n", aid, astart, aend;
            aonly ++;
            Read_A();
           }
         while  (! bdone)
           {
            printf "%24s  >  %-8s %7d %7d\n", "", bid, bstart, bend;
            bonly ++;
            Read_B();
           }

         print "";
         printf " A only: %6d  %5.1f%%\n", aonly, Percent(aonly, acount);
         printf " B only: %6d  %5.1f%%\n", bonly, Percent(bonly, bcount);
         printf "Matches: %6d  %5.1f%%  %5.1f%%\n", match_ct,
              Percent(match_ct, acount), Percent(match_ct, bcount);
         printf "  Exact: %6d  %5.1f%%  %5.1f%%\n", exact_ct,
              Percent(exact_ct, match_ct), Percent(exact_ct, acount);
         printf "AvgDiff: %8.1f\n", diff_sum / match_ct;
         printf "A count: %6d\n", acount;
         printf "B count: %6d\n", bcount;
        }



function  Percent  (x, y)
  {
   if  (y == 0)
       return  0.0;
     else
       return  (100.0 * x) / y;
  }



function  Read_A  ()
  {
   if  ((getline < afp) > 0)
       {
        aid = $1;
        astart = $2;
        aend = $3;
        acount ++;
       }
     else
       adone = 1;
  }



function  Read_B  ()
  {
   if  ((getline < bfp) > 0)
       {
        bid = $1;
        bstart = $2;
        bend = $3;
        bcount ++;
       }
     else
       bdone = 1;
  }



function  Usage_Exit  ()
  {
   print "# Usage:  glim-diff.awk  <a-pred> <b-pred>";
   print "#   Read gene predictions in <a-pred> and <b-pred>";
   print "#   and output them side by side.  Both must be";
   print "#   in sorted order by stop codon and the format for";
   print "#   each must be:";
   print "#     <id>  <start>  <stop>  [additional columns irrelevant]";
   print "#   Also print summary info at end.";

   exit;
  }
