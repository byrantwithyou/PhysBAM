#!/usr/bin/perl -w
require strict;
my $PHYSBAM=$ENV{'PHYSBAM'};
#for my $s ((0,1,2,3,4,5,6,8,9,10,11,12,13,14,16,17,18,19,20,21,24,28,33))
for my $s ((112))
{
for my $j ((2..2))
{
    my $mu = 2**($j-2);
    for my $i (0..9)
    {
        my $dt = 2**($i-13);
        for my $k ((3,6))
        {
            my $rho = 2**($k-4);
            for my $r (32,64)
            {
                my $pr=80-$s;
                print `$PHYSBAM/Tools/batch/slave -p $pr -- nice ./fluids_color -resolution $r -last_frame 200 -dt $dt $s -rho0 $rho -mu0 $mu -o shed-stab-$s-$j-$i-$k-$r-fix`;
                sleep(1);
            }
        }
    }
}
}


