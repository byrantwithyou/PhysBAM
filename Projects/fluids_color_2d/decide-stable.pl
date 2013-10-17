#!/usr/bin/perl -w
require strict;

my %hs=();

my %mu=();

for my $file (@ARGV)
{
    
    open F, "$file/common/log.txt";

    my $mx=0;
    my $mu=0;
    my $rho=0;
    my $dt=0;
    while(<F>)
    {
        if(/max u ([-+0-9.eE]+)/)
        {
            if($1>$mx){$mx=$1;}
        }
        if(/fluids_color/)
        {
            if(/-dt (\S+)/){$dt=$1;}
            if(/-rho0 (\S+)/){$rho=$1;}
            if(/-mu0 (\S+)/){$mu=$1;}
        }
    }

    close F;
    my $class='';
    if(! -e "$file/1000" || $mx>=1000){$class='U';}
    elsif($mx<10){$class='S';}
    else{$class='U';}
    $hs{"$class-$mu.txt"}.="$mu $rho $dt\n";
    $mu{$mu}=1;
}

for $k (keys %hs)
{
    open F, ">$k";

    print F $hs{$k};

    close F;
}

for $m (keys %mu)
{
    my $cmd = "set terminal postscript eps color size 2,2 ;";
    $cmd .= "set output \"stability-plot-$m.eps\" ;";
    $cmd .= "set size square ;";
    $cmd .= "plot \"S-$m.txt\" u (log10(\$2)):(log10(\$3)) ls 31 notitle , \"U-$m.txt\" u (log10(\$2)):(log10(\$3)) ls 197 notitle";
    `gnuplot -e '$cmd' ; sed -i 's/LC3 {1 0 1} def/LC3 {0 1 0} def/' stability-plot-$m.eps`;
}
