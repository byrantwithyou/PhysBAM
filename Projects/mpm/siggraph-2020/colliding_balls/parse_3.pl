#!/usr/bin/perl -w

require strict;

my $inf=1e300;

my $valid=0;
my $cfl_v=0;
my $cfl_f=0;
my $cfl_c=0;
my $max_v=0;
my $max_c=0;
my $max_f=0;
my $min_dt=0;
my $substeps=0;

my $best_cfl_v=$inf;
my $best_cfl_f=$inf;
my $best_cfl_c=$inf;
my $best_max_v=0;
my $best_max_c=0;
my $best_max_f=0;
my $best_min_dt=$inf;

sub min {return $_[0]<$_[1]?$_[0]:$_[1];}
sub max {return $_[0]<$_[1]?$_[1]:$_[0];}

my $time=0;

print "time best_cfl_v best_cfl_f best_cfl_c best_max_v best_max_c best_max_f best_min_dt\n";

sub commit_substep
{
    if($valid)
    {
        $best_cfl_v=&min($best_cfl_v,$cfl_v);
        $best_cfl_f=&min($best_cfl_f,$cfl_f);
        $best_cfl_c=&min($best_cfl_c,$cfl_c);
        $best_max_v=&max($best_max_v,$max_v);
        $best_max_c=&max($best_max_c,$max_c);
        $best_max_f=&max($best_max_f,$max_f);
        $best_min_dt=&min($best_min_dt,$min_dt);
    }
    $cfl_v=$inf;
    $cfl_f=$inf;
    $cfl_c=$inf;
    $max_v=0;
    $max_c=0;
    $max_f=0;
    $min_dt=$inf;
    $valid=1;        
}

while(<STDIN>)
{
    if(/<print>substep dt: (.+)<\/print>/){$cfl_v=$1;next;}
    if(/<print>F CFL (.+)<\/print>/){$cfl_f=$1;next;}
    if(/<print>SOUND CFL .+ (.+) \(.+\)<\/print>/){$cfl_c=$1;next;}
    if(/<print>PLOT (.+) (.+) (.+) (.+)<\/print>/)
    {
        $min_dt=$1;
        $max_v=$2;
        $max_c=$3;
        $max_f=$4;
        $substeps++;
        next;
    }
    if(/<scope id="FRAME"/ || /^</)
    {
        if($substeps==1)
        {
            $valid=1;
            &commit_substep();
        }
        if($valid)
        {
            $time += 1./24;
            print "$time $best_cfl_v $best_cfl_f $best_cfl_c $best_max_v $best_max_c $best_max_f $best_min_dt\n";
        }
        $valid=0;
        $substeps=0;
        $best_cfl_v=$inf;
        $best_cfl_f=$inf;
        $best_cfl_c=$inf;
        $best_max_v=0;
        $best_max_c=0;
        $best_max_f=0;
        $best_min_dt=$inf;
        next;
    }
    if(/<scope id="SUBSTEP"/)
    {
        &commit_substep();
    }
}
