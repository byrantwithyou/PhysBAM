#!/usr/bin/perl -w

require strict;
use File::Path "remove_tree";

# Usage: ./parse_3.pl out_ok

my $SIM=$ARGV[0];
$SIM=~/out_(.*)/ or die "expected out_*\n";
my $TEX="tex_$1";
my $RENDER="render_$1";

my $pt_0=5.8;
my $pt_1=12.55;

remove_tree($TEX);
mkdir($TEX);

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

my $clip_dt=.001;
my $clip_max_v=40;
my $clip_max_f=10;
my $clip_max_c=10;

my $frame = 0;


sub min {return $_[0]<$_[1]?$_[0]:$_[1];}
sub max {return $_[0]<$_[1]?$_[1]:$_[0];}

my $time=0;

my $all_data="time best_cfl_v best_cfl_f best_cfl_c best_max_v best_max_c best_max_f best_min_dt\n";

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

open I, "template_3.tex";
my $template='';
while(<I>){$template.=$_;}
close I;

open L, "$SIM/common/log.txt";
while(<L>)
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
            $time += 1.0/24;
            
            $best_cfl_v=&min($best_cfl_v/$clip_dt,.98);
            $best_cfl_f=&min($best_cfl_f/$clip_dt,.98);
            $best_cfl_c=&min($best_cfl_c/$clip_dt,.98);
            $best_min_dt=&min($best_min_dt/$clip_dt,.98);
            $best_max_v=&min($best_max_v/$clip_max_v,.98);
            $best_max_f=&min($best_max_f/$clip_max_f,.98);
            $best_max_c=&min($best_max_c/$clip_max_c,.98);

            my $n=sprintf("%04d",$frame);
            
            $all_data.="$time $best_cfl_v $best_cfl_f $best_cfl_c $best_max_v $best_max_c $best_max_f $best_min_dt\n";
            open O, ">$TEX/data-$n.txt";
            print O $all_data;
            close O;

            my $dtv=$best_cfl_v*($pt_1-$pt_0)+$pt_0;
            my $dtc=$best_cfl_c*($pt_1-$pt_0)+$pt_0;
            my $dtf=$best_cfl_f*($pt_1-$pt_0)+$pt_0;
            my $dtb=&min(&min($best_cfl_v,$best_cfl_c),$best_cfl_f);
            $_=$template;
            s@XXXX@../$RENDER/colliding_balls_mantra_ipr_$n.png@g;
            s/DATA/data-$n.txt/g;
            s/DTV/$dtv/g;
            s/DTF/$dtf/g;
            s/DTC/$dtc/g;
            my $colv="blue!50!white";
            my $colf="blue!50!white";
            my $colc="blue!50!white";
            if($best_cfl_f==$dtb){$colf="blue";}
            if($best_cfl_v==$dtb){$colv="blue";}
            if($best_cfl_c==$dtb){$colc="blue";}
            s/COLV/$colv/g;
            s/COLF/$colf/g;
            s/COLC/$colc/g;
            open O, ">$TEX/comp-$n.tex";
            print O $_;
            close O;
            system("( cd $TEX ; pdflatex -halt-on-error comp-$n.tex > /dev/null ; convert comp-$n.pdf comp-$n.png ; echo $n ) &");
            $frame++;
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
close L;
wait;
