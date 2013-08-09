#!/usr/bin/perl -w
require strict;

my @logs=@ARGV;

my $num_sims=@logs;
my $table1='';
my $table2='';
my $table3='';
my $table4='';

$table1 .= "id & command \\\\\n";
$table1 .= "\\hline\n";

$table2 .= sprintf("id &  time  & \\multicolumn{3}{c|}{Newton} & \\multicolumn{3}{c|}{Krylov} & \\multicolumn{3}{c|}{\$\\alpha\$} \\\\\n");
$table2 .= sprintf("   &  (\$s\$) &  total &   ave  &  max   & total  &  ave   &  max  & \$=1\$ & ave & min \\\\\n");
$table2 .= "\\hline\n";

$table3 .= sprintf("id & num  & num  & zoom & zoom & zoom &   max    &  num \\\\\n");
$table3 .= sprintf("   & grad & back &  lo  & full & bad  &  resid   & indef \\\\\n");
$table3 .= "\\hline\n";

my $agg_time=1;
my $agg_num_newton=1;
my $agg_ave_newton=1;
my $agg_max_newton=1;
my $agg_tot_cg_iter=1;
my $agg_ave_cg_iter=1;
my $agg_max_cg=1;
my $agg_num_alpha_one=1;
my $agg_ave_alpha=1;
my $agg_min_alpha=1;
my $agg_num_grad=0;
my $agg_num_back=0;
my $agg_zoom_lo=0;
my $agg_zoom_full=0;
my $agg_zoom_bad=0;
my $agg_max_grad=0;
my $agg_num_indef=0;

my $i=1;
for my $log (@logs)
{
    open F, "$log/common/log.txt";

    my $command='';

    my $num_steps=0;
    my $time=0;

    my $num_cg=0;
    my $tot_cg_iter=0;
    my $max_cg=0;

    my $num_alpha=0;
    my $tot_alpha=0;
    my $min_alpha=1;
    my $num_alpha_one=0;

    my $max_grad=0;
    my $last_grad=0;

    my $zoom_lo=0;
    my $zoom_full=0;
    my $zoom_bad=0;
    my $exit_zoom=0;

    my $num_indef=0;

    my $newton_this_dt=0;
    my $max_newton=0;

    my $num_grad=0;
    my $num_back=0;

    while(<F>)
    {
        if(/<print>command = (.*) * <\/print>/){$command=$1;next;}

        if(/<stat name=".*iteration.*" value="(.*)"\/>/)
        {
            if($1>$max_cg){$max_cg=$1;}
            $tot_cg_iter+=$1;
            $num_cg++;
            next;
        }

        if(/<print>RAW +(\S+) +(\S+) +(\S+) +(\S+) +(\S+) +(\S+)<\/print>/)
        {
            my $alpha=$6;
            $num_alpha++;
            $num_alpha_one+=($alpha == 1);
            $tot_alpha+=$alpha;
            if($alpha<$min_alpha){$min_alpha=$alpha;}
            if($exit_zoom)
            {
                $exit_zoom=0;
                if($alpha == 1){$zoom_full++;}
                else{$zoom_bad++;}
            }
            $newton_this_dt++;
            next;
        }

        if(/finish_before_indefiniteness/){$num_indef++;next;}

        if(/END Advance_One_Time_Step_Position +(\S+) s/)
        {
            if($last_grad>$max_grad){$max_grad=$last_grad;}
            $num_steps++;
            if($newton_this_dt>$max_newton){$max_newton=$newton_this_dt;}
            $newton_this_dt=0;
            next;
        }

        if(/GRAD STATS +(\S+) +(\S+) +(\S+) *<\/print>/)
        {
            $last_grad=$3;
            next;
        }

        if(/take decrease on zoom with (\S+) *<\/print>/){$zoom_lo++;next;}
        
        if(/exit zoom on/){$exit_zoom=1;next;}

        if(/END Simulation +(\S+) s/){$time=$1;next;}

        if(/LOOK BACKWARDS/){$num_back++;next;}

        if(/GRADIENT/){$num_grad++;next;}
    }
    close F;

    my $num_newton=$num_alpha;
    $table1 .= sprintf("%2d & \\verb|%s| \\\\\n", $i, $command);

    $table2 .= sprintf("%2d & %6.1f & %6d & %6.1f & %6d & %6d & %6.1f & %6d & %10d & %#8.3g & %#8.3g \\\\\n", $i, $time,
                          $num_newton, $num_newton/$num_steps, $max_newton,
                          $tot_cg_iter, $tot_cg_iter/$num_cg, $max_cg,
                          $num_alpha_one, $tot_alpha/$num_alpha, $min_alpha);

    $table3 .= sprintf("%2d & %4d & %4d & %4d & %4d & %4d & %#8.3g & %5d \\\\\n", $i,
                          $num_grad, $num_back, $zoom_lo, $zoom_full, $zoom_bad, $max_grad, $num_indef);

    $agg_time*=$time;
    $agg_num_newton*=$num_newton;
    $agg_ave_newton*=$num_newton/$num_steps;
    $agg_max_newton*=$max_newton;
    $agg_tot_cg_iter*=$tot_cg_iter;
    $agg_ave_cg_iter*=$tot_cg_iter/$num_cg;
    $agg_max_cg*=$max_cg;
    $agg_num_alpha_one*=$num_alpha_one;
    $agg_ave_alpha*=$tot_alpha/$num_alpha;
    $agg_min_alpha*=$min_alpha;
    $agg_num_grad+=$num_grad;
    $agg_num_back+=$num_back;
    $agg_zoom_lo+=$zoom_lo;
    $agg_zoom_full+=$zoom_full;
    $agg_zoom_bad+=$zoom_bad;
    if($agg_max_grad<$max_grad){$agg_max_grad=$max_grad;}
    $agg_num_indef+=$num_indef;

    $i++;
}

$table2 .= "\\hline\n";
$table2 .= sprintf("   & %6.1f & %6.1f & %6d & %6d & %6d & %6.1f & %6d & %10d & %#8.3g & %#8.3g \\\\\n", $agg_time**(1/$num_sims),
                   $agg_num_newton**(1/$num_sims), $agg_ave_newton**(1/$num_sims), $agg_max_newton**(1/$num_sims),
                   $agg_tot_cg_iter**(1/$num_sims), $agg_ave_cg_iter**(1/$num_sims), $agg_max_cg**(1/$num_sims),
                   $agg_num_alpha_one**(1/$num_sims), $agg_ave_alpha**(1/$num_sims), $agg_min_alpha**(1/$num_sims));

$table3 .= "\\hline\n";
$table3 .= sprintf("   & %4d & %4d & %4d & %4d & %4d & %#8.3g & %5d \\\\\n",
                   $agg_num_grad, $agg_num_back, $agg_zoom_lo, $agg_zoom_full, $agg_zoom_bad, $agg_max_grad, $agg_num_indef);

open O, ">results.tex";

print O <<EOF;
\\documentclass[10pt]{article}
\\usepackage{fullpage}
\\begin{document}

\\begin{tabular}{|rl|}
\\hline
$table1\\hline
\\end{tabular}
\\clearpage

\\begin{tabular}{|rr|rrr|rrr|rrr|}
\\hline
$table2\\hline
\\end{tabular}
\\clearpage

\\begin{tabular}{|rrrrrrrr|}
\\hline
$table3\\hline
\\end{tabular}
\\clearpage

\\end{document}
EOF

close O;

system("pdflatex -halt-on-error results.tex");

