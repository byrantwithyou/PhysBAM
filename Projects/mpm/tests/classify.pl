#!/usr/bin/perl
open(my $F,'<',$ARGV[0]);

my $finished=0;
my $fp_inv_op=0;
my $fp_div_zero=0;
my $assert_line=0;
my $fail_strain=0;
my $precompute=0;
my $sigpfe=0;
my $max_vel=0;
my $min_dt=1;
my $vel_limit=10;
my $dt_limit=1e-6;
my $fail=0;
my $vel_too_big=0;

while(<$F>)
{
    if(/END frame 100/){$finished=1;}
    if(/Floating point exception: reason 7 = "FP invalid operation"/){$fp_inv_op=1;}
    if(/Floating point exception: reason 7 = "FP divide by zero"/){$fp_div_zero=1;}
    if(/Line_Search_Wolfe_Conditions.*Assertion failed, condition = a0.dfa<0/){$assert_line=1;}
    if(/P_From_Strain/){$fail_strain=1;}
    if(/Compute_Precompute_Data/){$precompute=1;}
    if(/SIGNAL SIGFPE/){$sigpfe=1;}
    if(/max velocity: (\S+)/){if($1>$max_vel){$max_vel=$1;}}
    if(/substep dt: (\S+)/){if($1<$min_dt){$min_dt=$1;}}
    if(/VELOCITY TOO BIG/){$vel_too_big=1;}
}

print "$ARGV[0]:";
$finished && print " FINISHED";
$fp_inv_op && print " INVALID-OP";
$fp_div_zero && print " DIV-ZERO";
$assert_line && print " ASSERT-LS";
$fail_strain && print " FAIL-STRAIN";
$precompute && print " FAIL-PRECOMPUTE";
$sigpfe && print " SIGPFE";
($max_vel>$vel_limit || $vel_too_big) && print " BIG-V";
$min_dt<$dt_limit && print " SMALL-DT";
printf(" (%.1f)",$max_vel);
print "\n";
