#!/usr/bin/perl -w

require strict;

my $prep = 0;
my $solving = 0;
while(<STDIN>)
{
    if(/args *(.*) ms/){$prep+=$1}
    if(/parse input *(.*) ms/){$prep+=$1}
    if(/setup viewing *(.*) ms/){$prep+=$1}
    if(/masters *(.*) ms/){$prep+=$1}
    if(/merge blocks *(.*) ms/){$prep+=$1}

    if(/compute matrix *(.*) ms/){$solving+=$1}
    if(/elim irreg *(.*) ms/){$solving+=$1}
    if(/elim non sep *(.*) ms/){$solving+=$1}
    if(/elim 3 *(.*) ms/){$solving+=$1}
    if(/back solve *(.*) ms/){$solving+=$1}
    if(/exec jobs *(.*) ms/)
    {
        $solving+=$1;
        printf("%.3f %.3f %.3f\n",$prep/1000,$solving/1000,($prep+$solving)/1000);
    }
}
