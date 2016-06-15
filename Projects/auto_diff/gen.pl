#!/usr/bin/perl -w

require strict;

my $s=21;
my $u=10;
my $v=15;
my $m=0;
my $n=0;
my $p=0;
my $q=0;

for my $i (0..$s){print "id(s$i,\"s$i\");\n";}
for my $i (0..$u){print "id(u$i,\"u$i\");\n";}
for my $i (0..$v){print "id(v$i,\"v$i\");\n";}
for my $i (0..$m){print "id(m$i,\"m$i\");\n";}
for my $i (0..$n){print "id(n$i,\"n$i\");\n";}
for my $i (0..$p){print "id(p$i,\"p$i\");\n";}
for my $i (0..$q){print "id(q$i,\"q$i\");\n";}

die "";

for my $i (0..$s){print "OP_s(s$i);\n";}

for my $i (0..$s){for my $j (0..$s){print "OP_ss(s$i,s$j);\n";}}
# for my $i (0..$s){for my $j (0..$s){for my $k (0..$s){print "OP_sss(s$i,s$j,s$k);\n";}}}

for my $i (0..$u){print "OP_u(u$i);\n";}
for my $i (0..$u){for my $j (0..$u){print "OP_uu(u$i,u$j);\n";}}
for my $i (0..$u){for my $j (0..$s){print "OP_us(u$i,s$j);\n";}}

for my $i (0..$v){print "OP_v(v$i);\n";}
for my $i (0..$v){for my $j (0..$v){print "OP_vv(v$i,v$j);\n";}}
for my $i (0..$v){for my $j (0..$s){print "OP_vs(v$i,s$j);\n";}}

for my $i (0..$m){print "OP_m(m$i);\n";}
for my $i (0..$n){print "OP_n(n$i);\n";}
for my $i (0..$p){print "OP_p(p$i);\n";}
for my $i (0..$q){print "OP_q(q$i);\n";}

for my $i (0..$m){for my $j (0..$u){print "OP_mu(m$i,u$j);\n";}}

