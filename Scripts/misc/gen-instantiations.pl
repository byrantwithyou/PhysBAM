#!/usr/bin/perl -w
my %func=();
my %class=();

while(<>)
{
    if(!/undefined reference to \`(.*)\'/){next}
    $_=$1;
    s/PhysBAM:://g;
    s/, /,/g;
    $func{$_}=1;
    s/::.*//;
    $class{$_}=1;
}

for(sort keys %class)
{
    print "template class $_;\n";
}

print "\n------------------------------------------------------------------\n";

for(sort keys %func)
{
    print "template $_;\n";
}
