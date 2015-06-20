#!/usr/bin/perl -w
require strict;

my $num=1000;

my @ops=(
    ['s','(-#s)'],
    ['s','(+#s)'],
    ['s','(#s-#s)'],
    ['s','(#s+#s)'],
    ['s','(#s*#s)'],
    ['s','(#s/(sqr(#s)+1))'],
    ['s','sqrt(abs(#s)+1)'],
    ['s','sqr(#s)'],
    ['s','cube(#s)'],
    ['s','max(#s,#s)'],
    ['s','min(#s,#s)'],
    ['s','abs(#s)'],
    ['s','log(abs(#s)+1)'],
    ['s','exp(#s)'],
    ['s','sin(#s)'],
    ['s','cos(#s)'],
    ['s','tan(#s)'],
    ['s','hypot(#s,#s)'],
    ['s','hypot(#s,#s,#s)'],
    ['s','atan2(#s,#s)'],
    ['v','(-#v)'],
    ['v','(+#v)'],
    ['v','(#v-#v)'],
    ['v','(#v+#v)'],
    ['v','(#v*#s)'],
    ['v','(#v/(abs(#s)+1))'],
    ['v','(#s*#v)'],
    ['v','(#v.Cross(#v))'],
    ['s','(#v.Dot(#v))'],
    ['s','(#v.Magnitude_Squared())'],
    ['s','(#v.Magnitude())']);

my @s=('a','b','#r');
my @v=('u','v','TV(#r,#r,#r)');

for(1..$num)
{
    my @op=@{$ops[rand()*@ops]};
    my $var=$op[0];
    $_=$op[1];
    my $vs=int(rand()*2);
    while(/#v/)
    {
        my $v=$v[int(rand()*@v)];
        s/#v/$v/;
    }
    while(/#s/)
    {
        my $v=$s[int(rand()*@s)];
        s/#s/$v/;
    }
    while(/#r/)
    {
        my $v=int(rand()*200)/100-1;
        s/#r/$v/;
    }
    if(/.{80}/){next;}
    print "    TEST($_);\n";
    if($var eq 's'){push @s, $_;}
    if($var eq 'v'){push @v, $_;}
}
