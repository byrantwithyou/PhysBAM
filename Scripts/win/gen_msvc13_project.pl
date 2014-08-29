#!/usr/bin/perl -w
require strict;

# file.pl dir name lib -- list-file list-file ... -- token token ...

my $argi = 0;

my $dir = $ARGV[$argi++];
my $name = $ARGV[$argi++];
my $lib = $ARGV[$argi++];

use Digest::MD5 qw(md5_hex);
my $gid = md5_hex("$dir\\msvc13\\$name");
$gid=uc($gid);
$gid=~s/(.{8})(.{4})(.{4})(.{4})(.{12})/$1-$2-$3-$4-$5/;

$ARGV[$argi++] eq "--" or die "expected -- in args.";

my @rep_file = ();
for(;$ARGV[$argi] ne "--";$argi++)
{
    push @rep_file, $ARGV[$argi];
}

$ARGV[$argi++] eq "--" or die "expected -- in args.";

my @replacements = ();
for(;$argi<@ARGV;$argi++)
{
    push @replacements, $ARGV[$argi];
}

my @special = ($name, $gid);

my $pb = $ENV{"PHYSBAM"};

my @sources = ();
my @headers = ();

my $ind = 4;
for(; $ind < @ARGV && $ARGV[$ind] ne '--'; $ind++)
{
    push @sources, $ARGV[$ind];
}
for($ind++; $ind < @ARGV; $ind++)
{
    push @headers, $ARGV[$ind];
}

my $state = 'x';
my $accum = '';

sub flush_accum
{
    if($state=~/f([0-9]+)/)
    {
        my $file = $rep_file[$1];
        open RF, "$file";
        while(<RF>)
        {
            chomp;
            my $a = $accum;
            $a =~ s/#1/$_/g;
            print $a;
        }
        close RF;
    }
    elsif($state=~/s([0-9]+)/)
    {
        my $val = $special[$1];
        $accum =~ s/#1/$val/g;
        print $accum;
    }
    elsif($state=~/r([0-9]+)/)
    {
        my $val = $replacements[$1];
        $accum =~ s/#1/$val/g;
        print $accum;
    }
    elsif($state=~/l([0-9]+)/)
    {
        if($1 == $lib)
        {
            print $accum;
        }
    }
    elsif($state=~/p([0-9]+)/)
    {
        my $file = $rep_file[$1];
        open RF, "$file";
        while(<RF>)
        {
            /(\S+) (\S+)/ or die "failed to parse project info: $_";
            my $proj_name = $1;
            my $proj_dir = $2;
            my $proj_gid = md5_hex("$proj_dir\\msvc13\\$proj_name");
            $proj_gid=uc($proj_gid);
            $proj_gid=~s/(.{8})(.{4})(.{4})(.{4})(.{12})/$1-$2-$3-$4-$5/;
            my $a = $accum;
            $a =~ s/#1/$proj_name/g;
            $a =~ s/#2/$proj_dir/g;
            $a =~ s/#3/$proj_gid/g;
            print $a;
        }
        close RF;
    }
    else
    {
        die "envalid state $state\n";
    }
    $accum = '';
    $state = 'x';
}
while(<STDIN>)
{
    if(/^@([a-z][0-9]*)@/)
    {
        my $newstate=$1;
        my $rest="$'";
        if($newstate eq $state)
        {
            $accum.=$rest;
            next;
        }
        elsif($state ne 'x')
        {
            &flush_accum();
        }
        $state=$newstate;
        $accum=$rest;
        next;
    }
    my $line = $_;
    if($state ne 'x')
    {
        &flush_accum();
    }
    print $line;
}
if($state ne 'x')
{
    &flush_accum();
}
