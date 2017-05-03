#!/usr/bin/perl -w
require strict;
my $top_dir=`git rev-parse --show-toplevel`;
chomp $top_dir;
my $message=`git log --format=%B -n 1 HEAD`;
chomp $message;
chomp $message;
chomp $message;
`git reset --soft HEAD^`;
my $files = `git diff --name-only HEAD`;

my @working = map {[$_,$_]} split /\s/, $files;

my @done = ();

my %hash_done = ();

while(@working)
{
    my %hash = ();
    my @L = ();
    for(@working)
    {
        my ($full,$short)=@$_;
        if(!($short =~ /\//))
        {
            push @done, [$full, $short];
            $hash_done{$short}=1;
            $hash{$short}=1;
        }
        else
        {
            my $s = "$'";
            if(defined $hash_done{$s} || defined $hash{$s})
            {
                $hash{$s}=1;
                push @done, [$full, $short];
                $hash_done{$short}=1;
                $hash{$short}=1;
            }
            else
            {
                $hash{$s}=0;
                push @L, [$full, $short, $s];
            }
        }
    }
    @working=();
    for(@L)
    {
        my ($full,$short,$new_short)=@$_;
        if($hash{$new_short}==1)
        {
            push @done, [$full, $short];
            $hash_done{$short}=1;
            $hash{$short}=1;
        }
        else
        {
            push @working, [$full, $new_short];
        }
    }
}
chdir($top_dir);
for(@done)
{
    my ($full,$short)=@$_;
    system("git","commit",$full,"-m","$message $short");
}
