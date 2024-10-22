#!/usr/bin/perl -w

require strict;

my %lookup;
my %local_lookup;

open F, shift @ARGV;

while(<F>){
    chomp;
    if(/.*\/(.*)/)
    {
        (defined $lookup{$1}) and die "Duplicate file '$1' encountered\n.";
        $lookup{$1} = $_;
    }
    elsif(!defined $local_lookup{$_})
    {
        $local_lookup{$_} = $_;
    }
}
close F;

my @order = qw( Core Tools Particles_Tools Grid_Tools Grid_PDE Geometry Rigids Deformables Solids Hybrid_Methods Incompressible Compressible Fluids Dynamics OpenGL Ray_Tracing );
my %category = ();

my @local_includes = ();
my @system = ();
my %system = ("algorithm" => 1, "streambuf" => 1, "iomanip" => 1, "list" => 1, "ostream" => 1,
              "bitset" => 1, "ios" => 1, "locale" => 1, "queue" => 1, "string" => 1, 
              "complex" => 1, "typeinfo" => 1, "iosfwd" => 1, "map" => 1, "set" => 1, 
              "deque" => 1, "iostream" => 1, "memory" => 1, "sstream" => 1, "utility" => 1, 
              "exception" => 1, "valarray" => 1, "istream" => 1, "new" => 1, "stack" => 1, 
              "fstream" => 1, "iterator" => 1, "numeric" => 1, "stdexcept" => 1, "vector" => 1, 
              "functional" => 1, "limits" => 1, "cassert" => 1, "ciso646" => 1, "csetjmp" => 1, 
              "cstdio" => 1, "ctime" => 1, "cctype" => 1, "climits" => 1, "csignal" => 1, 
              "cstdlib" => 1, "cwchar" => 1, "cerrno" => 1, "clocale" => 1, "cstdarg" => 1, 
              "random" => 1,
              "cstring" => 1, "cwctype" => 1, "cfloat" => 1, "cmath" => 1, "cstddef" => 1);

my $cf = "";

my $new;
my $filename;

my $last = "";

sub dump_headers_set
{
    my $L=$_[0];
    for(sort {$x=lc $$a[0];$y=lc $$b[0];$x=~s/_/0/g;$y=~s/_/0/g;$x cmp $y;} @$L)
    {
#        print "try $$_[0] in $filename (vs $last)\n";
        if($$_[0] eq $last){print  "DUPLICATE HEADER: $last in $filename\n";next;}
        $last = $$_[0];
        my $end_comment=$$_[2]?" $$_[2]":"";
        $new.="$$_[1]#include $$_[0]$end_comment\n";
    }
}

sub dump_headers
{
    if(! scalar keys %category && !@system && !@other && !@local_includes){return;}
    $filename eq $cf or die "Trying to insert includes from $cf into $filename.\n";
    for my $o (@order)
    {
        if(defined $category{$o})
        {
            &dump_headers_set($category{$o});
            delete $category{$o};
        }
    }
    for my $o (sort keys %category)
    {
        &dump_headers_set($category{$o});
        delete $category{$o};
    }
    &dump_headers_set(\@system);
    &dump_headers_set(\@local_includes);
    &dump_headers_set(\@other);
    @system = ();
    @local_includes = ();
    @other = ();
    $cf = "";
    $last = "";
}
sub try_find_file
{
    my ($file,$pre,$post)=@_;
    if(defined $lookup{$file})
    {
        $lookup{$file}=~/^(.*?)\//;
        if(!defined $category{$1}){$category{$1}=[];}
        my $L=$category{$1};
        push @$L, ["<$lookup{$file}>",$pre,$post];
        return 1;
    }
    if(defined $local_lookup{$file})
    {
        push @local_includes, ["\"$local_lookup{$file}\"",$pre,$post];
        return 1;
    }
    return 0;
};
for $filenamex (@ARGV)
{
    $filename=$filenamex;
    my $orig='';
    $new='';

    open F, $filename;
    while(<F>)
    {
        $orig.=$_;
        if(!m@^\s*(//|)\s*\#\s*include\s*(([<\"])(.*/)?(.*)([>\"]))\s*(//.*\S|)@){&dump_headers();$new.=$_;next;}
        if($cf eq ""){$cf = $filename;}
        $filename eq $cf or die "Trying to mix includes from $cf into $filename.\n";
        my $ldelim = $3;
        my $rdelim = $6;
        my $full = $2;
        my $path = $4;
        my $file = $5;
        my $pre = $1;
        my $post = $7;
        if(&try_find_file($file,$pre,$post)){next;}
        if(/([a-zA-Z_]+)\.h/)
        {
            if(&try_find_file(uc($1).".h",$pre,$post)){next;}
        }
        if(!$path){push @other, [$full,$pre,$post];next;}
        print  "include not recognized: '$file' in '$filename'\n";
        push @other, [$full,$pre,$post];
    }
    &dump_headers();
    close F;

    if($new ne $orig)
    {
        open F, ">$filename";
        print F $new;
        close F;
    }
    $filename='qq';
}
