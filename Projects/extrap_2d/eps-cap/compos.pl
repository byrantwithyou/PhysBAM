#!/usr/bin/perl -w
require strict;

my $mg=10;
my $h=480;
my $w=640;
my $sr=3;
my $s=($h-2*$mg)/(2*$sr);
my $sh=$mg + $sr * $s;
my $arlen=.135;

my $d="([0-9.eE+-]+)";
my $nd="[^0-9.eE+-]+";

sub pt
{
    my $x=$_[0] * $s + $sh;
    my $y=$_[1] * $s + $sh;
    return sprintf("(%.2f,%.2f)", $x, $y);
}

sub spt
{
    $_[0]=~/$d$nd$d/;
    return &pt($1,$2);
}

my $arrows='';
my $arrow='';
my $vertline='';
my @trails=();
my @triangles=();
my @cc=qw(coltri1 coltri2 coltri3);
my %contour=();
my $tristretchstr='';
my $tristretchpts='';
$"=",";

my $tristretch=0;
while(<>)
{
    if(/p $nd$d$nd$d$nd$d$nd$d\]/)
    {
        if($1 + $2 > 0)
        {
            my $xx=$3*$arlen+$1;
            my $yy=$4*$arlen+$2;
            $arrows.="\\psline[linewidth=1px,linecolor=colarrow]{->}@{[&pt($1,$2)]}@{[&pt($xx,$yy)]}\n";
            $arrows.="\\psline[linewidth=1px,linecolor=colarrowref]{->}@{[&pt(-$1,-$2)]}@{[&pt(-$xx,-$yy)]}\n";
        }
    }
    if(/c /)
    {
        my $col=$cc[@trails];
        my $trail='';
        chomp;
        my @pts=split /[\[\]()] ?[\[\]()]/, "$'";
        shift @pts;
        pop @pts;
        my $p=-1;
        for(@pts)
        {
            /$d $d/;
            my $q=($2>0);
            if($p!=$q)
            {
                $p=$q;
                my $c=$p?$col:'invtri';
                $trail.="\n\\psline[linewidth=5px,linecolor=$c]{c-c}".&spt($_);
            }
            $trail.=&spt($_);
        }
        my $c=$p?$col:'invtri';
        push @trails, "$trail\n\\pscircle[fillcolor=$c,linestyle=none,fillstyle=solid]@{[&spt($pts[$#pts])]} {6}\n";
        if($tristretch)
        {
            $vertline = "\\psline[linecolor=arvelocity,fillcolor=arvelocity,fillstyle=solid,linewidth=2px]" . &pt($pts[$#pts],-$sr) . &pt($pts[$#pts],$sr) . "\n";
        }
    }
    if(/u $nd$d$nd$d$nd$d$nd$d\]/)
    {
        $contour{"$1,$2"}=1;
        $contour{"$3,$4"}=1;
    }
    if(/t /)
    {
        my @pts=split /[\[\]()] ?[\[\]()]/, "$'";
        shift @pts;
        pop @pts;
        for(my $i=0;$i+2<@pts;$i+=3)
        {
            my $col=$cc[@triangles];
            my ($a,$b,$c,$d,$e,$f)=((split ' ',$pts[$i]), (split ' ',$pts[$i+1]), (split ' ',$pts[$i+2]));
            $c-=$a;
            $d-=$b;
            $e-=$a;
            $f-=$b;
            if($c*$f<$d*$e){$col='invtri'}
            $triangle = "\\pspolygon[linecolor=black,fillcolor=$col,fillstyle=solid]";
            $triangle .= &spt($pts[$i]) . &spt($pts[$i+1]) . &spt($pts[$i+2]) . "\n";
            push @triangles, $triangle;
            if($tristretch)
            {
                my @p=split ' ',$pts[$i];
                #$tristretchstr .= "\\psline[linecolor=black,fillcolor=black,fillstyle=solid,linewidth=2px]{->}" . &pt(@p) . &pt($p[0], $p[1]+.3) . "\n";
                $arrow = "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{c->}" . &pt(@p) . &pt($p[0], $p[1]+.5) . "\n";
                my @q=split ' ',$pts[$i+1];
                $tristretchstr .= "\\psline[linecolor=coltriline,fillcolor=coltriline,fillstyle=solid,linewidth=1pt]" . &pt(3.3, $q[1]) . &pt(4.9, $q[1]) . "\n";
                $tristretchpts .= "\\pscircle[fillcolor=black,linestyle=none,fillstyle=solid]".&spt($pts[$i+0])."{3}\n";
                $tristretchpts .= "\\pscircle[fillcolor=black,linestyle=none,fillstyle=solid]".&spt($pts[$i+1])."{3}\n";
                $tristretchpts .= "\\pscircle[fillcolor=black,linestyle=none,fillstyle=solid]".&spt($pts[$i+2])."{3}\n";
            }
        }
    }
    if(/tristretch/){$tristretch=1;}
}

my $orig="\\pscircle[fillcolor=colcontour,linestyle=none,fillstyle=solid]@{[&pt(1,1)]} {11}\n";
my $contour = join '', map {&spt($_)} sort {$a=~/$d/;my $A=$1;$b=~/$d/;my $B=$1;$A<=>$B;} keys %contour;
$contour="\\psline[linewidth=5px,linecolor=colcontour]{c-c}$contour\n";

$"='';
my $wm1=$w-1;
my $hm1=$h-1;
print <<EOS;
\\documentclass{article}
\\usepackage{pstricks}
\\usepackage{color}

\\usepackage[margin=0cm,papersize={${wm1}px,${hm1}px}]{geometry}
\\definecolor{bg}{rgb}{0.75,0.75,0.765}
\\definecolor{coltri1}{rgb}{0,0.25,1}
\\definecolor{coltri2}{rgb}{0,0.9,0}
\\definecolor{coltri3}{rgb}{.95,0,0}
\\definecolor{invtri}{rgb}{.95,0,0}
\\definecolor{backtri}{rgb}{0.65,0.65,0.66}
\\definecolor{ltbacktri}{rgb}{0.75,0.75,0.765}
\\definecolor{vertline}{rgb}{0.575,0.575,0.585}
\\definecolor{colcontour}{rgb}{0.95,0.95,0}
\\definecolor{colarrow}{rgb}{1,0,1}
\\definecolor{colarrowref}{rgb}{0.5,0.5,0.51}
\\definecolor{coltriline}{rgb}{0.5,0.5,0.51}
\\definecolor{arvelocity}{rgb}{0,0.65,0}
\\begin{document}
\\noindent
\\psset{unit=1px}
\\begin{pspicture}(0,0)($wm1,$hm1)
\\psframe[fillcolor=bg,fillstyle=solid,linecolor=bg](0,0)($w,$h)
\\pspolygon[linecolor=backtri,fillcolor=backtri,fillstyle=solid]@{[&pt(-$sr,$sr)]}@{[&pt(-$sr,-$sr)]}@{[&pt($sr,-$sr)]}
\\pspolygon[linecolor=ltbacktri,fillcolor=ltbacktri,fillstyle=solid]@{[&pt(-$sr,$sr)]}@{[&pt($sr,$sr)]}@{[&pt($sr,-$sr)]}

$vertline
\\psline[linewidth=2px]{->}@{[&pt(-$sr,0)]}@{[&pt($sr,0)]}
\\psline[linewidth=2px]{->}@{[&pt(0,-$sr)]}@{[&pt(0,$sr)]}
$orig
$contour
$tristretchstr
$arrow
@triangles
$tristretchpts
@trails
$arrows

\\psframe[fillstyle=solid,linestyle=none,fillcolor=backtri]@{[&pt(-3,.02)]}@{[&pt(-2.5,.3)]}
\\psframe[fillstyle=solid,linestyle=none,fillcolor=backtri]@{[&pt(-.4,-3)]}@{[&pt(-.02,-2.7)]}
\\uput[ur]@{[&pt(-3,0)]} {{\\Huge\$\\sigma_1\$}}
\\uput[ul]@{[&pt(0,-3)]} {{\\Huge\$\\sigma_2\$}}

\\end{pspicture}
\\end{document}
EOS
