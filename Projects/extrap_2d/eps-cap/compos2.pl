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
my @invtriangles=();
my %contour=();
my $svddots='';
$"=",";
my $tritop=-100;
my $tribot=100;
my $numtri=0;
while(<>)
{
    if(0 && /p $nd$d$nd$d$nd$d$nd$d\]/)
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
        my $col="coltri".((scalar @trails)%8);
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
                $trail.="\n\\psline[linewidth=1px,linecolor=arvelocity]{c-c}".&spt($_);
            }
            if($1>3){next;}
            $trail.=&spt($_);
        }
        $pts[$#pts]=~/$d $d/;
        my $c=$p?'coltri0':'invtri';
        push @trails, "$trail\n";
        if($1>3){next;}
        $svddots.="\\pscircle[fillcolor=$c,linestyle=none,fillstyle=solid]@{[&spt($pts[$#pts])]} {2}\n";
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
            my $col="coltri".($numtri++%8);
            my ($a,$b,$c,$d,$e,$f)=((split ' ',$pts[$i]), (split ' ',$pts[$i+1]), (split ' ',$pts[$i+2]));
            my @r=($b+4.1,$a,$f+4.1,$e,$d+4.1,$c);
            $c-=$a;
            $d-=$b;
            $e-=$a;
            $f-=$b;
            if($c*$f<$d*$e){$col='invtri'}
            $triangle = "\\psline[linecolor=black,fillcolor=$col,fillstyle=solid]{c-c}";
            $triangle .= &pt(@r[0..1]) . &pt(@r[2..3]) . &pt(@r[4..5]) . &pt(@r[0..1]) . "\n";
            if($r[1]>$tritop){$tritop=$r[1];}
            if($r[1]<$tribot){$tribot=$r[1];}
            if($r[3]>$tritop){$tritop=$r[3];}
            if($r[3]<$tribot){$tribot=$r[3];}
            if($r[5]>$tritop){$tritop=$r[5];}
            if($r[5]<$tribot){$tribot=$r[5];}
            if($c*$f<$d*$e){push @invtriangles, $triangle;}
            else{push @triangles, $triangle;}
        }
    }
}

my $orig="\\pscircle[fillcolor=colcontour,linestyle=none,fillstyle=solid]@{[&pt(1,1)]} {11}\n";
my $contour = join '', map {&spt($_)} sort {$a=~/$d/;my $A=$1;$b=~/$d/;my $B=$1;$A<=>$B;} keys %contour;
$contour="\\psline[linewidth=5px,linecolor=colcontour]{c-c}$contour\n";
my $extratrilines="\\psline[linewidth=5px,linecolor=black]{c-c}" . &pt(3.2,$tribot) . &pt(5.0,$tribot) . "\n";
$extratrilines .= "\\psline[linewidth=5px,linecolor=black]{c-c}" . &pt(3.2,$tritop) . &pt(5.0,$tritop) . "\n";

$"='';
my $wm1=$w-1;
my $hm1=$h-1;
print <<EOS;
\\documentclass{article}
\\usepackage{pstricks}
\\usepackage{color}

\\usepackage[margin=0cm,papersize={${wm1}px,${hm1}px}]{geometry}
\\definecolor{bg}{rgb}{0.75,0.75,0.765}
\\definecolor{coltri0}{rgb}{0,0.25,1}
\\definecolor{coltri1}{rgb}{0,0.9,0}
\\definecolor{coltri2}{rgb}{0,1,1}
\\definecolor{coltri3}{rgb}{1,1,0}
\\definecolor{coltri4}{rgb}{1,1,0}
\\definecolor{coltri5}{rgb}{0,1,1}
\\definecolor{coltri6}{rgb}{0,0.9,0}
\\definecolor{coltri7}{rgb}{0,0.25,1}
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
$arrow
@triangles
@invtriangles
@trails
$arrows
$svddots
$extratrilines

\\end{pspicture}
\\end{document}
EOS
