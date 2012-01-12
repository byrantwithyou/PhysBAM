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

$ARGV[0]=~/$d/;
my $frame=$&+0;

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
my $velarrows='';
my $vertline='';
my @trails=();
my %contour=();
my $svddots='';
$"=",";
my $tritop=-100;
my $tribot=100;
my $numtri=0;
my $originrad=11;
my $quasidot='';
my $sv=0;
while(<STDIN>)
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
        $svddots.="\\pscircle[fillcolor=$c,linestyle=none,fillstyle=solid]@{[&spt($pts[$#pts])]} {6}\n";
    }
    if(/u $nd$d$nd$d$nd$d$nd$d\]/)
    {
        $contour{"$1,$2"}=1;
        $contour{"$3,$4"}=1;
    }
    if(/s $d/)
    {
        $sv=$1;
        my $r=($1-1)*$s;
        if($r>=$originrad){$originrad=0;}
        else{$originrad=sqrt($originrad*$originrad-$r*$r);}
    }
    if(/m $d$nd$d/)
    {
        my $a=$1/(2*$2);
        my $x=(1-$a*$sv+3*$a)/(1+2*$a);
        my $par=$frame&1;
        if($x<0 && !$par){$x=-$x;}
        $quasidot="\\pscircle[fillcolor=coltri0,linestyle=none,fillstyle=solid]@{[&pt($x,$x)]} {6}\n";
    }
}

my $orig="\\pscircle[fillcolor=colcontour,linestyle=none,fillstyle=solid]@{[&pt(1,1)]} {$originrad}\n";
my $contour = join '', map {&spt($_)} sort {$a=~/$d/;my $A=$1;$b=~/$d/;my $B=$1;$A<=>$B;} keys %contour;
my $contcolor=$frame>=40?'grcolcontour':'colcontour';
$contour="\\psline[linewidth=5px,linecolor=$contcolor]{c-c}$contour\n";
my $extratrilines="\\psline[linewidth=5px,linecolor=black]{c-c}" . &pt(3.2,$tribot) . &pt(5.0,$tribot) . "\n";
$extratrilines .= "\\psline[linewidth=5px,linecolor=black]{c-c}" . &pt(3.2,$tritop) . &pt(5.0,$tritop) . "\n";
$velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(4.1,$tritop) . &pt(4.1, $tritop+.5) . "\n";
$velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(4.7,$tritop) . &pt(4.7, $tritop+.5) . "\n";
$velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(3.5,$tritop) . &pt(3.5, $tritop+.5) . "\n";
$velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(4.1,$tribot) . &pt(4.1, $tribot-.5) . "\n";
$velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(4.7,$tribot) . &pt(4.7, $tribot-.5) . "\n";
$velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(3.5,$tribot) . &pt(3.5, $tribot-.5) . "\n";

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
\\definecolor{grcolcontour}{rgb}{0.5,0.5,0.51}
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
$quasidot
$arrows
\\psframe[fillstyle=solid,linestyle=none,fillcolor=backtri]@{[&pt(-3,.02)]}@{[&pt(-2.5,.3)]}
\\psframe[fillstyle=solid,linestyle=none,fillcolor=backtri]@{[&pt(-.4,-3)]}@{[&pt(-.02,-2.7)]}
\\uput[ur]@{[&pt(-3,0)]} {{\\Huge\$\\sigma_1\$}}
\\uput[ul]@{[&pt(0,-3)]} {{\\Huge\$\\sigma_2\$}}

\\end{pspicture}
\\end{document}
EOS
