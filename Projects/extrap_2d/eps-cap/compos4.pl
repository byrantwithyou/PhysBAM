#!/usr/bin/perl -w
require strict;

my $mg=10;
my $h=480;
my $w=640;
my $sr=3;
my $s=($h-2*$mg)/(2*$sr);
my $shx=$mg + $sr * $s;
my $shy=$shx;
my $arlen=.135;
my $dispx=160;
my $firstdisp=40;

my $d="([0-9.eE+-]+)";
my $nd="[^0-9.eE+-]+";

sub pt
{
    my $x=$_[0] * $s + $shx;
    my $y=$_[1] * $s + $shy;
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
my $alltri='';
my $alltrails='';
my @invtriangles=();
my %contour=();
my $svddots='';
$"=",";
my $tritop=-100;
my $tribot=100;
my $numtri=0;
$"="";
my $extratrilines='';

my $dir=$ARGV[0];
shift @ARGV;
my $number=@ARGV;

my $tristretch=0;
for $i (0..$#ARGV)
{
    my @triangles=();
    my @trails=();
    my $last=($i==$#ARGV);
    open F, "$dir/$ARGV[$i]";

    while(<F>)
    {
        if(/p $nd$d$nd$d$nd$d$nd$d\]/ && $last)
        {
            if($1 + $2 > 0)
            {
                my $xx=$3*$arlen+$1;
                my $yy=$4*$arlen+$2;
                $arrows.="\\psline[linewidth=.5px,linecolor=colarrow]{->}@{[&pt($1,$2)]}@{[&pt($xx,$yy)]}\n";
                $arrows.="\\psline[linewidth=.5px,linecolor=colarrowref]{->}@{[&pt(-$1,-$2)]}@{[&pt(-$xx,-$yy)]}\n";
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
                    if($last)
                    {
                        $trail.="\n\\psline[linewidth=1px,linecolor=arvelocity]{c-c}".&spt($_);
                    }
                }
                if($1>3){next;}
                if($last)
                {
                    $trail.=&spt($_);
                }
            }
            $pts[$#pts]=~/$d $d/;
            my $c=$p?'coltri0':'invtri';
            push @trails, "$trail\n";
            if($1>3){next;}
            if($last)
            {
                $svddots.="\\pscircle[fillcolor=$c,linestyle=none,fillstyle=solid]@{[&spt($pts[$#pts])]} {3}\n";
            }
        }
        if(/u $nd$d$nd$d$nd$d$nd$d\]/ && $i==$#ARGV)
        {
            $contour{"$1,$2"}=1;
            $contour{"$3,$4"}=1;
        }
        if(/t /)
        {
            $shx = $shy + $i * $dispx + $firstdisp;
            my @pts=split /[\[\]()] ?[\[\]()]/, "$'";
            shift @pts;
            pop @pts;
            for(my $i=0;$i+2<@pts;$i+=3)
            {
                my $col="coltri0";#"coltri".($numtri++%8);
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

            $velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(4.1,$tritop) . &pt(4.1, $tritop+.5) . "\n";
            $velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(4.7,$tritop) . &pt(4.7, $tritop+.5) . "\n";
            $velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(3.5,$tritop) . &pt(3.5, $tritop+.5) . "\n";
            $velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(4.1,$tribot) . &pt(4.1, $tribot-.5) . "\n";
            $velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(4.7,$tribot) . &pt(4.7, $tribot-.5) . "\n";
            $velarrows .= "\\psline[linecolor=arvelocity,fillcolor=arvelocity,arrowinset=0,arrowlength=0.8,fillstyle=solid,linewidth=13px]{->}" . &pt(3.5,$tribot) . &pt(3.5, $tribot-.5) . "\n";
            $extratrilines .= "\\psline[linewidth=5px,linecolor=black]{c-c}" . &pt(3.2,$tribot) . &pt(5.0,$tribot) . "\n";
            $extratrilines .= "\\psline[linewidth=5px,linecolor=black]{c-c}" . &pt(3.2,$tritop) . &pt(5.0,$tritop) . "\n";

            $shx = $shy;
        }
    }
    close F;

    $alltri.="@triangles";
    $alltrails.="@trails";
}

my $orig="\\pscircle[fillcolor=colcontour,linestyle=none,fillstyle=solid]@{[&pt(1,1)]} {11}\n";
my $contour = join '', map {&spt($_)} sort {$a=~/$d/;my $A=$1;$b=~/$d/;my $B=$1;$A<=>$B;} keys %contour;
$contour="\\psline[linewidth=5px,linecolor=colcontour]{c-c}$contour\n";

$"='';
my $H=$h;
my $W=$w+$dispx*($number-1)+$firstdisp;
my $wm1=$W-1;
my $hm1=$H-1;
print <<EOS;
\\documentclass{article}
\\usepackage{pstricks}
\\usepackage{color}

\\usepackage[margin=0cm,papersize={${wm1}px,${hm1}px}]{geometry}
\\definecolor{bg}{rgb}{1,1,1}
\\definecolor{coltri0}{rgb}{0,0.25,1}
\\definecolor{coltri1}{rgb}{0,0.9,0}
\\definecolor{coltri2}{rgb}{0,1,1}
\\definecolor{coltri3}{rgb}{1,1,0}
\\definecolor{coltri4}{rgb}{1,1,0}
\\definecolor{coltri5}{rgb}{0,1,1}
\\definecolor{coltri6}{rgb}{0,0.9,0}
\\definecolor{coltri7}{rgb}{0,0.25,1}
\\definecolor{invtri}{rgb}{1,0,0}
\\definecolor{backtri}{rgb}{0.90,0.90,0.90}
\\definecolor{ltbacktri}{rgb}{1,1,1}
\\definecolor{vertline}{rgb}{0.575,0.575,0.585}
\\definecolor{colcontour}{rgb}{.8,.8,0}
\\definecolor{colarrow}{rgb}{1,.3,0}
\\definecolor{colarrowref}{rgb}{0.5,0.5,0.51}
\\definecolor{coltriline}{rgb}{0.5,0.5,0.51}
\\definecolor{arvelocity}{rgb}{0,0.65,0}
\\begin{document}
\\noindent
\\psset{unit=1px}
\\begin{pspicture}(0,0)($wm1,$hm1)
\\psframe[fillcolor=bg,fillstyle=solid,linecolor=bg](0,0)($W,$H)
\\pspolygon[linecolor=backtri,fillcolor=backtri,fillstyle=solid]@{[&pt(-$sr,$sr)]}@{[&pt(-$sr,-$sr)]}@{[&pt($sr,-$sr)]}
\\pspolygon[linecolor=ltbacktri,fillcolor=ltbacktri,fillstyle=solid]@{[&pt(-$sr,$sr)]}@{[&pt($sr,$sr)]}@{[&pt($sr,-$sr)]}

$vertline
\\psline[linewidth=2px]{->}@{[&pt(-$sr,0)]}@{[&pt($sr,0)]}
\\psline[linewidth=2px]{->}@{[&pt(0,-$sr)]}@{[&pt(0,$sr)]}
$orig
$contour
$alltri
@invtriangles
$alltrails
$arrows
$svddots
$velarrows
$extratrilines
$arrows
\\psframe[fillstyle=solid,linestyle=none,fillcolor=backtri]@{[&pt(-3,.02)]}@{[&pt(-2.5,.3)]}
\\psframe[fillstyle=solid,linestyle=none,fillcolor=backtri]@{[&pt(-.4,-3)]}@{[&pt(-.02,-2.7)]}
\\uput[ur]@{[&pt(-3,0)]} {{\\Huge\$\\sigma_1\$}}
\\uput[ul]@{[&pt(0,-3)]} {{\\Huge\$\\sigma_2\$}}

\\end{pspicture}
\\end{document}
EOS
