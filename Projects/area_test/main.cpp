//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRAPEZOID_INTERSECTION.h>
using namespace PhysBAM;

template<class TV>
bool Inside_Trapezoid(TV a,TV b,const TV& p)
{
    if(a.x>b.x) exchange(a,b);
    PHYSBAM_ASSERT(a.y>=0);
    PHYSBAM_ASSERT(b.y>=0);
    PHYSBAM_ASSERT(p.y>=0);

    if(p.x<a.x) return false;
    if(p.x>b.x) return false;
    if((-a.y*b.x-b.y*p.x+a.x*b.y+a.y*p.x+p.y*b.x-p.y*a.x)>0) return false;
    return true;
}

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,2> TV;

    RANDOM_NUMBERS<T> random;

    TV a,b,c,d;
    random.Fill_Uniform(a,0,1);
    random.Fill_Uniform(b,0,1);
    random.Fill_Uniform(c,0,1);
    random.Fill_Uniform(d,0,1);

    VECTOR<TV,4> dA;
    VECTOR<VECTOR<MATRIX<T,2>,4>,4> H;

    fprintf(stderr, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(stderr, "%%%%BoundingBox: 0 0 1000 1000\n");
    fprintf(stderr, "1 0 0 setrgbcolor %g %g moveto %g %g lineto stroke\n", a.x*1000, a.y*1000, b.x*1000, b.y*1000);
    fprintf(stderr, "0 0 1 setrgbcolor %g %g moveto %g %g lineto stroke\n", c.x*1000, c.y*1000, d.x*1000, d.y*1000);
    fprintf(stderr, "0 0 0 setrgbcolor %g 0 moveto %g 1000 lineto stroke\n", a.x*1000, a.x*1000);
    fprintf(stderr, "0 1 0 setrgbcolor %g 0 moveto %g 1000 lineto stroke\n", b.x*1000, b.x*1000);
    fprintf(stderr, "0 0 0 setrgbcolor %g 0 moveto %g 1000 lineto stroke\n", c.x*1000, c.x*1000);
    fprintf(stderr, "0 1 0 setrgbcolor %g 0 moveto %g 1000 lineto stroke\n", d.x*1000, d.x*1000);

    T A = Trapezoid_Intersection_Area(a,b,c,d,dA,H);
    printf("%g\n", A);

    int I=0,N=1000000;
    for(int i=1;i<=N;i++){
        TV p;
        random.Fill_Uniform(p,0,1);
        if(Inside_Trapezoid(a,b,p) && Inside_Trapezoid(c,d,p))
            I++;
    }

    printf("%g\n", (T)I/N);

    return 0;
}
