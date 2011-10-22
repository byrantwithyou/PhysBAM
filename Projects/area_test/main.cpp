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

RANDOM_NUMBERS<double> rn;

bool Test()
{
    typedef double T;
    typedef VECTOR<double,2> TV;
    VECTOR<TV,4> a;
    rn.Fill_Uniform(a(1),0,1);
    rn.Fill_Uniform(a(2),0,1);
    rn.Fill_Uniform(a(3),0,1);
    rn.Fill_Uniform(a(4),0,1);

    T e=1e-4;
    VECTOR<TV,4> da;
    rn.Fill_Uniform(da(1),-e,e);
    rn.Fill_Uniform(da(2),-e,e);
    rn.Fill_Uniform(da(3),-e,e);
    rn.Fill_Uniform(da(4),-e,e);

    fprintf(stderr, "1 1 1 setrgbcolor newpath 0 0 moveto 1000 0 lineto 1000 1000 lineto 0 1000 lineto closepath fill\n");
    VECTOR<TV,4> dA1;
    VECTOR<VECTOR<MATRIX<T,2>,4>,4> H1;
    T A1 = Trapezoid_Intersection_Area(a(1),a(2),a(3),a(4),dA1,H1);

    VECTOR<TV,4> dA2;
    VECTOR<VECTOR<MATRIX<T,2>,4>,4> H2;
    T A2 = Trapezoid_Intersection_Area(a(1)+da(1),a(2)+da(2),a(3)+da(3),a(4)+da(4),dA2,H2);

    VECTOR<T,8>& Va=(VECTOR<T,8>&)da;
    VECTOR<T,8>& V1=(VECTOR<T,8>&)dA1;
    VECTOR<T,8>& V2=(VECTOR<T,8>&)dA2;

    T dA=Va.Dot_Product(Va,(V1+V2)/(T)2);
    T aa=dA/e;
    T bb=(A2-A1)/e;
    if((!aa != !bb) || fabs((aa-bb)/bb)>1e-6) printf("AG %g %g %g\n", aa, bb, fabs((aa-bb)/bb));

    MATRIX<T,8> M1,M2;
    for(int i=1;i<=8;i++) for(int j=1;j<=8;j++) M1(i,j)=H1((i+1)/2)((j+1)/2)((i+1)%2+1,(j+1)%2+1);
    for(int i=1;i<=8;i++) for(int j=1;j<=8;j++) M2(i,j)=H2((i+1)/2)((j+1)/2)((i+1)%2+1,(j+1)%2+1);

    VECTOR<T,8> dG1=(M1+M2)/(T)2*Va;
    VECTOR<T,8> dG2=V2-V1;
    T cc=dG2.Magnitude()/e;
    T dd=(dG2-dG1).Magnitude()/e;
    if((!cc != !dd) || dd/cc>1e-6) printf("AH %g %g %g\n", cc, dd, dd/cc);

/*
    T dG=0;
    for(int i=1;i<=4;i++) dG+=TV::Dot_Product(da(i),(dA1(i)+dA2(i))/2);
    T aa=dA/e;
    T bb=(A2-A1)/e;
    if((!aa != !bb) || fabs((aa-bb)/bb)>1e-6) printf("AG %g %g %g\n", aa, bb, fabs((aa-bb)/bb));
*/

/*
    printf("%g %g  (%g)\n", A, A2, A2-A);
    fprintf(stderr, "1 0 0 setrgbcolor %g %g moveto %g %g lineto stroke\n", a.x*1000, a.y*1000, b.x*1000, b.y*1000);
    fprintf(stderr, "0 0 1 setrgbcolor %g %g moveto %g %g lineto stroke\n", c.x*1000, c.y*1000, d.x*1000, d.y*1000);
    fprintf(stderr, "0 0 0 setrgbcolor %g 0 moveto %g 1000 lineto stroke\n", a.x*1000, a.x*1000);
    fprintf(stderr, "0 1 0 setrgbcolor %g 0 moveto %g 1000 lineto stroke\n", b.x*1000, b.x*1000);
    fprintf(stderr, "0 0 0 setrgbcolor %g 0 moveto %g 1000 lineto stroke\n", c.x*1000, c.x*1000);
    fprintf(stderr, "0 1 0 setrgbcolor %g 0 moveto %g 1000 lineto stroke\n", d.x*1000, d.x*1000);
*/
    return true;
}

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    typedef VECTOR<T,2> TV;


    fprintf(stderr, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(stderr, "%%%%BoundingBox: 0 0 1000 1000\n");


    for(int k=0;k<1000;k++){
        Test();
    }

    return 0;
}
