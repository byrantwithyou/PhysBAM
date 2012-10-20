//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Symbolics/MULTIVARIATE_POLYNOMIAL.h>

using namespace PhysBAM;

typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> TV_INT;
typedef VECTOR<T,2> TV2;
typedef VECTOR<int,2> TV_INT2;

template<class TV>
void Integral_Test()
{
    typedef typename TV::SCALAR T;
    RANDOM_NUMBERS<T> random;
    VECTOR<TV,3> pts;
    random.Fill_Uniform(pts(0),-(T)1,(T)1);
    random.Fill_Uniform(pts(1),-(T)1,(T)1);
    random.Fill_Uniform(pts(2),-(T)1,(T)1);

    MULTIVARIATE_POLYNOMIAL<TV> P;
    for(int i=0;i<10;i++){
        TV power;
        random.Fill_Uniform(power,0,(T)5);
        P.terms.Append(MULTIVARIATE_MONOMIAL<TV>(VECTOR<int,TV::m>(power),random.Get_Uniform_Number(-1,1)));}
    P.Simplify();
    std::cout<<P<<std::endl;

    int N=1000000;
    T i0=P.Integrate_Over_Primitive(pts),i1=0;
    for(int i=0;i<N;i++){
        VECTOR<T,2> uv;
        random.Fill_Uniform(uv,0,1);
        if(uv.Sum()>1) uv=(T)1-uv;
        TV pt=pts(0)+(pts(1)-pts(0))*uv.x+(pts(2)-pts(0))*uv.y;
        i1+=P.Value(pt);}
    i1*=TV::Cross_Product(pts(1)-pts(0),pts(2)-pts(0)).Magnitude()/(2*N);
    LOG::cout<<"integrals c "<<i0<<"  r "<<i1<<"  e "<<abs(i1-i0)<<std::endl;
}

int main(int argc, char* argv[])
{
    LOG::cout<<"### LINEAR POLYNOMIAL OVER A TRIANGLE"<<std::endl;
    {
        RANDOM_NUMBERS<T> random;
        MULTIVARIATE_POLYNOMIAL<TV> F;
        
        F.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(1,0,0),random.Get_Uniform_Number(-1,1)));
        F.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(0,1,0),random.Get_Uniform_Number(-1,1)));
        F.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(0,0,1),random.Get_Uniform_Number(-1,1)));
        F.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(0,0,0),random.Get_Uniform_Number(-1,1)));
        
        VECTOR<TV,3> triangle;
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) triangle(i)(j)=random.Get_Uniform_Number(-1,1);
        TV center=triangle.Sum()/3;
        
        LOG::cout<<"triangle "<<triangle<<std::endl;
        LOG::cout<<"center "<<center<<std::endl;
        LOG::cout<<"poly "<<F<<std::endl;
        LOG::cout<<"integral "<<F.Integrate_Over_Primitive(triangle)<<std::endl;
        LOG::cout<<"averaged "<<(center(0)*F.terms(0).coeff+center(1)*F.terms(1).coeff+center(2)*F.terms(2).coeff+F.terms(3).coeff)*
            TV::Cross_Product(triangle(1)-triangle(0),triangle(2)-triangle(0)).Magnitude()/2<<std::endl;
    }
    
    LOG::cout<<"### LINEAR POLYNOMIAL OVER AN INTERVAL"<<std::endl;
    {
        RANDOM_NUMBERS<T> random;
        MULTIVARIATE_POLYNOMIAL<TV> F;
        
        F.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(1,0,0),random.Get_Uniform_Number(-1,1)));
        F.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(0,1,0),random.Get_Uniform_Number(-1,1)));
        F.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(0,0,1),random.Get_Uniform_Number(-1,1)));
        F.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(0,0,0),random.Get_Uniform_Number(-1,1)));
        
        VECTOR<TV,2> interval;
        for(int i=0;i<2;i++) for(int j=0;j<3;j++) interval(i)(j)=random.Get_Uniform_Number(-10,10);
        TV center=interval.Sum()/2;
        
        LOG::cout<<"interval "<<interval<<std::endl;
        LOG::cout<<"center "<<center<<std::endl;
        LOG::cout<<"poly "<<F<<std::endl;
        LOG::cout<<"integral "<<F.Integrate_Over_Primitive(interval)<<std::endl;
        LOG::cout<<"averaged "<<(center(0)*F.terms(0).coeff+center(1)*F.terms(1).coeff+center(2)*F.terms(2).coeff+F.terms(3).coeff)*
            (interval(1)-interval(0)).Magnitude()<<std::endl;
            }
    
    LOG::cout<<"### POLYNOMIAL OVER A BOX"<<std::endl;
    {
        RANDOM_NUMBERS<T> random;
        MULTIVARIATE_POLYNOMIAL<TV> F;
        
        for(int ix=0;ix<5;ix++)
        for(int iy=0;iy<5;iy++)
        for(int iz=0;iz<5;iz++)
            if(ix+iy+iz<=5) F.terms.Append(MULTIVARIATE_MONOMIAL<TV>(TV_INT(ix,iy,iz),random.Get_Uniform_Number(-1,1)));

        TV min_corner(random.Get_Uniform_Number(-1,1),random.Get_Uniform_Number(-1,1),random.Get_Uniform_Number(-1,1));
        TV max_corner(random.Get_Uniform_Number(0,1),random.Get_Uniform_Number(0,1),random.Get_Uniform_Number(0,1));
        max_corner+=min_corner;
        RANGE<TV> range(min_corner,max_corner);

        LOG::cout<<"range "<<range<<std::endl;
        LOG::cout<<"poly "<<F<<std::endl;
        LOG::cout<<"integral   "<<F.Definite_Integral(range)<<std::endl;
        F.Integrate(1);

        VECTOR<TV,3> triangle1;
        VECTOR<TV,3> triangle2;
        VECTOR<TV,3> triangle3;
        VECTOR<TV,3> triangle4;
        
        triangle1(0)=TV(min_corner(0),max_corner(1),min_corner(2));
        triangle1(1)=TV(min_corner(0),max_corner(1),max_corner(2));
        triangle1(2)=TV(max_corner(0),max_corner(1),max_corner(2));

        triangle2(0)=TV(min_corner(0),max_corner(1),min_corner(2));
        triangle2(1)=TV(max_corner(0),max_corner(1),max_corner(2));
        triangle2(2)=TV(max_corner(0),max_corner(1),min_corner(2));

        triangle3(0)=TV(min_corner(0),min_corner(1),min_corner(2));
        triangle3(1)=TV(min_corner(0),min_corner(1),max_corner(2));
        triangle3(2)=TV(max_corner(0),min_corner(1),max_corner(2));

        triangle4(0)=TV(min_corner(0),min_corner(1),min_corner(2));
        triangle4(1)=TV(max_corner(0),min_corner(1),max_corner(2));
        triangle4(2)=TV(max_corner(0),min_corner(1),min_corner(2));

        LOG::cout<<"divergence "<<
            F.Integrate_Over_Primitive(triangle1)+
            F.Integrate_Over_Primitive(triangle2)-
            F.Integrate_Over_Primitive(triangle3)-
            F.Integrate_Over_Primitive(triangle4)<<std::endl;
    }

    LOG::cout<<"### POLYNOMIAL OVER A BOX 2D"<<std::endl;
    {
        RANDOM_NUMBERS<T> random;
        MULTIVARIATE_POLYNOMIAL<TV2> F;
        
        for(int ix=0;ix<5;ix++)
        for(int iy=0;iy<5;iy++)
            if(ix+iy<=5) F.terms.Append(MULTIVARIATE_MONOMIAL<TV2>(TV_INT2(ix,iy),random.Get_Uniform_Number(-1,1)));

        TV2 min_corner(random.Get_Uniform_Number(-1,1),random.Get_Uniform_Number(-1,1));
        TV2 max_corner(random.Get_Uniform_Number(0,1),random.Get_Uniform_Number(0,1));
        max_corner+=min_corner;
        RANGE<TV2> range(min_corner,max_corner);

        LOG::cout<<"range "<<range<<std::endl;
        LOG::cout<<"poly "<<F<<std::endl;
        LOG::cout<<"integral   "<<F.Definite_Integral(range)<<std::endl;
        F.Integrate(1);

        VECTOR<TV2,2> interval1;
        VECTOR<TV2,2> interval2;
        
        interval1(0)=TV2(min_corner(0),max_corner(1));
        interval1(1)=TV2(max_corner(0),max_corner(1));

        interval2(0)=TV2(min_corner(0),min_corner(1));
        interval2(1)=TV2(max_corner(0),min_corner(1));

        LOG::cout<<"divergence "<<
            F.Integrate_Over_Primitive(interval1)-
            F.Integrate_Over_Primitive(interval2)<<std::endl;
    }

    Integral_Test<VECTOR<T,2> >();
    Integral_Test<VECTOR<T,3> >();
    
    return 0;
}
