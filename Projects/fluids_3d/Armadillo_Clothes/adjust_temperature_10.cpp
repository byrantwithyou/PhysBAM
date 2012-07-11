//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN  
//##################################################################### 
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <Heat_Flows/HEAT_3D.h>
using namespace PhysBAM;

// for i in `seq 0 \`cat last_frame\``; do ../adjust_temperature_10 $i; done

template<class T> inline T Spline(const T a,const T b,const T t)
{
    T s=(t-a)/(b-a);return s<0?0:s>1?1:s*s*(3-2*s);
}

template<class T> inline T Blur_Temperature(const T temperature,const T phi)
{
    //T t=3000*Spline((T)-3,(T)0,phi);
    //T w=Spline((T)-3,(T)-1,phi);
    T t=3000*Spline((T)-3,(T)0,phi);
    T w=Spline((T)-3,(T)-.5,phi);
    return w*max(temperature,t)+(1-w)*temperature;
}

template<class T> inline T Adjust_Temperature(const T t,const int frame,const VECTOR<T,3> X)
{
    //T a=39./20,b=-1./2400,c=1/3e7;
    //T a=1.7,b=-.00043333,c=6.6666e-8;
    T a=1.25,b=.0000166667,c=-3.3333333333e-8; // seq_1.png, f[3000] == 3000, f[1500] == 1800, f[2000] == 2300
    T t_fix=t*(a+t*(b+t*c));
    T y=X.y,r=(X-VECTOR<T,3>(.5,y,.34)).Magnitude();
    if(y>.08 || frame>110) return t_fix;
    else{
        T reduction=Spline((T).1,(T).075,r);
        if(frame<90) return reduction*t_fix;
        else{
            T tau=(T)(frame-90)/20;
            return (tau+(1-tau)*reduction)*t_fix;}}
}

template<class T> inline T Adjust_Density(const T d,const int frame,const VECTOR<T,3> X)
{
    T y=X.y,r=(X-VECTOR<T,3>(.5,y,.34)).Magnitude();
    if(y>.08 || frame>110) return d;
    else{
        T reduction=Spline((T).1,(T).075,r);
        if(frame<90) return reduction*d;
        else{
            T tau=(T)(frame-90)/20;
            return (tau+(1-tau)*reduction)*d;}}
}

template<class T,class RW> void Process(int argc,char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Set_Extra_Arguments(1, "<frame>");
    parse_args.Parse();

    int frame=-10;
    if(parse_args.Num_Extra_Args() >= 1) frame=atoi(parse_args.Extra_Arg(0).c_str());
    else{std::cout<<"Incorrect.\n";exit(1);}

    std::cout<<"Adjust temperature and density for frame "<<frame<<std::endl;

    std::string prefix="";
    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);

    ARRAY<T,VECTOR<int,3> > temperature,density,phi;
    GRID<TV> grid;LEVELSET_3D<GRID<TV> > levelset(grid,phi);
    FILE_UTILITIES::Read_From_File<RW>("levelset"+f,levelset);
    FILE_UTILITIES::Read_From_File<RW>("temperature"+f,temperature);
    FILE_UTILITIES::Read_From_File<RW>("density"+f,density);

    HEAT_3D<T>::Smooth(grid,phi,2);
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)temperature(i,j,ij)=Adjust_Temperature<T>(Blur_Temperature(temperature(i,j,ij),phi(i,j,ij)*grid.one_over_dx),frame,grid.X(i,j,ij));
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)density(i,j,ij)=Adjust_Density<T>(density(i,j,ij),frame,grid.X(i,j,ij));
    FILE_UTILITIES::Write_To_File<RW>("adjusted_temperature"+f,temperature);
    FILE_UTILITIES::Write_To_File<RW>("adjusted_density"+f,density);
    FILE_UTILITIES::Write_To_File<RW>("adjusted_levelset"+f,levelset);
}

int main(int argc,char* argv[])
{
    Process<float,float>(argc,argv);
    return 0;
}
//#####################################################################

