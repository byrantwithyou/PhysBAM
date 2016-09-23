//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN  
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
using namespace PhysBAM;

// for i in `seq 0 \`cat last_frame\``; do ../adjust_temperature $i; done
// for i in 48 100 172 350 450 466; do ../adjust_temperature $i; done

template<class T> inline T Adjust_Temperature(const T t,const int frame,const T y)
{
    //T a=39./20,b=-1./2400,c=1/3e7;
    //T a=1.7,b=-.00043333,c=6.6666e-8;
    T a=1.25,b=.0000166667,c=-3.3333333333e-8; // seq_1.png, f[3000] == 3000, f[1500] == 1800, f[2000] == 2300
    T t_fix=t*(a+t*(b+t*c));
    if(y>.2) return t_fix;
    else if(frame>100) return t_fix;
    else{
        T reduction=y<.15?0:(y-.15)/.05;
        if(frame<90) return reduction*t_fix;
        else{
            T tau=(T)(frame-90)/10;
            return (tau+(1-tau)*reduction)*t_fix;}}
}

template<class T> inline T Adjust_Density(const T d,const int frame,const VECTOR<T,3> X)
{
    T y=X.y,r=(X-VECTOR<T,3>(.5,y,.5)).Magnitude();
    if(y>.28 || frame>150) return d;
    else{
        T reduction=r>.1?0:1;
        if(frame<132) return reduction*d;
        else{
            T tau=(T)(frame-132)/18;
            return (tau+(1-tau)*reduction)*d;}}
}

template<class T,class RW> void Process(int argc,char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    int frame=-10;
    parse_args.Extra(&frame,"frame","frame");
    parse_args.Parse();

    LOG::cout<<"Adjust temperature and density for frame "<<frame<<std::endl;

    std::string prefix="";
    std::string f=LOG::sprintf(".%d",frame);

    GRID<TV> grid;FILE_UTILITIES::Read_From_File<RW>("grid",grid);
    ARRAY<T,VECTOR<int,3> > temperature,density;
    FILE_UTILITIES::Read_From_File<RW>("temperature"+f,temperature);
    FILE_UTILITIES::Read_From_File<RW>("density"+f,density);
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)temperature(i,j,ij)=Adjust_Temperature<T>(temperature(i,j,ij),frame,grid.y(j));
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)density(i,j,ij)=Adjust_Density<T>(density(i,j,ij),frame,grid.X(i,j,ij));
    FILE_UTILITIES::Write_To_File<RW>("adjusted_temperature"+f,temperature);
    FILE_UTILITIES::Write_To_File<RW>("adjusted_density"+f,density);
}

int main(int argc,char* argv[])
{
    Process<float,float>(argc,argv);
    return 0;
}
//#####################################################################

