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

// for i in `seq 0 \`cat last_frame\``; do ../animate_translation_7 $i; done

template<class T> VECTOR<T,3> Translation_Vector(const int frame)
{
    T y;
    if(frame<70) y=0;
    else y=(T)-.09*(frame-70)/30;
    return VECTOR<T,3>(0,y,0);
}

template<class T,class RW> void Process(int argc,char* argv[])
{
    int frame=-10;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Extra(&frame,"frame","frame");
    parse_args.Parse();

    std::string f=LOG::sprintf(".%d",frame);
    std::cout<<"frame = "<<frame<<", translation = "<<Translation_Vector<T>(frame)<<std::endl;

    std::ostream* output=Safe_Open_Output("animated_translation"+f,false);    
    *output<<"Transform{\n"
           <<"  Type=Translate\n"
           <<"  Vector="<<Translation_Vector<T>(frame)<<"\n"
           <<"}\n";
    delete output;
}

int main(int argc,char* argv[])
{
    Process<float,float>(argc,argv);
    return 0;
}
//#####################################################################

