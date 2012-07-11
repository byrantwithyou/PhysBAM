//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
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
    PARSE_ARGS parse_args;
    parse_args.Set_Extra_Arguments(1, "<frame>");
    parse_args.Parse();

    int frame=-10;
    if(parse_args.Num_Extra_Args() >= 1) frame=atoi(parse_args.Extra_Arg(0).c_str());
    else{std::cout<<"Incorrect.\n";exit(1);}

    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);
    std::cout<<"frame = "<<frame<<", translation = "<<Translation_Vector<T>(frame)<<std::endl;

    std::ostream* output=FILE_UTILITIES::Safe_Open_Output("animated_translation"+f,false);    
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

