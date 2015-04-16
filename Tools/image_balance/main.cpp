//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

int main(int argc, char* argv[])
{
    bool annotate_input=false;
    std::string in_file,out_file,param_file;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-i",&in_file,"file","input filename");
    parse_args.Add("-o",&out_file,"file","output filename");
    parse_args.Add("-a",&annotate_input,"output annotated image");
    parse_args.Add("-p",&param_file,"file","parameter filename");
    parse_args.Parse();

    ARRAY<TV,VECTOR<int,2> > image;
    PNG_FILE<T>::Read(in_file,image);

    ARRAY<TV,VECTOR<int,2> > annotate_image=image;
    FILE* P=fopen(param_file.c_str(),"r");
    PHYSBAM_ASSERT(P);

    RANGE<VECTOR<int,2> > range;
    TV target;
    T weight;
    TV A,B,C;
    while(fscanf(P,"%d %d %d %d %lg %lg %lg %lg\n",
            &range.min_corner.x,&range.min_corner.y,&range.max_corner.x,&range.max_corner.y,
            &weight,&target.x,&target.y,&target.z)==8)
    {
        range=range.Intersect(image.domain);
        LOG::printf("RANGE: %P  WEIGHT: %P  TARGET: %P\n",range,weight,target);
        TV E;
        for(RANGE_ITERATOR<2> it(range);it.Valid();it.Next())
        {
            TV I=image(it.index);
            E+=I;
            A+=weight*I*I;
            B+=weight*I*target;
            C+=weight*target*target;
            annotate_image(it.index)=target;
        }
        E/=range.Size();
        LOG::printf("AVERAGE: %P\n",E);
    }

    fclose(P);

    TV s=B/A;
    LOG::printf("CORRECTION: %P\n",s);
    for(RANGE_ITERATOR<2> it(image.domain);it.Valid();it.Next())
        image(it.index)=s*image(it.index);

    if(annotate_input) PNG_FILE<T>::Write(out_file,annotate_image);
    else PNG_FILE<T>::Write(out_file,image);
    

    return 0;
}


