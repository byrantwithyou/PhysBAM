//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Images/IMAGE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    std::string file_input1,file_input2,file_output;
    double gamma=1,ratio=0.5;
    PARSE_ARGS args;
    args.Add("-gamma",&gamma,"gamma","gamma correction");
    args.Add("-ratio",&ratio,"ratio","blend ratio");
    args.Extra(&file_input1,"image in1","image to read");
    args.Extra(&file_input2,"image in2","image to read");
    args.Extra(&file_output,"image out","image to write");
    args.Parse(argc,argv);

    ARRAYS<VECTOR<VECTOR<float,3> ,2> > image1;IMAGE<float>::Read(file_input1,image1);
    ARRAYS<VECTOR<VECTOR<float,3> ,2> > image2;IMAGE<float>::Read(file_input2,image2);
    
    if(image1.domain.min_corner.x!=image2.domain.min_corner.x || image1.domain.min_corner.y!=image2.domain.min_corner.y || image1.domain.max_corner.x!=image2.domain.max_corner.x || image1.domain.max_corner.y!=image2.domain.max_corner.y){
        printf("Image sizes are not the same. Exiting\n");exit(0);}
    for(int i=image1.domain.min_corner.x;i<=image1.domain.max_corner.x;++i){
      for(int j=image1.domain.min_corner.y;j<=image1.domain.max_corner.y;++j){
        VECTOR<float,3> out_vector=(float)ratio*image1(i,j)+(1.0f-(float)ratio)*image2(i,j);
        image2(i,j)=out_vector;}}
    IMAGE<float>::Write(file_output,image2,gamma);
}
