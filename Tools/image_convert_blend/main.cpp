//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    PARSE_ARGS args;
    args.Add_Double_Argument("-gamma",1,"-gamma","gamma correction");
    args.Add_Double_Argument("-ratio",0.5,"-ratio","blend ratio");
    args.Set_Extra_Arguments(3,"<image in1> <image in2> <image out>","images to read and write");
    args.Parse(argc,argv);
    std::string file_input1=args.Extra_Arg(1);
    std::string file_input2=args.Extra_Arg(2);
    std::string file_output=args.Extra_Arg(3);
    double gamma=args.Get_Double_Value("-gamma");
    double ratio=args.Get_Double_Value("-ratio");

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
