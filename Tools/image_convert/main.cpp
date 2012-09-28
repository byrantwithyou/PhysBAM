//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
using namespace PhysBAM;

template<class T> void
Bloom(ARRAY<VECTOR<T,3>,VECTOR<int,2> >& image,const T bloom_radius,const T bloom_weight)
{
    LOG::SCOPE("bloom","Bloom Filter");
    int bloom_support=(int)(bloom_radius*max(image.counts.x,image.counts.y)),bloom_width=bloom_support/2,bloom_width_squared=sqr(bloom_width);
    //LOG::cout<<"bloom_support "<<bloom_support<<" width "<<bloom_width<<" squared="<<bloom_width_squared<<std::endl;
    LOG::Time("Building bloom filter table");
    ARRAY<T,VECTOR<int,1> > bloom_table(0,sqr(bloom_width)-1);
    for(int i=0;i<=bloom_table.counts.x;i++){T distance=sqrt((T)i)/bloom_width;bloom_table(i)=pow(max((T)0,(T)1-distance),(T)4);}
    //LOG::cout<<"table "<<bloom_table<<std::endl;
    LOG::Time("Filtering");
    ARRAY<VECTOR<T,3>,VECTOR<int,2> > bloom_image(1,image.counts.x,1,image.counts.y);
    PROGRESS_INDICATOR progress(image.counts.x*image.counts.y);
    for(int i=0;i<image.counts.x;i++) for(int j=0;j<image.counts.y;j++){
        T weight_sum=0;
        RANGE<VECTOR<int,2> > bloom_box(clamp_min(i-bloom_width,1),clamp_max(i+bloom_width,image.counts.x),clamp_min(j-bloom_width,1),clamp_max(j+bloom_width,image.counts.y));
        //std::cout<<"pixel "<<i<<","<<j<<" box="<<bloom_box<<std::endl;
        for(int ii=bloom_box.min_corner.x;ii<=bloom_box.max_corner.x;ii++) for(int jj=bloom_box.min_corner.y;jj<=bloom_box.max_corner.y;jj++) if(ii!=i || jj!=j){
            int dx=ii-i,dy=jj-j;
            int distance_squared=sqr(dx)+sqr(dy);
            if(distance_squared<bloom_width_squared){
                T weight=bloom_table(distance_squared);
                weight_sum+=weight;
                bloom_image(i,j)+=(T)weight*image(ii,jj);}}
        bloom_image(i,j)/=weight_sum;
        progress.Progress();}
    IMAGE<T>::Write("bloom.png",bloom_image,(T)1);
    LOG::Time("Blending");
    for(int i=0;i<image.counts.x;i++) for(int j=0;j<image.counts.y;j++) image(i,j)=(T)bloom_weight*bloom_image(i,j)+(1-bloom_weight)*image(i,j);
}

int main(int argc,char* argv[])
{
    float gamma=1,bloom_radius=0,bloom_weight=.5;
    bool use_bloom_radius=false;
    std::string file_input,file_output;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-gamma",&gamma,"-gamma","gamma correction");
    parse_args.Add("-bloom_weight",&bloom_weight,"-bloom_weight","bloom weight");
    parse_args.Add("-bloom_radius",&bloom_radius,&use_bloom_radius,"-bloom_radius","bloom radius"); 
    parse_args.Extra(&file_input,"image in","image to read");
    parse_args.Extra(&file_output,"image out","image to write");
    parse_args.Parse();
    
    ARRAY<VECTOR<float,3>,VECTOR<int,2> > image;
    LOG::Time("Reading Image");
    IMAGE<float>::Read(file_input,image);
    if(use_bloom_radius)
        Bloom(image,bloom_radius,bloom_weight);
    LOG::Time("Writing Image");
    IMAGE<float>::Write(file_output,image,(float)gamma);
}
