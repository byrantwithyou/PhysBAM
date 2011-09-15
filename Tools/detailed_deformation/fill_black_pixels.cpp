//#####################################################################
// Copyright 2005, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <cstring>
#include <fstream>

using namespace PhysBAM;
using namespace std;

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{

    if (argc!=2){cout<<"Usage: fill_black_pixels image"<<::endl;exit(-1);}
 
    ARRAYS<VECTOR<VECTOR_3D<float> ,2> > image;
    IMAGE<float>::Read(argv[1],image);

    ARRAYS<VECTOR<VECTOR_3D<float> ,2> > new_image(image);

    int m=image.m, n=image.n;
    for(int i=2;i<m;i++)for(int j=2;j<n;j++)if(image(i,j)==VECTOR_3D<float>()){
        int count=0;VECTOR_3D<float> new_pixel;
        for(int l=i-1;l<=i+1;l++)for(int k=j-1;k<=j+1;k++){
            if(image(l,k)==VECTOR_3D<float>())count++;new_pixel+=image(l,k);}
        new_image(i,j)=(count<=3?new_pixel/float(9.0-count):image(i,j));}

    cout<<"Writing file..."<<endl;
    IMAGE<float>::Write(string("filled_")+argv[1],new_image);
    return 0;
}
