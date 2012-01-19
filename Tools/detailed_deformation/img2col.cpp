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

    if (argc!=2){cout<<"Usage: img2col image"<<::endl;exit(-1);}

    ARRAYS<VECTOR<VECTOR_3D<float> ,2> > image;
    IMAGE<float>::Read(argv[1],image);

    int m=image.m, n=image.n;
    ARRAY<VECTOR_3D<float> > colors(m*n);
    for(int i=0;i<m;i++)for(int j=0;j<n;j++) colors((j-1)*m+i)=image(i,j);
    
    cout<<"Writing file..."<<endl;
    FILE_UTILITIES::Write_To_File<float>("image_colors.col",colors);
    return 0;
}
