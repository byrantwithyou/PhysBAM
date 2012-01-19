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

    if (argc!=2){cout<<"Usage: col2img vertex_colors"<<::endl;exit(-1);}

    ARRAY<VECTOR_3D<float> > colors;
    cout<<"Reading colors..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[1],colors);
 
    int m,n;for(m=(int)sqrt(colors.m);m>=1;m--)for(n=1;n<=colors.m;n++)if(m*n==colors.m)goto found;
    cout<<"Could not be converted..."<<endl;exit(-1);

  found:
    cout<<m<<" "<<n<<endl;
    ARRAYS<VECTOR<VECTOR_3D<float> ,2> > image(1,m,1,n);
    for(int i=0;i<m;i++)for(int j=0;j<n;j++) image(i,j)=colors((j-1)*m+i);
    
    cout<<"Writing file..."<<endl;
    IMAGE<float>::Write("vertex_colors.bmp",image);
    return 0;
}
