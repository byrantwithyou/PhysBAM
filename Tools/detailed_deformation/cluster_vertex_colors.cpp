//#####################################################################
// Copyright 2005, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <cstring>
#include <fstream>

using namespace PhysBAM;
using namespace std;

//#################################################################
// Function LuvApproxDistance from www.compupphase.com/cmetric.html
//#################################################################
inline float LuvApproxDistance(VECTOR_3D<float> color_one,VECTOR_3D<float> color_two) {
    float r_avg=color_one.x+color_two.x/2.0;
    float delta_r=color_one.x-color_two.x;float delta_g=color_one.y-color_two.y;float delta_b=color_one.z-color_two.z;
    return sqrt((2+r_avg)*delta_r*delta_r+4*delta_g*delta_g+(2+(1-r_avg))*delta_b*delta_b);
}

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{

    if (argc!=3){cout<<"Usage: cluster_vertex_colors input_colors num_clusters"<<endl;exit(-1);}

    int num_clusters=atoi(argv[2]);

    ARRAY<VECTOR_3D<float> > vertex_colors;
    cout<<"Reading color files..."<<endl;
    FILE_UTILITIES::Read_From_File<float>(argv[1],vertex_colors);

    ARRAY<int> vertex_clusters(vertex_colors.m);

    if(num_clusters!=3){cout<<"For now, restrict to three clusters for skin, lips, and hair."<<endl;exit(-1);}
    
    // Clusters has three elements: the number of assigned vertices, the sum of all those vertex colors, and the cluster center.
    ARRAY<TRIPLE<int,VECTOR_3D<float>,VECTOR_3D<float> > > clusters(num_clusters);
    clusters(1).z=VECTOR_3D<float>(0.95,0.8,0.8);clusters(2).z=VECTOR_3D<float>(0.83,0.65,0.60);clusters(3).z=VECTOR_3D<float>(0.48,0.32,0.25);

    int iterations=0;bool converged;
    
    do{ // Cluster vertex colors
        cout<<iterations<<endl;iterations++;
        converged=true;
        for(int i=1;i<=vertex_colors.m;i++){ // Assign vertices to new clusters
            float min_distance=LuvApproxDistance(vertex_colors(i),clusters(1).z);int cluster=1;
            for(int j=2;j<=num_clusters;j++){
                float temp_distance=LuvApproxDistance(vertex_colors(i),clusters(j).z);
                if(min_distance>temp_distance){min_distance=temp_distance;cluster=j;}}
            if(vertex_clusters(i)!=cluster){converged=false;vertex_clusters(i)=cluster;}
            clusters(cluster).x++;clusters(cluster).y+=vertex_colors(i);}
        for(int i=1;i<=num_clusters;i++){ // Update cluster centers
            clusters(i).z=clusters(i).y/clusters(i).x;
            cout<<clusters(i).x<<" "<<clusters(i).z<<endl;
            clusters(i).x=0;clusters(i).y=VECTOR_3D<float>();}
    } while(!converged);    

    // Update vertex colors for temporary viewing
    for(int i=1;i<vertex_colors.m;i++) vertex_colors(i)=clusters(vertex_clusters(i)).z;
        
    FILE_UTILITIES::Write_To_File<float>("vertex_cluster_colors.col",vertex_colors);

    return 0;
} 

