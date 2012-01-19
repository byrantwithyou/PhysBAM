#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
using namespace PhysBAM;

template<class T> void Remove_Degenerate_Triangles(const char* input_filename,const char* output_filename,T threshold)
{
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create_From_File(input_filename,true);
    ARRAY<int> nondegenerate_triangle_indices;
    ARRAY<int>& triangles = surface->triangle_mesh.triangles;
    int number_of_triangles=triangles.m;
    const ARRAYS<VECTOR<VECTOR_3D<T> ,1> >& X = surface->particles.X;
    for(int t=0;t<triangles.m;t++){
        if(TRIANGLE_3D<T>::Area(X(triangles(1,t)),X(triangles(2,t)),X(triangles(3,t)))>threshold)nondegenerate_triangle_indices.Append(t);}
    ARRAY<int> new_triangles(3,1,nondegenerate_triangle_indices.m);
    for(int t=0;t<new_triangles.m;t++)for(int k=0;k<3;k++)
        new_triangles(k,t)=triangles(k,nondegenerate_triangle_indices(t));
    ARRAY<int>::exchange_arrays(triangles,new_triangles);
    surface->triangle_mesh.number_nodes=ARRAY<int>::max(triangles);
    printf("%d total before, %d discarded\n",number_of_triangles,number_of_triangles-nondegenerate_triangle_indices.m);
    std::fstream output;output.open(output_filename,std::ios::out|std::ios::binary);
    if(output.is_open()){surface->Write(output);output.close();}
}

int main(int argc,char* argv[])
{
    PARSE_ARGS args;
    args.Add_String_Argument("-type","float","-type","float or double");
    args.Add_String_Argument("-input","<input>","-input","Input filename");
    args.Add_String_Argument("-output","<output>","-output","Output filename");
    args.Add_Double_Argument("-threshold",0,"-threshold","Area size of and below that will be discarded");
    args.Parse(argc,argv);
    const char* type=args.Get_String_Value("-type");
    const char* input_filename=args.Get_String_Value("-input");
    const char* output_filename=args.Get_String_Value("-output");
    double threshold=args.Get_Double_Value("-threshold");
    if(strcmp(type,"double")==0){Remove_Degenerate_Triangles<double>(input_filename,output_filename,threshold);}
    else if(strcmp(type,"float")==0){Remove_Degenerate_Triangles<float>(input_filename,output_filename,(float)threshold);}
    else {fprintf(stderr,"Invalid type specified.\n");exit(1);}
}