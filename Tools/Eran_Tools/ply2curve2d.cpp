#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <fstream>
#include <iostream>
#include "../../Personal_Libraries/Eran_Library/POLYGONAL_CURVE.h"
#include <Particles/SOLIDS_PARTICLES.h>
#include <string.h> // for gcc 2.96

using namespace PhysBAM;

bool closed=false,verbose=false;
double max_segment_length=0;

template<class T> void Convert(const char *input_filename, const char *output_filename, PARSE_ARGS& parse_args)
{
    bool closed = false;
    T max_segment_length = 0;

    if(closed) std::cout << "Closed" << std::endl;
    if(max_segment_length) std::cout << "Max segment length " << max_segment_length << std::endl;

    std::ifstream input(input_filename);
    if (!input)
    {
        std::cerr << "Input file " << input_filename << " does not exist!" << std::endl;
        exit(-1);
    }

    ERAN_LIB::POLYGONAL_CURVE polygonal_curve;
    polygonal_curve.Read(input);
    if(verbose) polygonal_curve.Print();

    SEGMENTED_CURVE_2D<T> *segmented_curve = SEGMENTED_CURVE_2D<T>::Create();
    SEGMENT_MESH &segment_mesh = segmented_curve->segment_mesh;
    SOLIDS_PARTICLES<T,VECTOR<T,2> > &particles = segmented_curve->particles;

    particles.Increase_Array_Size(polygonal_curve.vertices.m); // approximate since we may supersample
    for (int i = 1; i <= polygonal_curve.vertices.m; i++)
    {
        int idx = particles.Add_Particle();
        particles.X(idx) = polygonal_curve.vertices(i);
        if(verbose) std::cout << "Vertex " << i << ": " << particles.X(idx) << std::endl;
    }

    segment_mesh.Clean_Memory();
    for (int i = 1; i <= polygonal_curve.curves.m; i++) 
    {
        for (int j = 1; j <= polygonal_curve.curves(i).m; j++) {
            if(j==polygonal_curve.curves(i).m && !closed) break;
            int node1=polygonal_curve.curves(i)(j),node2=(j==polygonal_curve.curves(i).m)?polygonal_curve.curves(i)(1):polygonal_curve.curves(i)(j+1);
            T length=(particles.X(node2)-particles.X(node1)).Magnitude();
            int subdiv=1;
            if(max_segment_length){subdiv=max(1,(int)ceil(length/max_segment_length));}
            std::cout << "subdiv " << subdiv << std::endl;
            int last_node=node1;
            for(int s=0;s<subdiv;s++){int next_node;
                if(s==subdiv){next_node=node2;}
                else{next_node=particles.Add_Particle();
                    particles.X(next_node)=particles.X(node1)+(T)s/subdiv*(particles.X(node2)-particles.X(node1));}
                segment_mesh.segments.Append(last_node,next_node);
                if(verbose) std::cout << "Segment " << segment_mesh.segments.m << ": " << last_node << " " << next_node << std::endl;
                last_node=next_node;}
        }
    }
    segment_mesh.number_nodes=particles.number;

    std::ofstream output(output_filename, std::ios::binary);
    segmented_curve->template Write<T>(output);
}

int main(int argc, char *argv[])
{
    bool type_double = false;
    
    std::string input_filename;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-closed",&closed,"closed");
    parse_args.Add("-v",&verbose,"verbose");
    parse_args.Add("-max_segment_length",&max_segment_length,"max_segment_length");
    parse_args.Extra(&input_filename,"ply file", "ply file to convert");
    parse_args.Parse();

    if (!FILE_UTILITIES::Is_Ply2D_File(input_filename))
    {
        std::cerr << "Not a ply2d file: " << input_filename << std::endl;
        return -1;
    }

    char basename[256], output_filename[256];
    strcpy(basename, FILE_UTILITIES::Get_Basename(input_filename).c_str());
    sprintf(output_filename, "%s.curve2d", basename);

    std::cout << "Input filename: " << input_filename << std::endl;
    std::cout << "Output filename: " << output_filename << " [" << (type_double?"double":"float") << "]" << std::endl;

    if(!type_double) Convert<float>(input_filename, output_filename, parse_args);
    else{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Convert<double>(input_filename, output_filename, parse_args);
#else
        std::cerr << "Compiled without double support" << std::endl;
#endif
    }
}
