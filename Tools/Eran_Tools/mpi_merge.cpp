#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

PARSE_ARGS parse_args;
bool verbose = false;
std::string master_grid_filename, domain_filename_pattern, array_filename_pattern, output_filename;
int number_of_nodes;

template<class T,class T2> void Merge_Array_1D()
{
    std::cout << "1d not supported yet" << std::endl;
}

template<class T,class T2> void Merge_Array_2D()
{
    GRID_2D<T> master_grid;
    if (verbose) std::cout << "Reading master grid: " << master_grid_filename << std::endl;
    std::ifstream master_grid_input(master_grid_filename.c_str(),std::ios::binary);
    master_grid.template Read<T>(master_grid_input);

    if (verbose) std::cout << "Got master grid: m=" << master_grid.m << ", n=" << master_grid.n << ", domain=" << master_grid.Domain() << std::endl;

    ARRAYS<VECTOR<T2,2> > array;

    for (int node=0;node<number_of_nodes;node++)
    {
        char tmp[512];
        sprintf(tmp, domain_filename_pattern.c_str(), node);
        if (verbose) std::cout << "Reading node domain: " << tmp << std::endl;
        std::ifstream domain_input(tmp,std::ios::binary);
        BOX_2D<int> domain;
        domain.template Read<T>(domain_input);
        if (verbose) std::cout << "Got domain: " << domain << std::endl;

        sprintf(tmp, array_filename_pattern.c_str(), node);
        if (verbose) std::cout << "Reading node array: " << tmp << std::endl;
        std::ifstream array_input(tmp,std::ios::binary);
        ARRAYS<VECTOR<T2,2> > node_array;
        node_array.template Read<T>(array_input);

        if (node_array.m != (domain.xmax-domain.xmin+1) || node_array.n != (domain.ymax-domain.ymin+1)) {
            std::cout << "Array has size " << node_array.m << " x " << node_array.n << std::endl;
            std::cout << "But was expecting domain " << domain << std::endl;
            exit(-1);
        }

        if (node==0) {
            // Initialize size of array here since we need to know length and we get that after reading first array
            array.Resize(node_array.length, 1, master_grid.m, 1, master_grid.n);
            if (verbose) std::cout << "Setting array to length " << node_array.length << ", m=" << master_grid.m << ", n=" << master_grid.n << std::endl;
        }

        for (int i=domain.xmin;i<=domain.xmax;i++) for (int j=domain.ymin;j<=domain.ymax;j++) for (int k=1;k<=array.length;k++)
        {
            array(k,i,j)=node_array(k,i-domain.xmin+1,j-domain.ymin+1);
        }
    }

    std::ofstream output_stream(output_filename.c_str(),std::ios::binary);
    array.template Write<T>(output_stream);
}

template<class T,class T2> void Merge_Array_3D()
{
    std::cout << "3d not supported yet" << std::endl;
}

template<class T,class T2> void By_Composite_Type()
{
    if (parse_args.Get_Option_Value("-3d")) Merge_Array_3D<T,T2>();
    else if (parse_args.Get_Option_Value("-2d")) Merge_Array_2D<T,T2>();
    else Merge_Array_1D<T,T2>();
}

template<class T> void By_Base_Type()
{
    if (parse_args.Get_Option_Value("-vec3d")) By_Composite_Type<T, VECTOR<T,3> >();
    else if (parse_args.Get_Option_Value("-vec2d")) By_Composite_Type<T, VECTOR<T,2> >();
    else By_Composite_Type<T,T>();
}

int main(int argc, char *argv[])
{
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-int");
    parse_args.Add_Option_Argument("-1d");
    parse_args.Add_Option_Argument("-2d");
    parse_args.Add_Option_Argument("-3d");
    parse_args.Add_Option_Argument("-vec2d");
    parse_args.Add_Option_Argument("-vec3d");
    parse_args.Add_Option_Argument("-vec3d");
    parse_args.Add_Option_Argument("-v");
    parse_args.Add_Integer_Argument("-np", 1, "<num nodes>", "number of nodes");
    parse_args.Add_String_Argument("-o", "merged.out");
    parse_args.Set_Extra_Arguments(3, "<master grid> <domain filename pattern> <array filename pattern>");
    parse_args.Parse(argc,argv);

    if (parse_args.Num_Extra_Args() != 3) {
        std::cerr << "Missing arguments" << std::endl;
        return 1;
    }

    verbose = parse_args.Get_Option_Value("-v");

    master_grid_filename = parse_args.Extra_Arg(0);
    domain_filename_pattern = parse_args.Extra_Arg(1);
    array_filename_pattern = parse_args.Extra_Arg(2);
    output_filename = parse_args.Get_String_Value("-o");

    number_of_nodes = parse_args.Get_Integer_Value("-np");

    if (parse_args.Get_Option_Value("-double")) {
        if (verbose) std::cout << "Base type = double" << std::endl;
        By_Base_Type<double>();
    }
    else if (parse_args.Get_Option_Value("-float")) {
        if (verbose) std::cout << "Base type = float" << std::endl;
        By_Base_Type<float>();
    } else {
        if (verbose) std::cout << "Base type = int" << std::endl;
        By_Base_Type<int>();
    }
}
