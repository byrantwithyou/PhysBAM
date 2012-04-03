//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_CUTTING.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_FLUID_SYSTEM.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

using namespace PhysBAM;

typedef float RW;

template<class TV>
void Integration_Test(int argc,char* argv[])
{
    typedef typename TV::SCALAR T;
    static const int d=TV::m;
    typedef VECTOR<int,d> TV_INT;

    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-mu_i",1,"viscosity inside");
    parse_args.Add_Double_Argument("-mu_o",1,"viscosity outside");
    parse_args.Add_Double_Argument("-m",1,"meter scale");
    parse_args.Add_Double_Argument("-s",1,"second scale");
    parse_args.Add_Double_Argument("-kg",1,"kilogram scale");
    parse_args.Add_Integer_Argument("-test",1,"test number");
    parse_args.Add_Integer_Argument("-res",4,"resolution");
    parse_args.Add_Integer_Argument("-cgf",2,"coarse grid factor");
    parse_args.Parse(argc,argv);
    T m=parse_args.Get_Double_Value("-m");
    T s=parse_args.Get_Double_Value("-s");
    T kg=parse_args.Get_Double_Value("-kg");
    VECTOR<T,2> mu;
    mu(0)=parse_args.Get_Double_Value("-mu_i")*kg/(s*(d==3?m:1));
    mu(1)=parse_args.Get_Double_Value("-mu_o")*kg/(s*(d==3?m:1));
    int test=parse_args.Get_Integer_Value("-test");
    int res=parse_args.Get_Integer_Value("-res");
    int cgf=parse_args.Get_Integer_Value("-cgf");

    if(res%cgf){LOG::cerr<<"Resolution must be divisible by coarse grid factor."<<std::endl; exit(-1);}

    TV_INT counts=TV_INT()+res;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1)*m,true);
    GRID<TV> coarse_grid(grid.counts/cgf,grid.domain,true);
    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());
    
    // Setting the domain
    switch(test){
        case 1:{
            SPHERE<TV> sphere(TV()+(T).43*m,(T).31*m);
            for(UNIFORM_GRID_ITERATOR_NODE<TV> it(coarse_grid);it.Valid();it.Next())
                phi(it.index)=sphere.Signed_Distance(it.Location());
            
            break;}
        default:{
            LOG::cerr<<"Unknown test number."<<std::endl;
            exit(-1);
            break;}}

    INTERFACE_FLUID_SYSTEM<TV> ifs(grid,coarse_grid,phi);
    ifs.Set_Matrix(mu);

    ARRAY<TV,TV_INT> f_body[2];
    ARRAY<TV> f_interface;
    f_interface.Resize(ifs.object.mesh.elements.m);
    for(int s=0;s<2;s++) f_body[s].Resize(grid.Domain_Indices());

    // Setting the forces
    switch(test){
        case 1:{
            
            break;}
        default:{
            LOG::cerr<<"Unknown test number."<<std::endl;
            exit(-1);
            break;}}

    ifs.Set_RHS(f_body,f_interface);

    std::string output_directory="output";
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    printf("\n");
    for(int i=0;i<d;i++) printf("%c [%i %i) ", "uvw"[i], ifs.index_range_u[i].min_corner, ifs.index_range_u[i].max_corner);
    printf("p [%i %i) ", ifs.index_range_p.min_corner, ifs.index_range_p.max_corner);
    for(int i=0;i<d;i++) printf("%cq [%i %i) ", "uvw"[i], ifs.index_range_q[i].min_corner, ifs.index_range_q[i].max_corner);
    printf("\n");
}

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    Integration_Test<VECTOR<double,2> >(argc,argv);

    return 0;
}
