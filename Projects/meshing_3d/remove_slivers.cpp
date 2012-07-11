//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Program remove_slivers
//##################################################################### 
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
using namespace PhysBAM;
//#######################################################################
// Function main
//#######################################################################
int main(int argc,char **argv)
{
    typedef float T;
    typedef float RW;

    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-aspect_ratio",6);
    parse_args.Set_Extra_Arguments(2,"<input.tet> <output.tet>");
    parse_args.Parse();

    std::string input_file=parse_args.Extra_Arg(0);
    std::string output_file=parse_args.Extra_Arg(1);
    T aspect_ratio=parse_args.Get_Double_Value("-aspect_ratio"); 

    TETRAHEDRON_MESH tetrahedron_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> > particles;
    TETRAHEDRALIZED_VOLUME<T> tetrahedralized_volume(tetrahedron_mesh,particles);

    FILE_UTILITIES::Read_From_File<RW>(input_file,tetrahedralized_volume);
    tetrahedralized_volume.Update_Tetrahedron_List();

    LOG::cout<<"tetrahedra = "<<tetrahedron_mesh.tetrahedrons.m<<std::endl;
    for(int t=0;t<tetrahedron_mesh.tetrahedrons.m;t++){
        TETRAHEDRON<T> tet=(*tetrahedralized_volume.tetrahedron_list)(t);
        if(tet.Aspect_Ratio() > aspect_ratio){
            LOG::cout<<"removing tet "<<t<<" with aspect ratio "<<tet.Aspect_Ratio()<<std::endl;
            tetrahedron_mesh.tetrahedrons(1,t)=0;}}
    tetrahedron_mesh.Delete_Tetrahedrons_With_Missing_Nodes();
    tetrahedralized_volume.Discard_Valence_Zero_Particles_And_Renumber();

    FILE_UTILITIES::Write_To_File<RW>(output_file,tetrahedralized_volume);

    return 0;
}
//#######################################################################
