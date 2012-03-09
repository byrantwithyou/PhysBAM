//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN  
//##################################################################### 
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <Deformable_Objects/DEFORMABLE_OBJECT_LIST_3D.h>
#include <Fracture/TRIANGLES_OF_MATERIAL_3D.h>
#include <Geometry/EMBEDDED_TRIANGULATED_SURFACE.h>
#include <Level_Sets/LEVELSET_TRIANGULATED_OBJECT.h>
using namespace PhysBAM;

// for i in `seq 0 \`cat last_frame\``; do ../render_preprocess $i; done

template<class T,class RW> void Process(int argc,char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Set_Extra_Arguments(1, "<frame>");
    parse_args.Parse(argc,argv);

    int frame=-10;
    if(parse_args.Num_Extra_Args() >= 1) frame=atoi(parse_args.Extra_Arg(0).c_str());
    else{LOG::cout<<"Incorrect.\n";exit(1);}

    LOG::cout<<"Preprocessing frame "<<frame<<std::endl;

    std::string prefix="";
    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);
    std::string f_minus=STRING_UTILITIES::string_sprintf(".%d",frame-1);
    std::string o=prefix+STRING_UTILITIES::string_sprintf("melting_%d_",1);

    DEFORMABLE_OBJECT_LIST_3D<T> deformable_object_list;
    deformable_object_list.template Read_Static_Variables<RW>(prefix,frame);
    deformable_object_list.template Read_Dynamic_Variables<RW>(prefix,frame);
    assert(deformable_object_list.deformable_objects.m==1);
    DEFORMABLE_OBJECT_3D<T>& deformable_object=deformable_object_list(1);
    assert(deformable_object.embedded_triangulated_surface);
    deformable_object.triangles_of_material->Create_Material_Surface();
    DEFORMABLE_PARTICLES<T,VECTOR<T,3> >& particles=deformable_object.particles;
    DEFORMABLE_PARTICLES<T,VECTOR<T,3> >& material_particles=deformable_object.triangles_of_material->material_surface.particles;
    
    LEVELSET_TRIANGULATED_OBJECT<T,VECTOR<T,3> > levelset(*deformable_object.embedded_triangulated_surface);
    FILE_UTILITIES::Read_From_File<RW>(o+"levelset"+f,levelset);
    ARRAY<ARRAY<T>*> temperature,reaction;
    FILE_UTILITIES::Read_From_File<RW>(prefix+"melting_temperature"+f,temperature);
    FILE_UTILITIES::Read_From_File<RW>(prefix+"melting_reaction"+f,reaction);
    ARRAY<VECTOR_2D<T> >& node_locations=levelset.grid.Node_Locations();
    BOX_2D<T> box=levelset.grid.uniform_grid.Domain();

    ARRAY<T> maximum_temperature;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(prefix+"rendering_melting_grid_maximum_temperature"+f_minus,true,false);
    if(input){Read_Binary<RW>(*input,maximum_temperature);delete input;}
    else maximum_temperature.Resize(levelset.grid.number_of_nodes);

    ARRAY<VECTOR_2D<T> > texture_coordinates(material_particles.array_collection->Size());
    ARRAY<T> material_temperature(material_particles.array_collection->Size()),material_reaction(material_particles.array_collection->Size()),material_maximum_temperature(material_particles.array_collection->Size());

    for(int n=0;n<levelset.grid.number_of_nodes;n++)if(levelset.node_to_particle_mapping(n)) particles.V(levelset.node_to_particle_mapping(n))=VECTOR<T,3>(node_locations(n));
    deformable_object.triangles_of_material->Update_Particle_Velocities();
    for(int p=0;p<material_particles.array_collection->Size();p++)texture_coordinates(p)=(VECTOR_2D<T>(material_particles.V(p).x,material_particles.V(p).y)-box.Minimum_Corner())/box.Size();

    for(int n=0;n<levelset.grid.number_of_nodes;n++)if(levelset.node_to_particle_mapping(n)) particles.V(levelset.node_to_particle_mapping(n)).x=(*temperature(1))(n);
    deformable_object.triangles_of_material->Update_Particle_Velocities();
    for(int p=0;p<material_particles.array_collection->Size();p++)material_temperature(p)=material_particles.V(p).x;
    
    for(int n=0;n<levelset.grid.number_of_nodes;n++)if(levelset.node_to_particle_mapping(n)) particles.V(levelset.node_to_particle_mapping(n)).x=(*reaction(1))(n);
    deformable_object.triangles_of_material->Update_Particle_Velocities();
    for(int p=0;p<material_particles.array_collection->Size();p++)material_reaction(p)=material_particles.V(p).x;
    
    for(int n=0;n<levelset.grid.number_of_nodes;n++)maximum_temperature(n)=max(maximum_temperature(n),(*temperature(1))(n));
    for(int n=0;n<levelset.grid.number_of_nodes;n++)if(levelset.node_to_particle_mapping(n)) particles.V(levelset.node_to_particle_mapping(n)).x=maximum_temperature(n);
    deformable_object.triangles_of_material->Update_Particle_Velocities();
    for(int p=0;p<material_particles.array_collection->Size();p++)material_maximum_temperature(p)=material_particles.V(p).x;
    
    FILE_UTILITIES::Write_To_File<RW>(prefix+"rendering_texture_coordinates"+f,texture_coordinates);
    FILE_UTILITIES::Write_To_File<RW>(prefix+"rendering_melting_temperature"+f,material_temperature);
    FILE_UTILITIES::Write_To_File<RW>(prefix+"rendering_melting_reaction"+f,material_reaction);
    FILE_UTILITIES::Write_To_File<RW>(prefix+"rendering_melting_grid_maximum_temperature"+f,maximum_temperature);
    FILE_UTILITIES::Write_To_File<RW>(prefix+"rendering_melting_maximum_temperature"+f,material_maximum_temperature);
}

int main(int argc,char* argv[])
{
    Process<float,float>(argc,argv);
    return 0;
}
//#####################################################################

