//#####################################################################
// Copyright 2005, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <cstring>
#include <fstream>
#include "../../Projects/face/FACE_DRIVER.h"
#include "../../Public_Library/Deformable_Objects/DEFORMABLE_OBJECT_LIST_3D.h"

#define lower_cutoff .12

using namespace PhysBAM;
using namespace std;

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    if (argc!=8){cout<<"Usage: create_mesh_from_offsets_list list_directory frame condensation_file full_tetrahedral_volume triangle_mesh offset_file output_dir"<<::endl;exit(-1);}

    DEFORMABLE_OBJECT_LIST_3D<float> deformable_object_list;
    ARRAY<int> condensation_map;    
    TETRAHEDRON_MESH tet_mesh;
    TRIANGLE_MESH tri_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tet_particles;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > tri_particles;
    TETRAHEDRALIZED_VOLUME<float> tet_vol(tet_mesh,tet_particles);
    TRIANGULATED_SURFACE<float> tri_surf(tri_mesh,tri_particles);

    cout<<"Reading single-access files..."<<endl;string prefix(argv[1]);
    deformable_object_list.Read_Static_Variables<float>(prefix+"/");
    FILE_UTILITIES::Read_From_File<int>(argv[3],condensation_map);
    FILE_UTILITIES::Read_From_File<float>(argv[4],tet_vol);
    FILE_UTILITIES::Read_From_File<float>(argv[5],tri_surf);

    ARRAY<PAIR<int,VECTOR_3D<float> > > offset_particles(tri_particles.number);
    FILE_UTILITIES::Read_From_File<float>(argv[6],offset_particles);

    cout<<"Initializing tetrahderal triangulated surface..."<<endl;
    tet_vol.Initialize_Triangulated_Surface();
    cout<<"Updating tetrahedral triangle list..."<<endl;
    tet_vol.triangulated_surface->Update_Triangle_List();
    cout<<"Updating tetrahedral tetrahedron list..."<<endl;
    tet_vol.Update_Tetrahedron_List();

    ARRAY<int> particles_to_update;
    for(int i=1;i<=tet_particles.number;i++)if(condensation_map(i)!=0)particles_to_update.Append(i);

    int frame=atoi(argv[2]);

    char frame_c[5];sprintf(frame_c,"%d",frame);
    string frame_s(frame_c);
    cout<<"Processing frame "+frame_s+"..."<<endl;
    cout<<"Reading dynamic variables "<<endl;
    deformable_object_list.Read_Dynamic_Variables<float>(prefix+"/",frame);
        
    cout<<"Updating condensation based positions..."<<endl;
    // Use the condensation map to update tetrahedral particles appropriately
    for(int i=1;i<=particles_to_update.m;i++)tet_particles.X(particles_to_update(i))=deformable_object_list(1).tetrahedralized_volume->particles.X(condensation_map(particles_to_update(i)));

    cout<<"Initializing deformable object triangulated surface..."<<endl;
    deformable_object_list(1).tetrahedralized_volume->Initialize_Triangulated_Surface();
    cout<<"Initializing deformable object triangle list..."<<endl;
    deformable_object_list(1).tetrahedralized_volume->triangulated_surface->Update_Triangle_List();
    cout<<"Updating tetrahedral tetrahedron list..."<<endl;
    tet_vol.Update_Tetrahedron_List();
    
    cout<<"Transforming hires mesh via barycentric coords..."<<endl;
    // Reconstruct triangle mesh points from offsets        
    for(int p=1;p<=tri_particles.number;p++){if(tri_particles.X(p).y>-lower_cutoff)
        tri_particles.X(p)=(*tet_vol.tetrahedron_list)(offset_particles(p).x).Point_From_Barycentric_Coordinates(offset_particles(p).y);}
    
    /*
    FACE_CONTROL_PARAMETERS<float> face_control_parameters;
    ACTIVATION_CONTROL_SET<float>* act_control_set;
    ATTACHMENT_FRAME_CONTROL_SET<float>* att_control_set;

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(string(argv[1])+"/face_control_set_types");
    if(!input){cout<<"Face control set types not found."<<endl;exit(-1);}
    int control_sets;Read_Binary<float>(*input,control_sets);cout<<"Procesing "<<control_sets<<" control sets..."<<endl;
    for(int i=1;i<=control_sets;i++){
        int control_set_type;Read_Binary<float>(*input,control_set_type);
        if(control_set_type==FACE_CONTROL_SET<float>::ACTIVATION){
            act_control_set=new ACTIVATION_CONTROL_SET<float>;ARRAY<std::string> muscle_names;
            FILE_UTILITIES::Read_From_File<float>(STRING_UTILITIES::string_sprintf("%s/face_control_set_%d.muscle_names",argv[1],i),muscle_names);
            for(int j=1;j<=muscle_names.m;j++) act_control_set->Add_Activation(muscle_names(j));
            face_control_parameters.list.Append(act_control_set);}
        else if(control_set_type==FACE_CONTROL_SET<float>::ATTACHMENT_FRAME){
            ARRAY<ARRAY<int> > nodes_dummy;
            att_control_set=new ATTACHMENT_FRAME_CONTROL_SET<float>(deformable_object_list(1).particles.X.array,nodes_dummy,0);
            att_control_set->Read_Jaw_Joint_From_File<float>(STRING_UTILITIES::string_sprintf("%s/face_control_set_%d.jaw_joint_parameters",argv[1],i));
            face_control_parameters.list.Append(att_control_set);}}
    delete input;

    if(!att_control_set){cout<<"No rotations..."<<endl;exit(-1);}

    cout<<"Applying rotations..."<<endl;
    // Apply affine transformations (rotation) to flesh
    for(int p=1;p<=tri_particles.number;p++)
        tri_particles.X(p)=att_control_set->cranium_transform.affine_transform*tri_particles.X(p)+att_control_set->cranium_transform.translation;
    */   
    cout<<"Writing triangle mesh..."<<endl;
    FILE_UTILITIES::Write_To_File<float>(string(argv[7])+"/deformed_coord_mesh_"+frame_s+".tri",tri_surf);

    return 0;
}
//#################################################################
