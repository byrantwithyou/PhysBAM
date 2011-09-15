//#####################################################################
// Copyright 2006, Kevin Der, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURGICAL_EXAMPLE
//#####################################################################
#ifndef __SURGICAL_EXAMPLE__
#define __SURGICAL_EXAMPLE__

#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT.h"
#include "BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING.h"
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Deformable_Objects/DEFORMABLE_OBJECT_EVOLUTION_BACKWARD_EULER.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
#include <Solids_And_Fluids/SOLIDS_PARAMETERS_3D.h>
#include <string.h>

namespace PhysBAM{

template<class T,class RW>
class SURGICAL_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;

    SOLIDS_STANDARD_TESTS<TV,RW> tests;
    ARRAY<int> attached_nodes;
    ARRAY<VECTOR<T,3> > attached_node_locations,attached_node_velocities;
    SEGMENT_MESH hook_attachment_mesh;
    TETRAHEDRALIZED_VOLUME<T>* tissue_tet_vol;
    PARTICLES<T,VECTOR<T,3> > material_space_particles;
    TETRAHEDRALIZED_VOLUME<T>* material_space_tissue_tet_vol;
    TRIANGULATED_SURFACE<T>* tissue_surface;
    LINEAR_SPRINGS<T,TV>* ls;
    BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>* btsc;
    ARRAY<int> constrained_tets;
    ARRAY<VECTOR<T,3> > constrained_tet_location_weights;
    ARRAY<VECTOR<T,3> > constrained_tet_node_locations;
    ARRAY<VECTOR<int,2> > hook_start_and_end_frame;
    ARRAY<ARRAY<int> > hook_move_frames;
    ARRAY<VECTOR<T,3> > hook_material_positions;
    ARRAY<ARRAY<VECTOR<T,3> > > hook_moves;
    T constraint_stiffness,constraint_damping;
    ARRAY<T> hook_stiffnesses;
    BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>* btzls;
    ARRAY<VECTOR<int,2> > suture_constrained_tets;
    ARRAY<VECTOR<VECTOR<T,3>,2> > suture_constrained_tet_weights;
    ARRAY<VECTOR<int,2> > suture_start_and_end_frame;
    ARRAY<VECTOR<VECTOR<T,3>,2> > suture_material_positions;
    T suture_stiffness,suture_damping;
    ARRAY<T> suture_stiffnesses;
    bool run_from_history_file,waiting_to_initialize;
    std::string history_filename;
    int frame_delay_before_initializing_from_history;

    SURGICAL_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),tests(*this,solids_parameters),ls(0),frame_delay_before_initializing_from_history(30)
    {
        output_directory=STRING_UTILITIES::string_sprintf("Surgical_Example/output");
        last_frame=240;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        //solids_parameters.cfl*=10;
        //solids_parameters.Set_Evolution(new DEFORMABLE_OBJECT_EVOLUTION_BACKWARD_EULER<TV>(solids_parameters));
        solids_parameters.perform_self_collision=false;
        solids_parameters.perform_collision_body_collisions=false;
        verbose=false;
        verbose_dt=false;
        constraint_stiffness=(T)500;constraint_damping=(T)1;
        suture_stiffness=(T)2500;suture_damping=(T)0;
        //solids_parameters.min_dt=(T)1/30;
        run_from_history_file=true;
        waiting_to_initialize=true;
        history_filename="Surgical_Example/savedSurgery.hst";
        restart=false;
        restart_frame=20;
    }

    // Unused callbacks
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
    void Initialize_Bodies() PHYSBAM_OVERRIDE
    {
        Initialize_Bodies_With_Embedding();
        if(run_from_history_file)
            Parse_History_File(history_filename,false);
    }
//#####################################################################
// Function Add_Hook
//#####################################################################
    int Add_Hook(const TV hook_location)
    {
        DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
        PARTICLES<T,TV>& particles=deformable_object.particles;

        TV hook_location_changed=hook_location;
        hook_location_changed.y-=(T).01;

        VECTOR<T,3> weights;int closest_tet=0;
        for(int t=1;t<=tissue_tet_vol->mesh.elements.m;t++){
            int i,j,k,l;tissue_tet_vol->mesh.elements.Get(t,i,j,k,l);
            weights=TETRAHEDRON<T>::Barycentric_Coordinates(hook_location_changed,particles.X(i),particles.X(j),particles.X(k),particles.X(l));
            if(weights.x > -(T)1e-4 && weights.y > -(T)1e-4 && weights.z > -(T)1e-4 && (T)1-(weights.x+weights.y+weights.z) > -(T)1e-4){
                closest_tet=t;
                break;}}

        assert(closest_tet);
        constrained_tets.Append(closest_tet);
        constrained_tet_node_locations.Append(hook_location_changed);
        constrained_tet_location_weights.Append(weights);
        btsc->youngs_modulus.Append(constraint_stiffness);
        btsc->damping.Append(constraint_damping*constraint_stiffness);

        btsc->Initialize_Weights_Outer_Product();
        return constrained_tets.m;

    }
//#####################################################################
// Function Delete_Hook
//#####################################################################
    void Delete_Hook(int hook_number)
    {
        btsc->youngs_modulus(hook_number)=(T)0;
        btsc->damping(hook_number)=(T)0;
    }
//#####################################################################
// Function Reactivate_Hook
//#####################################################################
    void Reactivate_Hook(int hook_number)
    {
        btsc->youngs_modulus(hook_number)=constraint_stiffness;
        btsc->damping(hook_number)=constraint_damping;
    }
//#####################################################################
// Function Add_Suture
//#####################################################################
    int Add_Suture(const TV suture_end1,const TV suture_end2)
    {
        DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
        PARTICLES<T,TV>& particles=deformable_object.particles;

        TV suture_end1_changed=suture_end1;TV suture_end2_changed=suture_end2;
        suture_end1_changed.y-=(T).01;suture_end2_changed.y-=(T).01;

        VECTOR<T,3> weights1;int closest_tet1=0;
        for(int t=1;t<=tissue_tet_vol->mesh.elements.m;t++){
            int i,j,k,l;tissue_tet_vol->mesh.elements.Get(t,i,j,k,l);
            weights1=TETRAHEDRON<T>::Barycentric_Coordinates(suture_end1_changed,particles.X(i),particles.X(j),particles.X(k),particles.X(l));
            if(weights1.x > -(T)1e-4 && weights1.y > -(T)1e-4 && weights1.z > -(T)1e-4 && (T)1-(weights1.x+weights1.y+weights1.z) > -(T)1e-4){
                closest_tet1=t;
                break;}}

        VECTOR<T,3> weights2;int closest_tet2=0;
        for(int t=1;t<=tissue_tet_vol->mesh.elements.m;t++){
            int i,j,k,l;tissue_tet_vol->mesh.elements.Get(t,i,j,k,l);
            weights2=TETRAHEDRON<T>::Barycentric_Coordinates(suture_end2_changed,particles.X(i),particles.X(j),particles.X(k),particles.X(l));
            if(weights2.x > -(T)1e-4 && weights2.y > -(T)1e-4 && weights2.z > -(T)1e-4 && (T)1-(weights2.x+weights2.y+weights2.z) > -(T)1e-4){
                closest_tet2=t;
                break;}}

        assert(closest_tet1);
        assert(closest_tet2);
 
        suture_constrained_tets.Resize(suture_constrained_tets.m+1);
 
        suture_constrained_tets(1,suture_constrained_tets.m)=closest_tet1;
        suture_constrained_tets(2,suture_constrained_tets.m)=closest_tet2;
        suture_constrained_tet_weights.Resize(suture_constrained_tet_weights.m+1);
        suture_constrained_tet_weights(1,suture_constrained_tet_weights.m)=weights1;
        suture_constrained_tet_weights(2,suture_constrained_tet_weights.m)=weights2;
        btzls->youngs_modulus.Append(suture_stiffness);
        btzls->damping.Append(suture_damping*suture_stiffness);

        btzls->Initialize_Weights_Outer_Product();
        return constrained_tets.m;
    }
//#####################################################################
// Function Delete_Suture
//#####################################################################
    void Delete_Suture(int suture_number)
    {
        btzls->youngs_modulus(suture_number)=(T)0;
        btzls->damping(suture_number)=(T)0;
    }
//#####################################################################
// Function Reactivate_Suture
//#####################################################################
    void Reactivate_Suture(int suture_number)
    {
        btzls->youngs_modulus(suture_number)=suture_stiffness;
        btzls->damping(suture_number)=suture_damping;
    }
//#####################################################################
// Function Initialize_Bodies_With_Embedding
//#####################################################################
    void Initialize_Bodies_With_Embedding()
    {
        DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
        PARTICLES<T,TV>& particles=deformable_object.particles;

        TETRAHEDRALIZED_VOLUME<T>& tet_volume=tests.Create_Tetrahedralized_Volume("../../Personal_Libraries/Joey_Library/data/Output_Dup/embedding_volume.tet",RIGID_BODY_STATE<TV>(FRAME_3D<T>(TV(0,0,0))),false,true);
        tissue_tet_vol=&tet_volume;   
        tet_volume.Update_Bounding_Box();
        tet_volume.Set_Mass_Of_Particles(true);
        //for(int p=1;p<=particles.array_collection->Size();p++)
        // material_space_particles.X(material_space_particles.array_collection->Add_Element())=particles.X(p);
 
        //material_space_tissue_tet_vol=new TETRAHEDRALIZED_VOLUME<T>(tet_volume.mesh,material_space_particles);
        TRIANGULATED_SURFACE<T>& boundary_surface=*TRIANGULATED_SURFACE<T>::Create(tet_volume.particles);
        TRIANGLE_MESH boundary_mesh;
        tissue_surface=&boundary_surface;
        FILE_UTILITIES::Read_From_File<RW>("../../Personal_Libraries/Joey_Library/data/Output_Dup/boundary_mesh.tri",boundary_mesh);
        boundary_surface.mesh.Initialize_Mesh(boundary_mesh);
        boundary_surface.Set_Mass_Of_Particles(true);
        boundary_surface.Update_Bounding_Box();
        solid_body_collection.deformable_object.Add_Structure(&boundary_surface);
        //for(int p=1;p<=particles.array_collection->Size();p++) if(particles.X(p).y>tet_volume.bounding_box->ymax-(T).01) attached_nodes.Append(p);

        ARRAY<ARRAY<int> > embedding_map;
        FILE_UTILITIES::Read_From_File<RW>("../../Personal_Libraries/Joey_Library/data/Output_Dup/embedding_map",embedding_map);
        for(int p=1;p<=embedding_map.m;p++){
            ARRAY<int>& parents=embedding_map(p);
            if(parents.m==2){
                VECTOR<T,2> weights=SEGMENT_3D<T>::Barycentric_Coordinates(particles.X(p),particles.X(parents(1)),particles.X(parents(2)));
                if(weights.Min()<0) LOG::cout<<"Negative barycentric coordinates on particle "<<p<<" : "<<weights<<", parents : "<<parents<<std::endl;
                for(int i=1;i<=weights.m;i++) weights(i)=max((T)0,min(weights(i),(T)1));
                solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<T,TV,2>(particles,p,VECTOR<int,2>(parents(1),parents(2)),weights));}
            if(parents.m==3){
                VECTOR<T,3> weights=TRIANGLE_3D<T>::Barycentric_Coordinates(particles.X(p),particles.X(parents(1)),particles.X(parents(2)),particles.X(parents(3)));
                if(weights.Min()<0) LOG::cout<<"Negative barycentric coordinates on particle "<<p<<" : "<<weights<<", parents : "<<parents<<std::endl;
                for(int i=1;i<=weights.m;i++) weights(i)=max((T)0,min(weights(i),(T)1));
                solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<T,TV,3>(particles,p,VECTOR<int,3>(parents(1),parents(2),parents(3)),weights));}
            if(parents.m==4){
                VECTOR<T,3> weights=TETRAHEDRON<T>::Barycentric_Coordinates(particles.X(p),particles.X(parents(1)),particles.X(parents(2)),particles.X(parents(3)),particles.X(parents(4)));
                if(weights.Min()<0) LOG::cout<<"Negative barycentric coordinates on particle "<<p<<" : "<<weights<<", parents : "<<parents<<std::endl;
                for(int i=1;i<=weights.m;i++) weights(i)=max((T)0,min(weights(i),(T)1));
                solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<T,TV,4>(particles,p,VECTOR<int,4>(parents(1),parents(2),parents(3),parents(4)),weights));}}
        ////tests.Substitute_Soft_Bindings_For_Embedded_Nodes(boundary_surface,solid_body_collection.deformable_body_collection.soft_bindings);
        //tests.Add_Ground();
        //tests.Add_Rigid_Body("sphere",2,0.3);

        //// correct number nodes
        ////for(int i=1;i<=deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

        // setup binding springs
        ////ARRAY<bool>::copy(false,solid_body_collection.deformable_body_collection.soft_bindings.use_impulses_for_collisions);
        ////solid_body_collection.deformable_body_collection.soft_bindings.Initialize_Binding_Mesh();

        // add structures and rigid bodies to collisions
        deformable_object.collisions.collision_structures.Append(&boundary_surface);
        ////solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(&boundary_surface);
        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

        // correct number nodes
        for(int i=1;i<=deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

        // correct mass
        solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents(particles.mass.array);
        solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass.array);
        particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
        solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

        //solid_body_collection.Add_Force(new GRAVITY<T,TV>(tet_volume));
        //solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tet_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));
        solid_body_collection.Add_Force(Create_Quasistatic_Diagonalized_Finite_Volume(tet_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)1e4,(T).4,(T).1,(T).25),false,false));
        //ls=Create_Edge_Springs(hook_attachment_mesh,particles,(T)2000,(T)1,true,(T).1,true,(T)0,false);
        btsc= new BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>(particles,constrained_tets,constrained_tet_location_weights,constrained_tet_node_locations,tet_volume.mesh,constraint_stiffness,constraint_damping);
        solid_body_collection.Add_Force(btsc);

        btzls=new BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>(particles,suture_constrained_tets,suture_constrained_tet_weights,tet_volume.mesh,suture_stiffness,suture_damping);
        solid_body_collection.Add_Force(btzls);

        solid_body_collection.Update_Fragments();

        SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Initialize_Bodies();
    }
//#####################################################################
// Function Set_Hook_Position
//#####################################################################
    void Set_Hook_Position(const TV hook_position,const int hook_index)
    {
        constrained_tet_node_locations(hook_index)=hook_position;
        //attached_node_locations(hook_index)=hook_position;
    }
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
    {
        if(frame>frame_delay_before_initializing_from_history && run_from_history_file && waiting_to_initialize)
            Initialize_From_History_After_Delay();
        else if(run_from_history_file){
            for(int h=1;h<=hook_start_and_end_frame.m;h++){
                if(frame==hook_start_and_end_frame(1,h))
                    Reactivate_Hook(h);
                else if(frame==hook_start_and_end_frame(2,h))
                    Delete_Hook(h);}
            for(int s=1;s<=suture_start_and_end_frame.m;s++){
                if(frame==suture_start_and_end_frame(1,s))
                    Reactivate_Suture(s);
                else if(frame==suture_start_and_end_frame(2,s))
                    Delete_Suture(s);}}
    }
//#####################################################################
// Function Initialize_From_History_After_Delay
//#####################################################################
    void Initialize_From_History_After_Delay()
    {
        //add hooks
        for(int h=1;h<=hook_material_positions.m;h++){
            Add_Hook(hook_material_positions(h));
            Delete_Hook(h);//just setting the stiffness to zero until active
        }

        //add sutures
        for(int s=1;s<=suture_material_positions.m;s++){
            Add_Suture(suture_material_positions(1,s),suture_material_positions(2,s));
            Delete_Suture(s);//just setting the stiffness to zero until active
        }

        waiting_to_initialize=false;
    }
//#####################################################################
// Function Parse_History_File
//#####################################################################
    void Parse_History_File(std::string filename,const bool verbose=false)
    {
        char line[256];
        std::string command;
        char *pch;
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename,false);
 
        input->getline(line,256);
        while(!input->eof()){

            pch=0;pch=strtok(line," ");
            if(line[0])
                command=std::string(pch);

            if(command == "addHook"){
                pch=strtok(NULL," ");
                int hook_num=(int)atof(pch)+1;
                float x,y,z;
                pch=strtok(NULL," ");
                x=(float)atof(pch);
                pch=strtok(NULL," ");
                y=(float)atof(pch);
                pch=strtok(NULL," ");
                z=(float)atof(pch);
                int frame_num;
                pch=strtok(NULL," ");
                frame_num=(int)atof(pch);
                if(verbose)
                    LOG::cout<<"addHook "<<hook_num<<" at (x,y,z) = "<<"("<<x<<","<<y<<","<<z<<") at frame = "<<frame_num <<"\n"<<std::endl;
                //Add_Hook(VECTOR<float,3>(x,y,z));
                hook_material_positions.Append(VECTOR<float,3>(x,y,z));
                hook_start_and_end_frame.Resize(hook_start_and_end_frame.m+1);
                hook_start_and_end_frame(1,hook_start_and_end_frame.m)=frame_num;
                hook_start_and_end_frame(2,hook_start_and_end_frame.m)=-1;
                //Delete_Hook(hook_num);//not really deleting, only setting stiffness to 0 until active
            }
            else if(command == "moveHook"){
                pch=strtok(NULL," ");
                int hook_num=(int)atof(pch)+1;
                float x,y,z;
                pch=strtok(NULL," ");
                x=(float)atof(pch);
                pch=strtok(NULL," ");
                y=(float)atof(pch);
                pch=strtok(NULL," ");
                z=(float)atof(pch);
                int frame_num;
                pch=strtok(NULL," ");
                frame_num=(int)atof(pch);
                hook_move_frames(hook_num).Append(frame_num);
                hook_moves(hook_num).Append(VECTOR<float,3>(x,y,z));
                if(verbose)
                    LOG::cout<<"moveHook "<<hook_num<<" to (x,y,z) = "<<"("<<x<<","<<y<<","<<z<<") at frame = "<<frame_num <<"\n"<<std::endl;
            }
            else if(command == "deleteHook"){
                pch=strtok(NULL," ");
                int hook_num=(int)atof(pch)+1;
                int frame_num;
                pch=strtok(NULL," ");
                frame_num=(int)atof(pch);
                hook_start_and_end_frame(2,hook_num)=-1;
                if(verbose)
                    LOG::cout<<"deleteHook number "<< hook_num<<" at frame number "<<frame_num<<"\n"<<std::endl;
            }
            else if(command == "addSuture"){
                pch=strtok(NULL," ");
                int suture_num=(int)atof(pch)+1;
                float x,y,z;
                pch=strtok(NULL," ");
                x=(float)atof(pch);
                pch=strtok(NULL," ");
                y=(float)atof(pch);
                pch=strtok(NULL," ");
                z=(float)atof(pch);
                VECTOR<float,3> suture_location1=VECTOR<float,3>(x,y,z);
                pch=strtok(NULL," ");
                x=(float)atof(pch);
                pch=strtok(NULL," ");
                y=(float)atof(pch);
                pch=strtok(NULL," ");
                z=(float)atof(pch);
                VECTOR<float,3> suture_location2=VECTOR<float,3>(x,y,z);
                int frame_num;
                pch=strtok(NULL," ");
                frame_num=(int)atof(pch);
                suture_start_and_end_frame.Resize(suture_start_and_end_frame.m+1);
                suture_start_and_end_frame(1,suture_start_and_end_frame.m)=frame_num;
                suture_start_and_end_frame(2,suture_start_and_end_frame.m)=-1;
                suture_material_positions.Resize(suture_material_positions.m+1);
                suture_material_positions(1,suture_material_positions.m)=suture_location1;
                suture_material_positions(2,suture_material_positions.m)=suture_location2;
                //Add_Suture(suture_location1,suture_location2);
                //Delete_Suture(suture_start_and_end_frame.m);//only setting the stiffness to 0 until active
                if(verbose){
                    LOG::cout<<"addSuture num "<<suture_num<<" sp1 = ("<<suture_location1.x<<","<<suture_location1.y<<","<<suture_location1.z <<") ";
                    LOG::cout<<"sp2 = ("<<suture_location2.x<<","<<suture_location2.y<<","<<suture_location2.z<<" at frame number "<<frame_num <<"\n"<<std::endl;}}
            else if(command == "deleteSuture"){
                pch=strtok(NULL," ");
                int suture_num=(int)atof(pch)+1;
                int frame_num;
                pch=strtok(NULL," ");
                frame_num=(int)atof(pch);
                suture_start_and_end_frame(2,suture_num)=frame_num;
                if(verbose)
                    LOG::cout<<"deleteSuture number "<< suture_num<<" at frame number "<<frame_num <<"\n"<<std::endl;}

            input->getline(line,256);}
    }
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
    {
        TV v1(1,0,0),v2(0,1,0),v3(0,0,1);
        if(fragment_id==1) 
            V.Subset(solid_body_collection.deformable_object.dynamic_particles_of_fragment(fragment_id)).Fill(v1);
        if(fragment_id==2) 
            V.Subset(solid_body_collection.deformable_object.dynamic_particles_of_fragment(fragment_id)).Fill(v2);
        if(fragment_id==4) 
            V.Subset(solid_body_collection.deformable_object.dynamic_particles_of_fragment(fragment_id)).Fill(v3);
    }
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
    {
        int unconstrained_fragment_id=3;
        if(fragment_id!=unconstrained_fragment_id) 
            V.Subset(solid_body_collection.deformable_object.dynamic_particles_of_fragment(fragment_id)).Fill(TV());
    }
//#####################################################################
};
}
#endif
