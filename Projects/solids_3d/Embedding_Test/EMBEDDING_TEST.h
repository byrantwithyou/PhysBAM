//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// 1 full featured embedding test
// 2 lots of embedded stuff speed test
// 3 non embedded and simple hard embedded test
// 4 soft embeddings but with regular springs
//#####################################################################
#ifndef __EMBEDDING_TEST__
#define __EMBEDDING_TEST__
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>

namespace PhysBAM{

template<class T_input>
class EMBEDDING_TEST:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;

public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;using BASE::solid_body_collection;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual

    bool target_position;
    bool use_zero_length_springs;
    int constrained_point;
    bool test_3_use_bound;
    T stiffness;
//std::string data_directory;

    EMBEDDING_TEST(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),test_3_use_bound(true),stiffness(1)
    {
        parse_args.Add("-stiffen",&stiffness,"value","stiffness multiplier for various tests");
        parse_args.Parse();
    }

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    //RIGID_BODY_PARTICLES<TV>& rigid_body_particles=deformable_body_collection.rigid_body_particles;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    if(test_number==1){
        TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);
        TRIANGULATED_SURFACE<T>& e_surface=*TRIANGULATED_SURFACE<T>::Create(particles);
        TETRAHEDRALIZED_VOLUME<T>& volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
        SEGMENTED_CURVE<TV>& curve=*SEGMENTED_CURVE<TV>::Create(particles);
        SEGMENTED_CURVE<TV>& curve2=*SEGMENTED_CURVE<TV>::Create(particles);
        VECTOR<int,3> tri_parts(particles.Add_Element(),particles.Add_Element(),particles.Add_Element());
        VECTOR<int,3> etri_parts=tri_parts;
        VECTOR<int,4> tet_parts(particles.Add_Element(),particles.Add_Element(),particles.Add_Element(),particles.Add_Element());
        particles.mass.Subset(tri_parts).Fill((T)1./3);
        particles.mass.Subset(tet_parts).Fill((T)1./4);
        particles.X(tri_parts[0])=TV((T)-0.668156,(T)1.666096,(T)1.202544);
        particles.X(tri_parts[1])=TV((T)0.104949,(T)0.717676,(T)0.022739);
        particles.X(tri_parts[2])=TV((T)0.306855,(T)0.963611,(T)1.181754);
        particles.X(tet_parts[0])=TV((T)-0.433013,(T)-0.018620,(T)-0.750000);
        particles.X(tet_parts[1])=TV((T)-0.433013,(T)-0.018620,(T)0.750000);
        particles.X(tet_parts[2])=TV((T)0.866025,(T)-0.018620,(T)0.000000);
        particles.X(tet_parts[3])=TV((T)0.000000,(T)1.000000,(T)0.000000);
        surface.mesh.elements.Append(tri_parts);//deformable_body_collection.Add_Structure(&surface);
        TETRAHEDRON<T> tet(particles.X.Subset(tet_parts));
        for(int i=0;i<tri_parts.m;i++){
            TV point=particles.X(tri_parts(i));
            if(tet.Inside(point)){
                TV coordinates=tet.First_Three_Barycentric_Coordinates(point);
                binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,tri_parts(i),tet_parts,coordinates));
                etri_parts[i]=particles.Add_Element();
                LOG::cout<<"embedding particle "<<etri_parts[i]<<" to "<<tri_parts[i]<<std::endl;
                particles.X(etri_parts[i])=particles.X(tri_parts[i]);
                curve2.mesh.elements.Append(VECTOR<int,2>(etri_parts[i],tri_parts[i]));
                soft_bindings.Add_Binding(VECTOR<int,2>(etri_parts[i],tri_parts[i]),false);}}
        if(target_position){
            int target_particle=particles.Add_Element();
            particles.X(target_particle)=particles.X(7);
            particles.mass(target_particle)=1; // TODO: make this work for zero mass
            curve.mesh.elements.Append(VECTOR<int,2>(6,7));
            //if(!use_zero_length_springs) curve.mesh.elements.Append(VECTOR<int,2>(1,7));
            deformable_body_collection.Add_Structure(&curve);}
        deformable_body_collection.Add_Structure(&curve2);
        volume.mesh.elements.Append(tet_parts);deformable_body_collection.Add_Structure(&volume);
        e_surface.mesh.elements.Append(etri_parts);
        //tests.Substitute_Soft_Bindings_For_Embedded_Nodes(e_surface,soft_bindings);
        deformable_body_collection.Add_Structure(&e_surface);}
    else if(test_number==2){
        TRIANGULATED_SURFACE<T>& temporary_surface=*TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File<T>(data_directory+"/Rigid_Bodies/sphere.tri",temporary_surface);
        LOG::cout<<"polygon count was "<<temporary_surface.mesh.elements.m<<std::endl;
        //temporary_surface.Loop_Subdivide();
        //temporary_surface.Loop_Subdivide();
        //temporary_surface.Loop_Subdivide();
        LOG::cout<<"polygon count is "<<temporary_surface.mesh.elements.m<<std::endl;
        TRIANGULATED_SURFACE<T>* surface=(TRIANGULATED_SURFACE<T>*)temporary_surface.Append_Particles_And_Create_Copy(particles);
        TETRAHEDRALIZED_VOLUME<T>& temporary_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();surface->Update_Bounding_Box();
        temporary_volume.Initialize_Cube_Mesh_And_Particles(GRID<TV>(VECTOR<int,3>(2,2,2),*surface->bounding_box));
        constrained_point=particles.Size()+1;
        LOG::cout<<"particles.Size() "<<particles.Size()<<std::endl;
        TETRAHEDRALIZED_VOLUME<T>* volume=(TETRAHEDRALIZED_VOLUME<T>*)temporary_volume.Append_Particles_And_Create_Copy(particles);
        
        for(int i=0;i<particles.Size();i++) particles.mass(i)=0;
        //volume->density=float(1)/volume->Total_Volume();
        T density=TV::m==1?1:TV::m==2?100:1000;
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*surface,density);

        deformable_body_collection.Add_Structure(volume);
        //deformable_body_collection.Add_Structure(surface);

        ARRAY<int> tri_particles;
        TRIANGULATED_SURFACE<T>* esurface=(TRIANGULATED_SURFACE<T>::Create(particles));
        HASHTABLE<int,int> hard_to_soft;
        Get_Unique(tri_particles,surface->mesh.elements.Flattened());
        particles.Preallocate(tri_particles.m);
        for(int t=0;t<volume->mesh.elements.m;t++){
            const VECTOR<int,4>& nodes=volume->mesh.elements(t);
            TETRAHEDRON<T> tet(particles.X.Subset(nodes));
            for(int i=0;i<tri_particles.m;i++){int p=tri_particles(i);
                if(!hard_to_soft.Contains(p)  && tet.Inside(particles.X(p))){
                    TV coordinates=tet.First_Three_Barycentric_Coordinates(particles.X(p));
                    binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,p,nodes,coordinates));
                    int soft_particle=particles.Add_Element(); // TODO: make add more efficient
                    particles.mass(soft_particle)=particles.mass(p);
                    particles.X(soft_particle)=particles.X(p);
                    hard_to_soft.Insert(p,soft_particle);
                    soft_bindings.Add_Binding(VECTOR<int,2>(p,soft_particle),!use_zero_length_springs);}}}
        for(int t=0;t<surface->mesh.elements.m;t++){
            VECTOR<int,3> nodes=surface->mesh.elements(t);
            for(int k=0;k<3;k++) nodes[k]=hard_to_soft.Get(nodes[k]);
            esurface->mesh.elements.Append(nodes);}
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*esurface,density);
        deformable_body_collection.Add_Structure(esurface);
    }
    else if(test_number==3){
        TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);
        TETRAHEDRALIZED_VOLUME<T>& volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
        SEGMENTED_CURVE<TV>& curve=*SEGMENTED_CURVE<TV>::Create(particles);
        VECTOR<int,3> tri_parts(particles.Add_Element(),particles.Add_Element(),particles.Add_Element());
        VECTOR<int,4> tet_parts(particles.Add_Element(),particles.Add_Element(),particles.Add_Element(),particles.Add_Element());
        particles.mass.Subset(tri_parts).Fill((T)1./3);
        particles.mass.Subset(tet_parts).Fill((T)1./4);
        particles.X(tri_parts[0])=TV((T)-0.668156,(T)1.666096,(T)1.202544);
        particles.X(tri_parts[1])=TV((T)0.104949,(T)0.717676,(T)0.022739);
        particles.X(tri_parts[2])=TV((T)0.306855,(T)0.963611,(T)1.181754);
        particles.X(tet_parts[0])=TV((T)-0.433013,(T)-0.018620,(T)-0.750000);
        particles.X(tet_parts[1])=TV((T)-0.433013,(T)-0.018620,(T)0.750000);
        particles.X(tet_parts[2])=TV((T)0.866025,(T)-0.018620,(T)0.000000);
        particles.X(tet_parts[3])=TV((T)0.000000,(T)1.000000,(T)0.000000);
        particles.Add_Element(); // fake soft particle
        particles.Add_Element(); // fake constrained point
        particles.X(7)=particles.X(8)=particles.X(3);
        particles.mass(7)=particles.mass(8)=particles.mass(3);
        if(test_3_use_bound){
            particles.Add_Element(); // hard bound particle
            particles.X(9)=particles.X(1);
            particles.mass(9)=particles.mass(1);
            curve.mesh.elements.Append(VECTOR<int,2>(9,7));
            TETRAHEDRON<T> tet(particles.X.Subset(tet_parts));
            solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,10,tet_parts,tet.First_Three_Barycentric_Coordinates(particles.X(9))));
        }
        else curve.mesh.elements.Append(VECTOR<int,2>(3,7));
        curve.mesh.elements.Append(VECTOR<int,2>(7,8));
        deformable_body_collection.Add_Structure(&curve);
        surface.mesh.elements.Append(tri_parts);deformable_body_collection.Add_Structure(&surface);
        volume.mesh.elements.Append(tet_parts);deformable_body_collection.Add_Structure(&volume);
        constrained_point=9;
    }
    if(test_number==4){
        TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);
        TRIANGULATED_SURFACE<T>& e_surface=*TRIANGULATED_SURFACE<T>::Create(particles);
        TETRAHEDRALIZED_VOLUME<T>& volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
        SEGMENTED_CURVE<TV>& curve=*SEGMENTED_CURVE<TV>::Create(particles);
        VECTOR<int,3> tri_parts(particles.Add_Element(),particles.Add_Element(),particles.Add_Element());
        VECTOR<int,3> etri_parts=tri_parts;
        VECTOR<int,4> tet_parts(particles.Add_Element(),particles.Add_Element(),particles.Add_Element(),particles.Add_Element());
        particles.mass.Subset(tri_parts).Fill((T)1./3);
        particles.mass.Subset(tet_parts).Fill((T)1./4);
        particles.X(tri_parts[0])=TV((T)-0.668156,(T)1.666096,(T)1.202544);
        particles.X(tri_parts[1])=TV((T)0.104949,(T)0.717676,(T)0.022739);
        particles.X(tri_parts[2])=TV((T)0.306855,(T)0.963611,(T)1.181754);
        particles.X(tet_parts[0])=TV((T)-0.433013,(T)-0.018620,(T)-0.750000);
        particles.X(tet_parts[1])=TV((T)-0.433013,(T)-0.018620,(T)0.750000);
        particles.X(tet_parts[2])=TV((T)0.866025,(T)-0.018620,(T)0.000000);
        particles.X(tet_parts[3])=TV((T)0.000000,(T)1.000000,(T)0.000000);
        surface.mesh.elements.Append(tri_parts);//deformable_body_collection.Add_Structure(&surface);
        TETRAHEDRON<T> tet(particles.X.Subset(tet_parts));
        for(int i=0;i<tri_parts.m;i++){
            TV point=particles.X(tri_parts(i));
            if(tet.Inside(point)){
                LOG::cout<<"found point "<<tri_parts(i)<<" inside the tet "<<std::endl;
                TV coordinates=tet.First_Three_Barycentric_Coordinates(point);
                binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,tri_parts(i),tet_parts,coordinates));
                etri_parts[i]=particles.Add_Element();assert(etri_parts[i]==8); // 8
                particles.mass(7)=1;
                LOG::cout<<"embedding particle "<<etri_parts[i]<<" to "<<tri_parts[i]<<std::endl;
                particles.X(etri_parts[i])=particles.X(tri_parts[i]);
                soft_bindings.Add_Binding(VECTOR<int,2>(etri_parts[i],tri_parts[i]),false);
            }}
        particles.Add_Element();assert(particles.Size()==9); // 9
        particles.X(8)=particles.X(7);
        particles.mass(8)=1; // TODO: make this work for zero mass
        constrained_point=9;
        curve.mesh.elements.Append(VECTOR<int,2>(7,8));
        curve.mesh.elements.Append(VECTOR<int,2>(1,7));
            //if(!use_zero_length_springs) curve.mesh.elements.Append(VECTOR<int,2>(1,7));
        deformable_body_collection.Add_Structure(&curve);
        volume.mesh.elements.Append(tet_parts);deformable_body_collection.Add_Structure(&volume);
        e_surface.mesh.elements.Append(etri_parts);
        //tests.Substitute_Soft_Bindings_For_Embedded_Nodes(e_surface,soft_bindings);
        deformable_body_collection.Add_Structure(&e_surface);}

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();


    // correct mass
    LOG::printf("before distribute: %P\n",particles.mass);
    binding_list.Distribute_Mass_To_Parents();
    LOG::printf("after distribute before clear: %P\n",particles.mass);
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    LOG::printf("after clear: %P\n",particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
    LOG::printf("soft: %P\n",particles.mass);

    LOG::cout<<"total mass "<<particles.mass.Sum()<<std::endl;

}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    frame_rate=24;
    last_frame=(int)(64*frame_rate);
    solids_parameters.cfl=(T)1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    target_position=true;
    use_zero_length_springs=true;
    output_directory="Embedding_Test/output";

    Get_Initial_Data();
    
    if(test_number==1){
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,tetrahedralized_volume.mesh,0));
        solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)1e6,(T).45,(T).01,(T).25),true,(T).1));
        if(false && target_position){
            SEGMENTED_CURVE<TV>& segmented_curve=deformable_body_collection.template Find_Structure<SEGMENTED_CURVE<TV>&>();
            LINEAR_SPRINGS<TV>* spring_force=Create_Edge_Springs(segmented_curve,100/(1+sqrt((T)2)),(T)1);
            ARRAY<T> restlengths(segmented_curve.mesh.elements.m);restlengths.Fill((T).5);
            //spring_force->Set_Restlength(restlengths);
            spring_force->Clamp_Restlength((T).1);
            spring_force->Set_Overdamping_Fraction(3);
            solid_body_collection.Add_Force(spring_force);}
        soft_bindings.Initialize_Binding_Mesh();
        soft_bindings.binding_mesh->elements.Append(VECTOR<int,2>(7,8));
        soft_bindings.use_impulses_for_collisions.Fill(false);
        if(use_zero_length_springs) solid_body_collection.Add_Force(Create_Edge_Binding_Springs(particles,*soft_bindings.binding_mesh,(T)1e2,(T)1));
        else{
            LINEAR_SPRINGS<TV>* springs=Create_Edge_Springs(particles,*soft_bindings.binding_mesh,(T)1e1,(T)1);
            springs->Clamp_Restlength((T).1);
            springs->Set_Overdamping_Fraction(1);
            solid_body_collection.Add_Force(springs);}}
    else if(test_number==2){
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,tetrahedralized_volume.mesh,0));
        solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)1e6,(T).45,(T).01,(T).25),true,(T).1));
        soft_bindings.Initialize_Binding_Mesh();
        if(use_zero_length_springs){
            soft_bindings.use_impulses_for_collisions.Fill(false);
            LOG::cout<<"stiffness "<<stiffness<<std::endl;
            //solid_body_collection.Add_Force(Create_Edge_Binding_Springs(*soft_bindings.binding_mesh,particles,stiffness,(T)3));}
            LINEAR_SPRINGS<TV>* spring_force=Create_Edge_Springs(particles,*soft_bindings.binding_mesh,stiffness,(T)3);
            solid_body_collection.Add_Force(spring_force);
            spring_force->Clamp_Restlength((T).01);
            spring_force->Set_Overdamping_Fraction(2);}}
    if(test_number==3 || test_number==4){
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,tetrahedralized_volume.mesh,0));
        solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)1e6,(T).45,(T).01,(T).25),true,(T).1));
        SEGMENTED_CURVE<TV>& curve=deformable_body_collection.template Find_Structure<SEGMENTED_CURVE<TV>&>();
        LINEAR_SPRINGS<TV>* spring_force=Create_Edge_Springs(curve,100/(1+sqrt((T)2)),(T)1);
        spring_force->Clamp_Restlength((T).1);
        solid_body_collection.Add_Force(spring_force);}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    if(test_number==1) V(target_position?9:8)=TV();
    //else if(test_number==2) V(constrained_point)=TV();
    else if(test_number==3) V(constrained_point)=TV();
    else if(test_number==4) V(constrained_point)=TV();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    if(test_number==1) V(target_position?9:8)=TV();
    //else if(test_number==2) V(constrained_point)=TV();
    else if(test_number==3) V(constrained_point)=TV();
    else if(test_number==4) V(constrained_point)=TV();
}
//#####################################################################
};
}
#endif
