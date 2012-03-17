//#####################################################################
// Copyright 2006-2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINDING_SPRINGS_TEST
//#####################################################################
#ifndef __BINDING_SPRINGS_TEST__
#define __BINDING_SPRINGS_TEST__

#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_MATERIAL_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class BINDING_SPRINGS_TEST:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::parse_args;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions;using BASE::solid_body_collection; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;

    SEGMENT_MESH spring_segment_mesh;
    DEFORMABLES_FORCES<TV>* spring_force;
    T stiffness;
    T overdamping_fraction;

    BINDING_SPRINGS_TEST(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),stiffness(1100),overdamping_fraction(0)
    {
    }

    virtual ~BINDING_SPRINGS_TEST()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
    parse_args->Add_Double_Argument("-binding_stiffness",100,"");
    parse_args->Add_Double_Argument("-binding_overdamping_fraction",1,"");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options()
{
    BASE::Parse_Options();
    if(parse_args->Is_Value_Set("-binding_stiffness")) stiffness=(T)parse_args->Get_Double_Value("-binding_stiffness");
    if(parse_args->Is_Value_Set("-binding_overdamping_fraction")) overdamping_fraction=(T)parse_args->Get_Double_Value("-binding_overdamping_fraction");
} 
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    switch(test_number){
      case 0: Stability_Test();break;
      case 1: Segments();break;
      case 2: Embedded_Segments();break;
      case 3: Sphere();break;
      case 4: Embedded_Sphere();break;
      case 5: Falling_Sphere();break;
      case 6: Falling_Embedded_Sphere();break;
      default: LOG::cerr<<"Unknown test number "<<test_number<<std::endl;exit(1);}

    // correct number nodes
    spring_segment_mesh.Set_Number_Nodes(particles.Size());
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    output_directory=STRING_UTILITIES::string_sprintf("Binding_Springs_Test/Test_%d",test_number);
    frame_rate=24;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    last_frame=200;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
    solids_parameters.cfl=1;

    Get_Initial_Data();

//    spring_force=Create_Zero_Length_Linear_Springs(stiffness,overdamping_fraction);
    if(solid_body_collection.deformable_body_collection.soft_bindings.bindings.m) spring_force=Create_Edge_Binding_Springs(particles,spring_segment_mesh,stiffness,overdamping_fraction);
    else spring_force=Create_Edge_Zero_Length_Springs(particles,spring_segment_mesh,stiffness,overdamping_fraction);

    solid_body_collection.Add_Force(spring_force);

    solid_body_collection.Update_Simulated_Particles();

    // ------TEST
/*    LOG::cout<<"Effective masses="<<std::endl;
    for(int p=0;p<particles.Size();p++) LOG::cout<<"\teffective mass("<<p<<")="<<particles.mass.effective_mass(p)<<std::endl;

    // try a different effective mass
    ARRAY<TV> F(particles.Size());
    for(int i=0;i<solid_body_collection.deformable_body_collection.binding_list.bindings.m;i++) F(solid_body_collection.deformable_body_collection.binding_list.bindings(i)->particle_index)=TV::All_Ones_Vector();
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(1,F);

    for(int b=0;b<solid_body_collection.deformable_body_collection.soft_bindings.bindings.m;b++){VECTOR<int,2>& binding=solid_body_collection.deformable_body_collection.soft_bindings.bindings(b);
        particles.mass.effective_mass(binding.x)=particles.mass.array(binding.x)=particles.mass.effective_mass(binding.y)/solid_body_collection.deformable_body_collection.binding_list.Binding(binding.y)->Embedded_Force(F).x;
        particles.mass.one_over_effective_mass(binding.x)=particles.mass.one_over_mass(binding.x)=(T)1/particles.mass.array(binding.x);}

    LOG::cout<<"Effective masses="<<std::endl;
    for(int p=0;p<particles.Size();p++) LOG::cout<<"\teffective mass("<<p<<")="<<particles.mass.effective_mass(p)<<std::endl;

    // recompute damping
    ((IMPLICIT_ZERO_LENGTH_SPRINGS<T,TV>*)spring_force)->Set_Overdamping_Fraction(overdamping_fraction);
*/    // ------END TEST
}
//#####################################################################
// Function Create_Zero_Length_Linear_Springs
//#####################################################################
LINEAR_SPRINGS<TV>* Create_Zero_Length_Linear_Springs(const T youngs_modulus,const T overdamping_fraction)
{
    LINEAR_SPRINGS<TV>* spring_force=Create_Edge_Springs(solid_body_collection.deformable_body_collection.particles,solid_body_collection.rigid_body_collection,spring_segment_mesh,(T)2e1,
        overdamping_fraction);
    ARRAY<T>::copy(0,spring_force->visual_restlength);spring_force->Clamp_Restlength(1);
    spring_force->Set_Overdamping_Fraction(overdamping_fraction);
    spring_force->use_implicit_velocity_independent_forces=true;
    return spring_force;
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(LINEAR_SPRINGS<TV>* linear_springs=solid_body_collection.template Find_Force<LINEAR_SPRINGS<TV>*>())
        linear_springs->Print_Deformation_Statistics();

/*    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    LOG::cout<<"Frame "<<frame<<std::endl;
      for(int p=0;p<particles.Size();p++) LOG::cout<<"X("<<p<<")="<<particles.X(p)<<std::endl;*/
}
//#####################################################################
// Function Stability_Test
//#####################################################################
void Stability_Test()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    particles.Add_Elements(2);
    particles.mass.Fill((T)1);
    particles.X(1)=TV(-2,0,0);
    particles.X(2)=TV(2,0,0);

    SEGMENTED_CURVE<TV>& segmented_curve1=*SEGMENTED_CURVE<TV>::Create(particles);
    segmented_curve1.mesh.elements.Append(VECTOR<int,2>(1,2));
    spring_segment_mesh.elements.Append(VECTOR<int,2>(1,2));

    deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve1);
}
//#####################################################################
// Function Segments
//#####################################################################
void Segments()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    particles.Add_Elements(6);
    particles.mass.Fill((T)1);
    particles.X(1)=TV(-3,2,0);
    particles.X(2)=TV(-2,1,0);
    particles.X(3)=TV(-1,0,0);
    particles.X(4)=TV(3,2,0);
    particles.X(5)=TV(2,1,0);
    particles.X(6)=TV(1,0,0);

    SEGMENTED_CURVE<TV>& segmented_curve1=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& segmented_curve2=*SEGMENTED_CURVE<TV>::Create(particles);
    segmented_curve1.mesh.elements.Append(VECTOR<int,2>(1,2));
    segmented_curve1.mesh.elements.Append(VECTOR<int,2>(2,3));
    segmented_curve2.mesh.elements.Append(VECTOR<int,2>(4,5));
    segmented_curve2.mesh.elements.Append(VECTOR<int,2>(5,6));
    spring_segment_mesh.elements.Append(VECTOR<int,2>(1,4));
    spring_segment_mesh.elements.Append(VECTOR<int,2>(2,5));
    spring_segment_mesh.elements.Append(VECTOR<int,2>(3,6));

    deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve1);
    deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve2);
}
//#####################################################################
// Function Sphere
//#####################################################################
void Sphere()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    TRIANGULATED_SURFACE<T>& sphere1=tests.Create_Triangulated_Object(data_directory+"/Triangulated_Surfaces/sphere.tri",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-1,(T)3,0))),true,true);
    int sphere_particles=particles.Size();
    TRIANGULATED_SURFACE<T>& sphere2=tests.Create_Triangulated_Object(data_directory+"/Triangulated_Surfaces/sphere.tri",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(1,(T)3,0))),true,true);
    sphere1.Update_Number_Nodes();
    sphere2.Update_Number_Nodes();

    for(int i=0;i<sphere_particles;i++) spring_segment_mesh.elements.Append(VECTOR<int,2>(i,i+sphere_particles));

    deformable_body_collection.deformable_geometry.Add_Structure(&sphere1);
    deformable_body_collection.deformable_geometry.Add_Structure(&sphere2);
}
//#####################################################################
// Function Falling_Sphere
//#####################################################################
void Falling_Sphere()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",
        RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,1000);
    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));

    // duplicate surface
    tetrahedralized_volume.Update_Number_Nodes();tetrahedralized_volume.Initialize_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& drifted_surface=tests.Create_Drifted_Surface(*tetrahedralized_volume.triangulated_surface,solid_body_collection.deformable_body_collection.soft_bindings);
    deformable_body_collection.deformable_geometry.Add_Structure(&drifted_surface);
    spring_segment_mesh.elements.Append_Elements(solid_body_collection.deformable_body_collection.soft_bindings.bindings);

    tests.Add_Ground();

    // collisions
    deformable_body_collection.collisions.collision_structures.Append(&drifted_surface);
}
//#####################################################################
// Function Embedded_Sphere
//#####################################################################
void Embedded_Sphere()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    // create object and soft bindings
    EMBEDDED_MATERIAL_SURFACE<TV,3>& embedding=tests.Create_Embedded_Tetrahedralized_Volume(SPHERE<TV>(TV(),(T).9),RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true);
    embedding.Update_Binding_List_From_Embedding(solid_body_collection.deformable_body_collection);
//    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(embedding.material_surface,solid_body_collection.deformable_body_collection.binding_list,solid_body_collection.deformable_body_collection.soft_bindings);
    
    // duplicate surface
    TRIANGULATED_SURFACE<T>& drifted_surface=tests.Create_Drifted_Surface(embedding.material_surface,solid_body_collection.deformable_body_collection.soft_bindings);
    deformable_body_collection.deformable_geometry.Add_Structure(&drifted_surface);
    drifted_surface.particles.X.Subset(solid_body_collection.deformable_body_collection.soft_bindings.bindings.Project(1)).Project(1)+=1;
    spring_segment_mesh.elements.Append_Elements(solid_body_collection.deformable_body_collection.soft_bindings.bindings);

    embedding.Update_Number_Nodes();

    spring_segment_mesh.elements.Append_Elements(solid_body_collection.deformable_body_collection.soft_bindings.bindings);
}
//#####################################################################
// Function Falling_Embedded_Sphere
//#####################################################################
void Falling_Embedded_Sphere()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    // create object and soft bindings
    EMBEDDED_MATERIAL_SURFACE<TV,3>& embedding=tests.Create_Embedded_Tetrahedralized_Volume(SPHERE<TV>(TV(),(T).9),RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true);
    embedding.Update_Binding_List_From_Embedding(solid_body_collection.deformable_body_collection);
    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(embedding.material_surface,solid_body_collection.deformable_body_collection.soft_bindings);
    embedding.Update_Number_Nodes();

    // add forces
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedding.embedded_object.simplicial_object;
    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,tetrahedralized_volume.mesh,0));
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),
        true,(T).1));

    spring_segment_mesh.elements.Append_Elements(solid_body_collection.deformable_body_collection.soft_bindings.bindings);

    tests.Add_Ground();

    // collisions
    deformable_body_collection.collisions.collision_structures.Append(&embedding.material_surface);
}
//#####################################################################
// Function Embedded_Segments
//#####################################################################
void Embedded_Segments()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    particles.Add_Elements(7);

    // segment particles
    particles.X(1)=TV(-2,2,0);
    particles.X(2)=TV(0,0,0);
    particles.X(3)=TV(-2,-2,0);
    particles.mass(1)=1;
    particles.mass(2)=2;
    particles.mass(3)=1;

    // embedded particles
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,4,VECTOR<int,2>(1,2),VECTOR<T,2>(.25,.5)));
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,5,VECTOR<int,2>(2,3),VECTOR<T,2>(.5,.5)));
    for(int i=0;i<solid_body_collection.deformable_body_collection.binding_list.bindings.m;i++) solid_body_collection.deformable_body_collection.binding_list.bindings(i)->Clamp_To_Embedded_Position();
    particles.mass(4)=0;
    particles.mass(5)=0;

    // drifted particles
    TV offset(1,2,0);
    particles.X(6)=particles.X(4)+offset;
    particles.X(7)=particles.X(5)+offset;
    particles.mass(6)=1;
    particles.mass(7)=1;
    
    SEGMENTED_CURVE<TV>& segmented_curve1=*SEGMENTED_CURVE<TV>::Create(particles);
    segmented_curve1.mesh.elements.Append(VECTOR<int,2>(1,2));
    segmented_curve1.mesh.elements.Append(VECTOR<int,2>(2,3));
    deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve1);

    solid_body_collection.deformable_body_collection.soft_bindings.Add_Binding(VECTOR<int,2>(6,4),false);
    solid_body_collection.deformable_body_collection.soft_bindings.Add_Binding(VECTOR<int,2>(7,5),false);
    solid_body_collection.deformable_body_collection.soft_bindings.Initialize_Binding_Mesh();
    deformable_body_collection.deformable_geometry.Add_Structure(new SEGMENTED_CURVE<TV>(*solid_body_collection.deformable_body_collection.soft_bindings.binding_mesh,particles));
    
    spring_segment_mesh.elements.Append_Elements(solid_body_collection.deformable_body_collection.soft_bindings.bindings);
}
//#####################################################################
};
}
#endif
