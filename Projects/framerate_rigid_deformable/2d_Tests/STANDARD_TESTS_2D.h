//#####################################################################
// Copyright 2010, Jonathan Su, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_2D
//#####################################################################
// 1. Contrained spring
// 2. Single spring
// 3. Two springs symmetric
// 4. Two springs asymmetric
// 5. Single spring moving to the right
// 6. Spring mesh
// 7. Contrained spring with gravity
// 8. Two springs asymmetric constrained
// 9. Cloth panel hanging under gravity
// 10. Three springs, extended triangle
// 11. Three springs, compressed triangle
// 12. Ten springs
//#####################################################################
#ifndef __STANDARD_TESTS_2D__
#define __STANDARD_TESTS_2D__

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <climits>

namespace PhysBAM{

template<class T>
struct ERROR_COLOR_MAP
{
    T mn,mx;
    bool use_log,reverse_order,use_two_sets;
    INTERPOLATION_CURVE<T,VECTOR<T,3> >& colors;

    ERROR_COLOR_MAP(T min_value,T max_value,bool log_scale,bool reverse,bool two_sets)
        :mn(min_value),mx(max_value),use_log(log_scale),reverse_order(reverse),use_two_sets(two_sets),colors(*new INTERPOLATION_CURVE<T,VECTOR<T,3> >)
    {
        if(use_log)
        {
            mn=std::log(std::abs(mn));
            mx=std::log(std::abs(mx));
        }
        PHYSBAM_ASSERT(mn<mx);
        T a=mn,d=(mx-mn)/(6*(1+use_two_sets));
        if(reverse_order){a=mx,d=-d;}
        colors.Add_Control_Point(a+0*d,VECTOR<T,3>(1,1,1));
        colors.Add_Control_Point(a+1*d,VECTOR<T,3>(1,0,0));
        colors.Add_Control_Point(a+2*d,VECTOR<T,3>(1,1,0));
        colors.Add_Control_Point(a+3*d,VECTOR<T,3>(0,1,0));
        colors.Add_Control_Point(a+4*d,VECTOR<T,3>(0,1,1));
        colors.Add_Control_Point(a+5*d,VECTOR<T,3>(0,0,1));
        if(use_two_sets){
            colors.Add_Control_Point(a+6*d,VECTOR<T,3>(.3,0,.3));
            colors.Add_Control_Point(a+7*d,VECTOR<T,3>(.3,0,0));
            colors.Add_Control_Point(a+8*d,VECTOR<T,3>(.3,.3,0));
            colors.Add_Control_Point(a+9*d,VECTOR<T,3>(0,.3,0));
            colors.Add_Control_Point(a+10*d,VECTOR<T,3>(0,.3,.3));
            colors.Add_Control_Point(a+11*d,VECTOR<T,3>(0,0,.3));}
        colors.Add_Control_Point(a+6*(1+use_two_sets)*d,VECTOR<T,3>());
    }

    ~ERROR_COLOR_MAP()
    {
        delete &colors;
    }

    VECTOR<T,3> operator()(T x) const
    {
        if(use_log){if(x==0) x=mn;else x=std::log(std::abs(x));}
        return colors.Value(x);
    }
};

template<class T_input>
class STANDARD_TESTS_2D:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef typename FORCE_ELEMENTS::ITERATOR SEGMENT_ITERATOR;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    T stiffness_multiplier;
    T damping_multiplier;
    bool use_be,use_tr;
    GEOMETRY_PARTICLES<TV> residual_energy_particles;

    int number_side_panels;
    T aspect_ratio;
    T side_length;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;using BASE::frame_rate;

    STANDARD_TESTS_2D(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),use_be(true),use_tr(false),
        number_side_panels(40),aspect_ratio((T)1.7),side_length((T)1.0)
    {
    }

    ~STANDARD_TESTS_2D()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Fragments() PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Double_Argument("-stiffen",1,"","stiffness multiplier for various tests");
    parse_args->Add_Double_Argument("-dampen",1,"","damping multiplier for various tests");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Option_Argument("-use_be","use backward Euler");
    parse_args->Add_Option_Argument("-use_tr","use trapezoid rule");
    parse_args->Add_Integer_Argument("-side_panels",INT_MAX,"Cloth side panels");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("2d_Tests/Test_%d",test_number);
    stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen");
    damping_multiplier=(T)parse_args->Get_Double_Value("-dampen");
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    if(parse_args->Is_Value_Set("-use_tr")){use_tr=true;use_be=false;}
    if(parse_args->Is_Value_Set("-use_be")){use_tr=false;use_be=true;}
    if(parse_args->Is_Value_Set("-side_panels")){
        number_side_panels=parse_args->Get_Integer_Value("-side_panels");
        output_directory+=STRING_UTILITIES::string_sprintf("_sidepanels=%d",number_side_panels);}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    switch(test_number){
        case 1: Constrained_Spring_Test();break;
        case 2: Single_Spring_Test();break;
        case 3: Two_Spring_Test((T)0);break;
        case 4: Two_Spring_Test((T).1);break;
        case 5: Single_Spring_Moving_To_Right_Test();break;
        case 6: Spring_Mesh();break;
        case 7: Constrained_Spring_Gravity_Test();break;
        case 8: Constrained_Two_Spring_Test((T).1);break;
        case 9: Hanging_Cloth_Test();break;
        case 10: Three_Spring_Test((T)2);break;
        case 11: Three_Spring_Test((T)0.1);break;
        case 12: Ten_Spring_Test();break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    residual_energy_particles.array_collection->template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);

    if(use_tr) solids_parameters.use_trapezoidal_rule_for_velocities=true;
    else if(use_be) solids_parameters.use_trapezoidal_rule_for_velocities=false;
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(soft_bindings);

    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // add forces
    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 10:
        case 11:
        case 12:
            {Create_Spring_Forces((T).95);break;}
        case 8:
            {Create_Spring_Forces((T)1);break;}
        case 7:
            {Create_Spring_Forces((T)1);
            solid_body_collection.Add_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,true));
            break;}
        case 9:
            {Create_Spring_Forces();
            solid_body_collection.Add_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,true));
            break;}
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
}
//#####################################################################
// Function Constrained_Spring_Test
//#####################################################################
void Constrained_Spring_Test()
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=FLT_MAX;
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=(T)1;
    segmented_curve.particles.X(new_edge_node1)=TV((T)0,(T)0);segmented_curve.particles.X(new_edge_node2)=TV((T)1,(T)0);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    last_frame=1000;
}
//#####################################################################
// Function Constrained_Spring_Gravity_Test
//#####################################################################
void Constrained_Spring_Gravity_Test()
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=FLT_MAX;
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=(T)1;
    segmented_curve.particles.X(new_edge_node1)=TV((T)0,(T)0);segmented_curve.particles.X(new_edge_node2)=TV((T)0,(T)-1);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    last_frame=1000;
}
//#####################################################################
// Function Single_Spring_Test
//#####################################################################
void Single_Spring_Test()
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=(T)1;
    segmented_curve.particles.X(new_edge_node1)=TV((T)0,(T)0);segmented_curve.particles.X(new_edge_node2)=TV((T)1,(T)0);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    last_frame=1000;
}
//#####################################################################
// Function Two_Spring_Test
//#####################################################################
void Two_Spring_Test(const T midpoint_location)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    int new_edge_node3=segmented_curve.particles.array_collection->Add_Element();
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=(T)1;
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node3)=(T)1;
    segmented_curve.particles.X(new_edge_node1)=TV((T)-1,(T)0);segmented_curve.particles.X(new_edge_node2)=TV(midpoint_location,(T)0);segmented_curve.particles.X(new_edge_node3)=TV((T)1,(T)0);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node2,new_edge_node3));
    last_frame=1000;
}
//#####################################################################
// Function Three_Spring_Test
//#####################################################################
 void Three_Spring_Test(const T height)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    int new_edge_node3=segmented_curve.particles.array_collection->Add_Element();
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=(T)1;
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node3)=(T)1;
    segmented_curve.particles.X(new_edge_node1)=TV((T)0,(T)0);segmented_curve.particles.X(new_edge_node2)=TV((T)0.5, height);segmented_curve.particles.X(new_edge_node3)=TV((T)1,(T)0);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node2,new_edge_node3));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node3,new_edge_node1));
    last_frame=1000;
}
//#####################################################################
// Function Ten_Spring_Test
//#####################################################################
 void Ten_Spring_Test()
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    int new_edge_node3=segmented_curve.particles.array_collection->Add_Element();int new_edge_node4=segmented_curve.particles.array_collection->Add_Element();
    int new_edge_node5=segmented_curve.particles.array_collection->Add_Element();int new_edge_node6=segmented_curve.particles.array_collection->Add_Element();
    int new_edge_node7=segmented_curve.particles.array_collection->Add_Element();int new_edge_node8=segmented_curve.particles.array_collection->Add_Element();
    int new_edge_node9=segmented_curve.particles.array_collection->Add_Element();int new_edge_node10=segmented_curve.particles.array_collection->Add_Element();
    int new_edge_node11=segmented_curve.particles.array_collection->Add_Element();
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node3)=static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node4)=
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node5)=static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node6)=
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node7)=static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node8)=
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node9)=static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node10)=
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node11)=(T)1;
    segmented_curve.particles.X(new_edge_node1)=TV((T)0,(T)0);segmented_curve.particles.X(new_edge_node2)=TV((T)1,(T)0);segmented_curve.particles.X(new_edge_node3)=TV((T)2,(T)0);
    segmented_curve.particles.X(new_edge_node4)=TV((T)3,(T)0);segmented_curve.particles.X(new_edge_node5)=TV((T)4,(T)0);segmented_curve.particles.X(new_edge_node6)=TV((T)5,(T)0);
    segmented_curve.particles.X(new_edge_node7)=TV((T)6,(T)0);segmented_curve.particles.X(new_edge_node8)=TV((T)7,(T)0);segmented_curve.particles.X(new_edge_node9)=TV((T)8,(T)0);
    segmented_curve.particles.X(new_edge_node10)=TV((T)9,(T)0);segmented_curve.particles.X(new_edge_node11)=TV((T)10,(T)0);

    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node2,new_edge_node3));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node3,new_edge_node4));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node4,new_edge_node5));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node5,new_edge_node6));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node6,new_edge_node7));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node7,new_edge_node8));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node8,new_edge_node9));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node9,new_edge_node10));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node10,new_edge_node11));
    last_frame=1000;
}
//#####################################################################
// Function Constrained_Two_Spring_Test
//#####################################################################
void Constrained_Two_Spring_Test(const T midpoint_location)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(particles);solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
    int new_edge_node1=segmented_curve.particles.array_collection->Add_Element();int new_edge_node2=segmented_curve.particles.array_collection->Add_Element();
    int new_edge_node3=segmented_curve.particles.array_collection->Add_Element();
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node3)=FLT_MAX;
    static_cast<PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=(T)1;
    segmented_curve.particles.X(new_edge_node1)=TV((T)-1,(T)0);segmented_curve.particles.X(new_edge_node2)=TV(midpoint_location,(T)0);segmented_curve.particles.X(new_edge_node3)=TV((T)1,(T)0);
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
    segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node2,new_edge_node3));
    last_frame=1000;
}
//#####################################################################
// Function Single_Spring_Moving_To_Right_Test
//#####################################################################
void Single_Spring_Moving_To_Right_Test()
{
    Single_Spring_Test();
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    particles.V(1)=particles.V(2)=TV((T)1,(T)0);
}
//#####################################################################
// Function Spring_Mesh
//#####################################################################
void Spring_Mesh()
{
    int grid_m=100,grid_n=100;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    particles.array_collection->Add_Elements(grid_m*grid_n);

    SEGMENTED_CURVE<TV> *segmented_curve=SEGMENTED_CURVE<TV>::Create(particles);
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);    

    for(int i=0;i<grid_m;i++) for(int j=0;j<grid_n;j++){
        particles.X((i-1)*grid_n+j)=TV(i-2,j-2);
        if(i<grid_m) segmented_curve->mesh.elements.Append(VECTOR<int,2>((i-1)*grid_n+j,i*grid_n+j));
        if(j<grid_n) segmented_curve->mesh.elements.Append(VECTOR<int,2>((i-1)*grid_n+j,(i-1)*grid_n+j+1));}

    particles.mass.Fill((T)1);
}
//#####################################################################
// Function Hanging_Cloth_Test
//#####################################################################
void Hanging_Cloth_Test()
{
    int m=number_side_panels+1,n=(int)(aspect_ratio*number_side_panels)+1;
    PARTICLES<TV>& particles=*(new PARTICLES<TV>());
    particles.Store_Mass();

    TRIANGLE_MESH mesh;
    mesh.Initialize_Herring_Bone_Mesh(m,n);particles.array_collection->Add_Elements(mesh.number_nodes);
    T mass_node=aspect_ratio*sqr(side_length)/(m*n);particles.mass.Fill(mass_node);
    int i,j;
    i=1;j=n;particles.mass(i+m*(j-1))=FLT_MAX;i=m;j=n;particles.mass(i+m*(j-1))=FLT_MAX;
    T dx=side_length/(m-1),dy=aspect_ratio*side_length/(n-1);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) particles.X(i+m*(j-1))=TV((i-1)*dx,(j-1)*dy);
    mesh.Initialize_Segment_Mesh();

    SEGMENTED_CURVE<TV> *segmented_curve=SEGMENTED_CURVE<TV>::Create(particles);
    for(int s=1;s<=mesh.segment_mesh->elements.m;s++)
        segmented_curve->mesh.elements.Append(VECTOR<int,2>(mesh.segment_mesh->elements(s)(1),mesh.segment_mesh->elements(s)(2)));
    tests.Copy_And_Add_Structure(*segmented_curve,0);

    delete &particles;
}
//#####################################################################
// Function Create_Spring_Forces
//#####################################################################
LINEAR_SPRINGS<TV>* Create_Spring_Forces(T restlength=(T)0)
{
    SEGMENTED_CURVE_2D<T>& segmented_curve=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE_2D<T>&>();
    T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
    LINEAR_SPRINGS<TV>* spring_force=Create_Edge_Springs(segmented_curve,linear_stiffness,linear_damping);
    solid_body_collection.Add_Force(spring_force);
    if(restlength){
        ARRAY<T> restlengths(spring_force->segment_mesh.elements.m);restlengths.Fill(restlength);
        spring_force->Set_Restlength(restlengths);}
    return spring_force;
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==1) V(1)=TV();
    else if(test_number==7) V(1)=TV();
    else if(test_number==8) V(1)=V(3)=TV();
    else if(test_number==9){
        int i,j;int n=(int)(aspect_ratio*number_side_panels)+1,m=number_side_panels+1;
        i=1;j=n;V(i+m*(j-1))=TV();i=m;j=n;V(i+m*(j-1))=TV();}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==1) V(1)=TV();
    else if(test_number==7) V(1)=TV();
    else if(test_number==8) V(1)=V(3)=TV();
    else if(test_number==9){
        int i,j;int n=(int)(aspect_ratio*number_side_panels)+1,m=number_side_panels+1;
        i=1;j=n;V(i+m*(j-1))=TV();i=m;j=n;V(i+m*(j-1))=TV();}
}
//#####################################################################
// Function Set_Deformable_Particle_Is_Simulated
//#####################################################################
void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE
{
    if(test_number==1) particle_is_simulated(1)=false;
    else if(test_number==7) particle_is_simulated(1)=false;
    else if(test_number==8) particle_is_simulated(1)=particle_is_simulated(3)=false;
    else if(test_number==9){
        int i,j;int m=number_side_panels+1,n=(int)(aspect_ratio*number_side_panels)+1;
        i=1;j=n;particle_is_simulated(i+m*(j-1))=false;i=m;j=n;particle_is_simulated(i+m*(j-1))=false;}
}
//#####################################################################
};
}
#endif
