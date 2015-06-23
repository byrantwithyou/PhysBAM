//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//    1. What is it
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Implicit_Objects/SMOOTH_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE_PATCH.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_PENALTY.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
#include <Deformables/Forces/B_SPLINE_PATCH_CURVATURE_FORCE.h>
#include <Deformables/Forces/ELASTIC_ETHER_DRAG.h>
#include <Deformables/Forces/OPENSUBDIV_SURFACE_CURVATURE_FORCE.h>
#include <Deformables/Forces/RALEIGH_DAMPING_FORCE.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_OBJECTIVE.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <fstream>
#include <string>
#include "STANDARD_TESTS_BASE.h"
namespace PhysBAM{

template<class TV> class STANDARD_TESTS;
template<class T_input>
class STANDARD_TESTS<VECTOR<T_input,3> >:public STANDARD_TESTS_BASE<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<int,2> IV;
public:
    typedef STANDARD_TESTS_BASE<TV> BASE;
    using BASE::solid_body_collection;using BASE::stream_type;using BASE::solids_evolution;
    using BASE::test_number;using BASE::data_directory;using BASE::m;using BASE::s;using BASE::kg;
    using BASE::tests;using BASE::density;using BASE::Get_Initial_Data_After;using BASE::use_penalty_self_collisions;
    using BASE::Initialize_Bodies_After;using BASE::point_curves;using BASE::kinematic_points;
    using BASE::thickness_multiplier;using BASE::stiffness_multiplier;using BASE::resolution;
    using BASE::rand;using BASE::gauss_order;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args)
    {
        parse_args.Parse();
    }

    virtual ~STANDARD_TESTS()
    {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    bool automatically_add_to_collision_structures=true;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    switch(test_number){
        case 997:{
            OPENSUBDIV_SURFACE<TV>* surf=OPENSUBDIV_SURFACE<TV>::Create(particles);

            std::string filename=data_directory+"/OpenSubdiv_Surfaces/strip_40.dat.gz";
            int res=40;
            
            T thickness=1e-3;
            surf->Initialize(filename,thickness*thickness_multiplier);
            surf->Set_Mass(density);

            deformable_body_collection.Add_Structure(surf);
            deformable_body_collection.Add_Structure(Create_Triangulated_Surface(*surf,true));
            
            T c1=3.5e6;
            T c2=1.3e6;

            MOONEY_RIVLIN_CURVATURE<T> model=MOONEY_RIVLIN_CURVATURE<T>(c1*stiffness_multiplier,c2*stiffness_multiplier,thickness*thickness_multiplier);
            deformable_body_collection.Add_Force(new OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,3>(particles,*surf,model));
//            // 3 is the gauss order.


//            for(int i=0;i<res;i++){
//                int p=i*res;
//                kinematic_points.Append(p);
//                INTERPOLATION_CURVE<T,TV> c;
//                c.Add_Control_Point(0,particles.X(p));
//            point_curves.Append(c);}

            int M=res;
            int N=res;
            int first0=M/3;
            int last0=2*M/3;
            int boundary_thickness=ceil((T)N/10.0);
            int last1=N;
            for(int i=0;i<=last0-first0;i++){
                for(int j=0;j<boundary_thickness;j++){
                    int p=surf->control_points((i+first0)*M+j);
                    kinematic_points.Append(p);
                    INTERPOLATION_CURVE<T,TV> c;
                    c.Add_Control_Point(0,particles.X(p));
                    point_curves.Append(c);}
                for(int j=last1-boundary_thickness;j<last1;j++){
                    int p=surf->control_points((i+first0)*M+j);
                    kinematic_points.Append(p);
                    INTERPOLATION_CURVE<T,TV> c;
                    c.Add_Control_Point(0,particles.X(p));
//                    c.Add_Control_Point(0,particles.X(p)+m*TV(0,0,0.01));
                    c.Add_Control_Point(1,particles.X(p)+m*TV(0,0,0.01));
//                    c.Add_Control_Point(2,particles.X(p)-m*TV(0,0,.25));
                    
                    c.Add_Control_Point(5,particles.X(p)-m*TV(0,0,1));
////                    c.Add_Control_Point(0,particles.X(p)+m*TV(0,0,0.01));
////                    c.Add_Control_Point(1,particles.X(p)+m*TV(0,0,0.01));
////                    c.Add_Control_Point(5,particles.X(p)-m*TV(0,0,1));
                    point_curves.Append(c);}}
            
            
            Add_Gravity();
            automatically_add_to_collision_structures=false;
            break;}
        case 996:{// let part of a duck drape itself.
            OPENSUBDIV_SURFACE<TV>* surf=OPENSUBDIV_SURFACE<TV>::Create(particles);
            std::string filename=data_directory+"/OpenSubdiv_Surfaces/duck_1073f.dat.gz";
            T thickness=1e-3;
            surf->Initialize(filename,thickness*thickness_multiplier);
            surf->Set_Mass(density);
            deformable_body_collection.Add_Structure(surf);
            deformable_body_collection.Add_Structure(Create_Triangulated_Surface(*surf,true));
            T c1=3.5e6;
            T c2=1.3e6;
            MOONEY_RIVLIN_CURVATURE<T> model=MOONEY_RIVLIN_CURVATURE<T>(c1*stiffness_multiplier,c2*stiffness_multiplier,thickness*thickness_multiplier);
            deformable_body_collection.Add_Force(new OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,3>(particles,*surf,model));

            // let's hold the top part of the duck fixed. // now try moving the top part. // and stretching too?
            for(int p=0;p<particles.X.m;p++){
                int cp=surf->control_points(p);
                if(particles.X(cp)(1)>.25){
                    kinematic_points.Append(cp);
                    INTERPOLATION_CURVE<T,TV> c;
                    c.Add_Control_Point(0,particles.X(cp));
                    c.Add_Control_Point(1,particles.X(cp));
                    c.Add_Control_Point(5,particles.X(cp)-m*TV(1,0,0));
//                    c.Add_Control_Point(5,0.5*particles.X(cp)-m*TV(1,0,0));
                    point_curves.Append(c);}}
            Add_Gravity();
            automatically_add_to_collision_structures=false;
            break;}
        case 1:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)25,0)*m)),true,true,density,m);
            tests.Add_Ground(0,1.99*m);
            break;}
        case 2:{  //Flattened curled ribbon
            int M=10;
            if(resolution)
                M=resolution;
            int N=2*pi*M;
            ARRAY<TV,IV> X(IV(M+1,N));
            for(int i=0;i<=M;i++){
                for(int j=0;j<N;j++)
                    X(i,j)=m*TV((T)i/M,sin((2*j+1)*pi/N),-cos((2*j+1)*pi/N));}
            Add_Patch(X);
            for(int i=0;i<particles.X.m;i++){
                particles.X(i)(2)=atan2(particles.X(i)(1),particles.X(i)(2));
                particles.X(i)(1)*=0;}

            automatically_add_to_collision_structures=false;
            break;}
        case 3:{  //Random Test
            int M=10,N=10;
            if(resolution)
                M=N=resolution;
            Add_Patch(M,N);
            rand.Fill_Uniform(particles.X,-1,1);
            automatically_add_to_collision_structures=false;
            break;}
        case 4:{
            int M=10;
            if(resolution)
                M=resolution;
            int N=(pi*M)/2;
            ARRAY<TV,IV> X(IV(M+1,N));
            for(int i=0;i<=M;i++){
                for(int j=0;j<N;j++)
                    X(i,j)=m*TV((T)4*i/M-2,sin((2*j+1)*pi/N),-cos((2*j+1)*pi/N));}
            B_SPLINE_PATCH<TV,3>* patch=Add_Patch(X,false,true);

            Add_Gravity();
            
            int boundary_thickness=ceil(patch->control_points.domain.max_corner(0)/10.0);
            int last0=patch->control_points.domain.max_corner(0);
            int last1=patch->control_points.domain.max_corner(1)-3;
            for(int j=0;j<last1;j++){
                for(int i=0;i<boundary_thickness;i++){
                    int p=patch->control_points(i,j);
                    kinematic_points.Append(p);
                    INTERPOLATION_CURVE<T,TV> c;
                    c.Add_Control_Point(0,particles.X(p));
                    point_curves.Append(c);}
                for(int i=last0-boundary_thickness;i<last0;i++){
                    int p=patch->control_points(i,j);
                    kinematic_points.Append(p);
                    INTERPOLATION_CURVE<T,TV> c;
                    c.Add_Control_Point(0,particles.X(p)+m*TV(0.01,0,0));
                    c.Add_Control_Point(1,particles.X(p)+m*TV(0.01,0,0));
                    c.Add_Control_Point(5,particles.X(p)-m*TV(1,0,0));
                    point_curves.Append(c);}}

            automatically_add_to_collision_structures=false;
            break;}
        case 7:{ // Fold napkin
            int M=10,N=10;
            if(resolution)
                M=N=resolution;
            B_SPLINE_PATCH<TV,3>* patch=Add_Patch(M,N);
            
            int last0=patch->control_points.domain.max_corner(0)-1;
            int last1=patch->control_points.domain.max_corner(1)-1;
            int i00=patch->control_points(0,0);
            int i01=patch->control_points(0,last1);
            int i10=patch->control_points(last0,0);
            int i11=patch->control_points(last0,last1);
            kinematic_points.Append(i00);
            kinematic_points.Append(i01);
            kinematic_points.Append(i10);
            kinematic_points.Append(i11);
            point_curves.Resize(4);
            point_curves(0).Add_Control_Point(0,particles.X(i00));
            point_curves(0).Add_Control_Point(4,particles.X(i00)+m*TV(0.5,0.5,0.5));
            point_curves(1).Add_Control_Point(0,particles.X(i01));
            point_curves(1).Add_Control_Point(5,particles.X(i01)+m*TV(0.5,-0.5,-0.5));
            point_curves(2).Add_Control_Point(0,particles.X(i10));
            point_curves(2).Add_Control_Point(5,particles.X(i10)+m*TV(-0.5,-0.5,0.5));
            point_curves(3).Add_Control_Point(0,particles.X(i11));
            point_curves(3).Add_Control_Point(4,particles.X(i11)+m*TV(-0.5,0.5,-0.5));
            
            break;}
        case 8:{ // curl corners
            int M=10,N=10;
            if(resolution)
                M=N=resolution;
            B_SPLINE_PATCH<TV,3>* patch=Add_Patch(M,N);
            
            int last0=patch->control_points.domain.max_corner(0)-1;
            int last1=patch->control_points.domain.max_corner(1)-1;
            int i00=patch->control_points(0,0);
            int i01=patch->control_points(0,last1);
            int i10=patch->control_points(last0,0);
            int i11=patch->control_points(last0,last1);
            kinematic_points.Append(i00);
            kinematic_points.Append(i01);
            kinematic_points.Append(i10);
            kinematic_points.Append(i11);
            point_curves.Resize(4);
            point_curves(0).Add_Control_Point(0,particles.X(i00));
            point_curves(0).Add_Control_Point(8,particles.X(i00)+m*TV(0.45,0.1,0.45));
            point_curves(1).Add_Control_Point(0,particles.X(i01));
            point_curves(2).Add_Control_Point(0,particles.X(i10));
            point_curves(3).Add_Control_Point(0,particles.X(i11));
            point_curves(3).Add_Control_Point(8,particles.X(i11)+m*TV(-0.45,0.1,-0.45));
            
            automatically_add_to_collision_structures=false;
            break;}
        case 9:{ // bend
            int M=10;
            if(resolution)
                M=resolution;
            int N=2*M;
            B_SPLINE_PATCH<TV,3>* patch=Add_Patch(M,N);
            Add_Gravity();
            
            int first0=patch->control_points.domain.max_corner(0)/3;
            int last0=(2*patch->control_points.domain.max_corner(0))/3;
            int boundary_thickness=ceil(patch->control_points.domain.max_corner(1)/10.0);
            int last1=patch->control_points.domain.max_corner(1);
            for(int i=0;i<=last0-first0;i++){
                for(int j=0;j<boundary_thickness;j++){
                    int p=patch->control_points(i+first0,j);
                    kinematic_points.Append(p);
                    INTERPOLATION_CURVE<T,TV> c;
                    c.Add_Control_Point(0,particles.X(p));
                    point_curves.Append(c);}
                for(int j=last1-boundary_thickness;j<last1;j++){
                    int p=patch->control_points(i+first0,j);
                    kinematic_points.Append(p);
                    INTERPOLATION_CURVE<T,TV> c;
                    c.Add_Control_Point(0,particles.X(p)+m*TV(0,0,0.01));
                    c.Add_Control_Point(1,particles.X(p)+m*TV(0,0,0.01));
                    c.Add_Control_Point(5,particles.X(p)-m*TV(0,0,1));
                    point_curves.Append(c);}}
            
            automatically_add_to_collision_structures=false;
            break;}
        default:
            LOG::cerr<<"Initial Data: Unrecognized test number "<<test_number<<std::endl;exit(1);}

    Get_Initial_Data_After(automatically_add_to_collision_structures);

    switch(test_number){
        case 1: case 2: case 3: case 4: case 7: case 8: case 9: case 997: case 996:{
//                for(int i=0;i<kinematic_points.m;i++)
//                    particles.mass(kinematic_points(i))=FLT_MAX;
                 use_penalty_self_collisions=false;}
               // Fallthrough
        break;
        default:
            LOG::cerr<<"Missing bodies implementation for test number "<<test_number<<std::endl;exit(1);}

    Initialize_Bodies_After();
}

B_SPLINE_PATCH<TV,3>* Add_Patch(int M,int N,T c1=3.5e6,T c2=1.3e6,T thickness=1e-3){
    ARRAY<TV,IV> X(IV(M+1,N+1));
    for(int i=0;i<=M;i++){
        for(int j=0;j<=N;j++)
            X(i,j)=m*TV((T)i/M,0,(T)j/M);}

    return Add_Patch(X,c1,c2,thickness);
}

B_SPLINE_PATCH<TV,3>* Add_Patch(const ARRAY<TV,IV>& X,bool loop_s=false,bool loop_t=false,T c1=3.5e6,T c2=1.3e6,T thickness=1e-3)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    B_SPLINE_PATCH<TV,3>* patch=B_SPLINE_PATCH<TV,3>::Create(particles);
    Smooth_Fit<TV>(*patch,X,loop_s,loop_t);

    deformable_body_collection.Add_Structure(patch);
    MOONEY_RIVLIN_CURVATURE<T> model=MOONEY_RIVLIN_CURVATURE<T>(c1*stiffness_multiplier,c2*stiffness_multiplier,thickness*thickness_multiplier);
    // B_SPLINE_CURVATURE_FORCE will set the mass on the first Update_Position_Based_State call
    switch(gauss_order){
        case 1: deformable_body_collection.Add_Force(new B_SPLINE_PATCH_CURVATURE_FORCE<T,1>(particles,*patch,model,density)); break;
        case 2: deformable_body_collection.Add_Force(new B_SPLINE_PATCH_CURVATURE_FORCE<T,2>(particles,*patch,model,density)); break;
        case 3: deformable_body_collection.Add_Force(new B_SPLINE_PATCH_CURVATURE_FORCE<T,3>(particles,*patch,model,density)); break;
        case 4: deformable_body_collection.Add_Force(new B_SPLINE_PATCH_CURVATURE_FORCE<T,4>(particles,*patch,model,density)); break;
        case 5: deformable_body_collection.Add_Force(new B_SPLINE_PATCH_CURVATURE_FORCE<T,5>(particles,*patch,model,density)); break;
        case 6: deformable_body_collection.Add_Force(new B_SPLINE_PATCH_CURVATURE_FORCE<T,6>(particles,*patch,model,density)); break;
        case 7: deformable_body_collection.Add_Force(new B_SPLINE_PATCH_CURVATURE_FORCE<T,7>(particles,*patch,model,density)); break;
        default: LOG::cerr<<"Unsupported gauss order "<<gauss_order<<std::endl;exit(1);}

    return patch;
}

GRAVITY<TV>& Add_Gravity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    GRAVITY<TV>* g=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true);
    solid_body_collection.Add_Force(g);
    return *g;
}

};
}
#endif
