//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_2D
//#####################################################################
//    1. What is it
//#####################################################################
#ifndef __STANDARD_TESTS_2D__
#define __STANDARD_TESTS_2D__

#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Matrices/BANDED_MATRIX.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Implicit_Objects_Uniform/SMOOTH_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE_PATCH.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/LINEAR_BINDING.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_PENALTY.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Forces/B_SPLINE_CURVATURE_FORCE.h>
#include <Deformables/Forces/MOONEY_RIVLIN_B_SPLINE_CURVATURE_FORCE.h>
#include <Deformables/Forces/BEZIER_C2_FORCE.h>
#include <Deformables/Forces/BEZIER_CURVATURE_FORCE.h>
#include <Deformables/Forces/ELASTIC_ETHER_DRAG.h>
#include <Deformables/Forces/RALEIGH_DAMPING_FORCE.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_OBJECTIVE.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <fstream>
#include "STANDARD_TESTS_BASE.h"
namespace PhysBAM{

template<class TV> class STANDARD_TESTS;
template<class T_input>
class STANDARD_TESTS<VECTOR<T_input,2> >:public STANDARD_TESTS_BASE<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef STANDARD_TESTS_BASE<TV> BASE;
    using BASE::solid_body_collection;using BASE::stream_type;using BASE::solids_evolution;using BASE::parse_args;
    using BASE::test_number;using BASE::data_directory;using BASE::m;using BASE::s;using BASE::kg;
    using BASE::tests;using BASE::density;using BASE::Get_Initial_Data_After;using BASE::use_penalty_self_collisions;
    using BASE::Initialize_Bodies_After;using BASE::resolution;using BASE::stiffness_multiplier;
    using BASE::curvature_stiffness_multiplier;using BASE::point_curves;using BASE::kinematic_points;
    using BASE::thickness;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input)
        :BASE(stream_type_input)
    {
    }

    virtual ~STANDARD_TESTS()
    {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
//    parse_args->Add("-opt",&here,"meaning");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    bool automatically_add_to_collision_structures=true;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    switch(test_number){
        case 1:case 2:{
            if(resolution==0) resolution=20;
            ARRAY<TV> X(resolution+1);
            for(int i=0;i<=resolution;i++){
                T t=(T)i/resolution;
                X(i)=TV(t,sin(t*2*pi));}
            BEZIER_SPLINE<TV,3>* spline=BEZIER_SPLINE<TV,3>::Create(particles);
            Smooth_Fit<TV>(*spline,X);
            deformable_body_collection.Add_Structure(spline);
            deformable_body_collection.Add_Structure(Create_Segmented_Curve(*spline,true));
            deformable_body_collection.Add_Force(new BEZIER_CURVATURE_FORCE<TV>(particles,*spline,curvature_stiffness_multiplier,stiffness_multiplier));
            for(int i=0;i<spline->control_points.m;i++)
                particles.mass.Subset(spline->control_points(i))+=(T).25*density/spline->control_points.m;
            kinematic_points.Append(0);
            kinematic_points.Append(resolution*3);
            point_curves.Resize(2);
            point_curves(0).Add_Control_Point(0,particles.X(0));
            point_curves(0).Add_Control_Point(1,particles.X(0));
            point_curves(1).Add_Control_Point(0,particles.X(resolution*3));
            point_curves(1).Add_Control_Point(1,particles.X(resolution*3)+TV(test_number==1?1:-.9,0));
            for(int i=4;i<particles.X.m;i+=3)
                deformable_body_collection.binding_list.Add_Binding(
                    new LINEAR_BINDING<TV,2>(particles,i,VECTOR<int,2>(i-2,i-1),VECTOR<T,2>(-1,2)));
//            particles.mass(0)=FLT_MAX;
//            particles.mass(resolution*3)=FLT_MAX;
            break;}
        case 3:case 4:{
            if (resolution==0) resolution=20;
            ARRAY<TV> X(resolution+1);
            for(int i=0;i<=resolution;i++){
                T t=(T)i/resolution;
                X(i)=TV(t,sin(t*2*pi));}
            B_SPLINE<TV,3>* spline=B_SPLINE<TV,3>::Create(particles);
            Smooth_Fit<TV>(*spline,X);
            deformable_body_collection.Add_Structure(spline);
            deformable_body_collection.Add_Structure(Create_Segmented_Curve(*spline,true));
            deformable_body_collection.Add_Force(new B_SPLINE_CURVATURE_FORCE<TV>(particles,*spline,curvature_stiffness_multiplier,stiffness_multiplier));
            int lastpt=spline->control_points.m-1;
            for(int i=2;i<lastpt;i++){
                T toAdd=(T)density/(24*resolution);
                particles.mass(i-2)+=toAdd;
                particles.mass(i-1)+=11*toAdd;
                particles.mass(i)+=11*toAdd;
                particles.mass(i+1)+=toAdd;}
            kinematic_points.Append(0);
            kinematic_points.Append(lastpt);
            point_curves.Resize(2);
            point_curves(0).Add_Control_Point(0,particles.X(0));
            point_curves(0).Add_Control_Point(1,particles.X(0));
            point_curves(1).Add_Control_Point(0,particles.X(lastpt));
            point_curves(1).Add_Control_Point(1,particles.X(lastpt)+TV(test_number==3?1:4,0));
//            particles.mass(0)=FLT_MAX;
//            particles.mass(lastpt)=FLT_MAX;
            break;}
        case 13:case 14:{
            if (resolution==0) resolution=20;
            ARRAY<TV> X(resolution+1);
            for(int i=0;i<=resolution;i++){
                T t=(T)i/resolution;
                X(i)=TV(t,sin(t*2*pi));}
            B_SPLINE<TV,3>* spline=B_SPLINE<TV,3>::Create(particles);
            Smooth_Fit<TV>(*spline,X);
            deformable_body_collection.Add_Structure(spline);
            deformable_body_collection.Add_Structure(Create_Segmented_Curve(*spline,true));
            deformable_body_collection.Add_Force(new MOONEY_RIVLIN_B_SPLINE_CURVATURE_FORCE<TV>(particles,*spline,thickness/2,stiffness_multiplier));
            int lastpt=spline->control_points.m-1;
            for(int i=2;i<lastpt;i++){
                T toAdd=(T)density/(24*resolution);
                particles.mass(i-2)+=toAdd;
                particles.mass(i-1)+=11*toAdd;
                particles.mass(i)+=11*toAdd;
                particles.mass(i+1)+=toAdd;}
            kinematic_points.Append(0);
            kinematic_points.Append(lastpt);
            point_curves.Resize(2);
            point_curves(0).Add_Control_Point(0,particles.X(0));
            point_curves(0).Add_Control_Point(1,particles.X(0));
            point_curves(1).Add_Control_Point(0,particles.X(lastpt));
            point_curves(1).Add_Control_Point(1,particles.X(lastpt)+TV(test_number==13?1:-0.9,0));
//            particles.mass(0)=FLT_MAX;
//            particles.mass(lastpt)=FLT_MAX;
            break;}
        case 6:{
//            ARRAY<TV> curveX(5);
//            for(int i=0;i<5;i++)
//                curveX(i)=TV(i,0);
//            B_SPLINE<TV,3>* spline=B_SPLINE<TV,3>::Create(particles);
//            Smooth_Fit<TV>(*spline,curveX);
//            LOG::printf("curve version gives %P\n",particles.X);

            typedef VECTOR<T,3> TV3;
            int m=5; int n=4;
            ARRAY<ARRAY<TV3> > X(m);



            RANDOM_NUMBERS<T> random;random.Set_Seed(1823);

            for(int i=0; i<m; i++)
            {
                X(i).Resize(n);
                for(int j=0; j<n;j++)
                    X(i)(j)=TV3(random.Get_Uniform_Number(-(T)10,(T)10),random.Get_Uniform_Number(-(T)10,(T)10),random.Get_Uniform_Number(-(T)10,(T)10));
//                    X(i)(j)=TV3(i,j,i*j);
//                    X(i)(j)=TV3(i,0,0);
//                    X(i)(j)=TV3(1,2,3);
            }

            LOG::printf("targets: %P\n",X);

            B_SPLINE_PATCH<TV3,3>* patch=B_SPLINE_PATCH<TV3,3>::Create();
            Smooth_Fit<TV3>(*patch,X);

            for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                    LOG::printf("i: %P j: %P val: %P\n\n",i,j,patch->Evaluate((T)i/(m-1),(T)j/(n-1)));

            BEZIER_SPLINE_PATCH<TV3,3> bezpatch;
            Fill_Bezier(bezpatch,*patch);

            T s=.6;
            T t=.28;
            int d=3;
            LOG::printf("s: %P t: %P\n",s,t);
            LOG::printf("B_SPLINE value: %P\n",patch->Evaluate(s,t));
            
            ARRAY<ARRAY<TV3> > Xt(n); // Transpose of X.
            for(int i=0;i<n;i++)
            {
                Xt(i).Resize(m);
                for(int j=0;j<m;j++)
                    Xt(i)(j)=X(j)(i);
            }

            B_SPLINE_PATCH<TV3,3>* patcht = B_SPLINE_PATCH<TV3,3>::Create();
            Smooth_Fit<TV3>(*patcht,Xt);


            LOG::printf("B_SPLINE transpose's value: %P\n",patcht->Evaluate(t,s));





            
            t=clamp(t,patch->knots_t(d-1),patch->knots_t(patch->knots_t.m-d));
            int id_t=std::upper_bound(patch->knots_t.begin(),patch->knots_t.end()-d,t)-patch->knots_t.begin();
            s=clamp(s,patch->knots_s(d-1),patch->knots_s(patch->knots_s.m-d));
            int id_s=std::upper_bound(patch->knots_s.begin(),patch->knots_s.end()-d,s)-patch->knots_s.begin();
            id_t-=3;
            id_s-=3;
            int id=id_s*(patch->knots_t.m-5) + id_t;
            T bez_t=(t-patch->knots_t(id_t+2))/(patch->knots_t(id_t+3)-patch->knots_t(id_t+2));
            T bez_s=(s-patch->knots_s(id_s+2))/(patch->knots_s(id_s+3)-patch->knots_s(id_s+2));

            
            LOG::printf("BEZIER value: %P\n",bezpatch.Evaluate(id,bez_s,bez_t));

            LOG::printf("\n\nSecond derivative testing:\n");
            s=random.Get_Uniform_Number((T)0,(T)1);
            t=random.Get_Uniform_Number((T)0,(T)1);
            
            T hh = .1;
            for(int i=0;i<6;i++)
            {
                LOG::printf("h: %P\n",hh);
                LOG::printf("(s,0): %P\n",(patch->Evaluate(s,2*hh)-(T)2*patch->Evaluate(s,hh)+patch->Evaluate(s,0))/(hh*hh));
                LOG::printf("(0,t): %P\n",(patch->Evaluate(2*hh,t)-(T)2*patch->Evaluate(hh,t)+patch->Evaluate(0,t))/(hh*hh));
                LOG::printf("(s,1): %P\n",(patch->Evaluate(s,1-2*hh)-(T)2*patch->Evaluate(s,1-hh)+patch->Evaluate(s,1))/(hh*hh));
                LOG::printf("(1,t): %P\n",(patch->Evaluate(1-2*hh,t)-(T)2*patch->Evaluate(1-hh,t)+patch->Evaluate(1,t))/(hh*hh));
                hh /= (T)10;
            }

            LOG::printf("\n\n");
            s=random.Get_Uniform_Number((T)0,(T)1);
            t=random.Get_Uniform_Number((T)0,(T)1);
            
            hh = (T)0.1;
            for(int i=0;i<6;i++)
            {
                LOG::printf("h: %P\n",hh);
                LOG::printf("(s,0): %P\n",(patch->Evaluate(s,2*hh)-(T)2*patch->Evaluate(s,hh)+patch->Evaluate(s,0))/(hh*hh));
                LOG::printf("(0,t): %P\n",(patch->Evaluate(2*hh,t)-(T)2*patch->Evaluate(hh,t)+patch->Evaluate(0,t))/(hh*hh));
                LOG::printf("(s,1): %P\n",(patch->Evaluate(s,1-2*hh)-(T)2*patch->Evaluate(s,1-hh)+patch->Evaluate(s,1))/(hh*hh));
                LOG::printf("(1,t): %P\n",(patch->Evaluate(1-2*hh,t)-(T)2*patch->Evaluate(1-hh,t)+patch->Evaluate(1,t))/(hh*hh));
                hh /= (T)10;
            }
//            for(int i=0;i<m;i++)
//                for(int j=0;j<n;j++)
//                    LOG::printf("i: %P j: %P val: %P\n\n",i,j,bezpatch.Evaluate((T)i/(m-1),(T)j/(n-1)));

            break;}
//        case 7:{
//            SPARSE_MATRIX_FLAT_MXN<T> mat;
//            typedef MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<T>,T,KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > > SYSTEM;
//            KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > x,rhs;
//            mat.m=2;
//            mat.n=2;
//            x.v.Resize(2);
//            rhs.v.Resize(2);
//            
//            mat(0,0)=0;
//            mat(1,0)=1;
//            mat(0,1)=1;
//            mat(1,1)=0;
//
//            rhs.v(0)=3;
//            rhs.v(1)=5;
//// Fill rhs and mat!
//// use rhs.v(i), right?
//            
//            
//            
//            
//            mat.Construct_Incomplete_Cholesky_Factorization();
//            SYSTEM system(mat);
//            system.P=mat.C;
//            
//    
//            CONJUGATE_GRADIENT<T> cg;
//            cg.print_diagnostics=true;
//            bool result=cg.Solve(system,x,rhs,vectors,(T)1e-4,1,1000);
//            LOG::printf("solution is %P\n",rhs.v);
//            PHYSBAM_ASSERT(result);
//            
//            break;}
        default:
            LOG::cerr<<"Initial Data: Unrecognized test number "<<test_number<<std::endl;exit(1);}

    Get_Initial_Data_After(automatically_add_to_collision_structures);

    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 4:
        case 13:
        case 14:
            break;
        default:
            LOG::cerr<<"Missing bodies implementation for test number "<<test_number<<std::endl;exit(1);}

    Initialize_Bodies_After();
}
};
}
#endif
