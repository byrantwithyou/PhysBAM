//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
#include <Deformables/Forces/OPENSUBDIV_SURFACE_CURVATURE_FORCE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include "STANDARD_TESTS_3D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
STANDARD_TESTS(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type,parse_args)
{
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
~STANDARD_TESTS()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Write_Output_Files(const int frame)
{
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Read_Output_Files(const int frame)
{
    BASE::Read_Output_Files(frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Initialize()
{
    switch(test_number)
    {
        case 1:{ // rotating sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5,.5),.3);
            VECTOR<T,3> angular_velocity(TV(0.4,0,0));
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity);}
                ,density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            Add_Fixed_Corotated(1e3*scale_E,0.3);
        } break;
        case 2:{ // Oscillating sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,3>()+1.5);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
        } break;
        case 3:{ // Freefall sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
            Add_Gravity(TV(0,-9.8,0));
        } break;
        case 4:{ // subdivision surface - 40x40 strip
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-3,-3,-3),TV(3,3,3)),true);
            
            
            std::string filename=data_directory+"/OpenSubdiv_Surfaces/strip_40.dat.gz";
            OPENSUBDIV_SURFACE<TV>* surf=OPENSUBDIV_SURFACE<TV>::Create();
            T thickness=1e-3;
            surf->Initialize(filename,thickness);

            T density=1000*scale_mass;
//            surf->Set_Masses(density,thickness);

            T c1=3.5e6;
            T c2=1.3e6;
            T stiffness_multiplier=scale_E;
            T thickness_multiplier=1;

            OPENSUBDIV_SURFACE<TV>& new_surf=Seed_Lagrangian_Particles(*surf,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,3>();},density,false);
            MOONEY_RIVLIN_CURVATURE<T> model=MOONEY_RIVLIN_CURVATURE<T>(c1*stiffness_multiplier,c2*stiffness_multiplier,thickness*thickness_multiplier);
            Add_Force(*new OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,3>(particles,new_surf,model));

            Add_Walls(8,COLLISION_TYPE::separate,.3,.1,false);

            Add_Gravity(TV(0,-9.8,0));
        } break;
        case 5:{ // subdivision surface - duck
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-3,-3,-3),TV(3,3,3)),true);
            
            std::string filename=data_directory+"/OpenSubdiv_Surfaces/duck_1073f.dat.gz";
            OPENSUBDIV_SURFACE<TV>* surf=OPENSUBDIV_SURFACE<TV>::Create();
            T thickness=1e-3;
            surf->Initialize(filename,thickness);
            
            T density=1000*scale_mass;

            T c1=3.5e6;
            T c2=1.3e6;
            T stiffness_multiplier=scale_E;
            T thickness_multiplier=1;

            OPENSUBDIV_SURFACE<TV>& new_surf=Seed_Lagrangian_Particles(*surf,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,3>();},density,false);
            MOONEY_RIVLIN_CURVATURE<T> model=MOONEY_RIVLIN_CURVATURE<T>(c1*stiffness_multiplier,c2*stiffness_multiplier,thickness*thickness_multiplier);
            Add_Force(*new OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,3>(particles,new_surf,model));

//            for(int p=0;p<particles.X.m;p++) // squish the duck.
//                particles.X(p)(0)*=0.5;

            Add_Walls(8,COLLISION_TYPE::separate,.3,.1,false);

            Add_Gravity(TV(0,-9.8,0));
        } break;
        case 6:{ // subdivision surface - drop several ducks.
            int num_duckies=3;
            grid.Initialize(resolution*TV_INT(2,num_duckies,2),RANGE<TV>(TV(-6,-3,-6),TV(6,3+6*(num_duckies-1),6)),true);
            
            std::string filename=data_directory+"/OpenSubdiv_Surfaces/duck_1073f.dat.gz";

            T c1=3.5e6;
            T c2=1.3e6;
            T stiffness_multiplier=scale_E;
            T thickness_multiplier=1;
            T density=1000;
            T thickness=1e-3;

            MOONEY_RIVLIN_CURVATURE<T> model=MOONEY_RIVLIN_CURVATURE<T>(c1*stiffness_multiplier,c2*stiffness_multiplier,thickness*thickness_multiplier);

            OPENSUBDIV_SURFACE<TV>* surf=OPENSUBDIV_SURFACE<TV>::Create();
            surf->Initialize(filename,thickness);
            for(int i=0;i<num_duckies;i++){
                for(int p=0;p<particles.X.m;p++)
                    particles.X(p)(1)+=6;
                OPENSUBDIV_SURFACE<TV>& new_surf=Seed_Lagrangian_Particles(*surf,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,3>();},density,false,false);
                Add_Force(*new OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,3>(particles,new_surf,model));
            }

            Add_Walls(8,COLLISION_TYPE::separate,.3,.1,false);

            Add_Gravity(TV(0,-9.8,0));
            delete surf;
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
}
//#####################################################################
// Function Begin_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Begin_Frame(const int frame)
{
}
//#####################################################################
// Function End_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
End_Frame(const int frame)
{
}
//#####################################################################
// Function Begin_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Begin_Time_Step(const T time)
{
}
//#####################################################################
// Function End_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
End_Time_Step(const T time)
{
}
template class STANDARD_TESTS<VECTOR<float,3> >;
template class STANDARD_TESTS<VECTOR<double,3> >;
}
