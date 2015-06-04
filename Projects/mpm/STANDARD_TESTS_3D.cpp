//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TORUS.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
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
        case 7:{ // skew impact of two elastic spheres
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(),TV(30,30,30)),true);
            T density=5*scale_mass;
            SPHERE<TV> sphere1(TV(10,13,15),2);
            Seed_Particles_Helper(sphere1,[=](const TV& X){return TV(0.75,0,0);},[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            SPHERE<TV> sphere2(TV(20,15,15),2);
            Seed_Particles_Helper(sphere2,[=](const TV& X){return TV(-0.75,0,0);},[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            Add_Neo_Hookean(31.685*scale_E,0.44022); //solve({E/(2*(1+r))=11,E*r/((1+r)*(1-2*r))=81},{E,r});
        } break;

        case 22:{ // (fluid test) pool of water 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(TV(0,0,0),TV(1,0.25,1));
            T density=2*scale_mass;
            Seed_Particles(box,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            Add_Gravity(TV(0,-9.8,0));
        } break;
        case 24:{ // (fluid test) circle drop 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.75,.5),.2);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            Add_Gravity(TV(0,-9.8,0));
        } break;

        case 8:{ // torus into a box
            // TODO: fix crash ./mpm -3d 8 -affine -last_frame 500 -midpoint -resolution 40 -newton_tolerance 1e-3 
            grid.Initialize(TV_INT(resolution,resolution*2,resolution)+1,RANGE<TV>(TV(),TV(1,2,1)),true);

            // Add_Walls(8,COLLISION_TYPE::separate,.3,.1,false);

            T thickness=.1;
            Add_Collision_Object(RANGE<TV>(TV(.2,0,.2),TV(.8,.1,.8)),COLLISION_TYPE::separate,.3);
            Add_Collision_Object(RANGE<TV>(TV(.2-thickness,0,0.2-thickness),TV(.8+thickness,.4,.2)),COLLISION_TYPE::separate,.3);
            Add_Collision_Object(RANGE<TV>(TV(.2-thickness,0,0.8),TV(.8+thickness,.4,0.8+thickness)),COLLISION_TYPE::separate,.3);
            Add_Collision_Object(RANGE<TV>(TV(0.2-thickness,0,0.2),TV(.2,.4,.8)),COLLISION_TYPE::separate,.3);
            Add_Collision_Object(RANGE<TV>(TV(.8,0,0.2),TV(.8+thickness,.4,.8)),COLLISION_TYPE::separate,.3);

            T density=5*scale_mass;
            ANALYTIC_IMPLICIT_OBJECT<TORUS<T> > torus(TORUS<T>(TV(),TV(1,0,0),0.02*2,0.03*2));
            Seed_Particles(torus,0,0,density,particles_per_cell);
            int m=particles.number;
            GRID<TV> torus_grid(TV_INT(2,3,2),RANGE<TV>(TV(.4,1,.4),TV(.6,1.5,.6)));
            for(NODE_ITERATOR<TV> iterator(torus_grid);iterator.Valid();iterator.Next()){
                TV center=iterator.Location();
                T angle=random.Get_Uniform_Number((T)0,(T)pi*2);
                ROTATION<TV> rotation(angle,TV(0,1,0));
                for(int k=0;k<m;k++)
                    Add_Particle(center+rotation.Rotate(particles.X(k)),0,0,particles.mass(k),particles.volume(k));
                break;
            }
            for(int k=0;k<m;k++) particles.Add_To_Deletion_List(k);
            particles.Delete_Elements_On_Deletion_List(true);
            Add_Fixed_Corotated(150,0.4);
            Add_Gravity(TV(0,-9.8,0));
        } break;

        case 9:{ // torus into a bowl
            // TODO: fix crash in LEVELSET_MAKER_UNIFORM ./mpm -3d 9 -resolution 80
            grid.Initialize(TV_INT(resolution,resolution,resolution)+1,RANGE<TV>(TV(-2,-2,-2),TV(2,2,2)),true);

            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/Rigid_Bodies/bowl.tri.gz",*surface);
            LOG::cout<<"Read mesh "<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh "<<surface->particles.number<<std::endl;
            surface->mesh.Initialize_Adjacent_Elements();    
            surface->mesh.Initialize_Neighbor_Nodes();
            surface->mesh.Initialize_Incident_Elements();
            surface->Update_Bounding_Box();
            surface->Initialize_Hierarchy();
            surface->Update_Triangle_List();
            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Initialize_Implicit_Surface(*surface,300);
            Add_Collision_Object(levelset,COLLISION_TYPE::separate,.3);

            T density=5*scale_mass;
            ANALYTIC_IMPLICIT_OBJECT<TORUS<T> > torus(TORUS<T>(TV(),TV(1,0,0),0.02*2,0.03*2));
            Seed_Particles(torus,0,0,density,particles_per_cell);
            int m=particles.number;
            GRID<TV> torus_grid(TV_INT(2,3,2),RANGE<TV>(TV(0,1.5,0),TV(.6,2.5,.6)));
            for(NODE_ITERATOR<TV> iterator(torus_grid);iterator.Valid();iterator.Next()){
                TV center=iterator.Location();
                T angle=random.Get_Uniform_Number((T)0,(T)pi*2);
                ROTATION<TV> rotation(angle,TV(0,1,0));
                for(int k=0;k<m;k++)
                    Add_Particle(center+rotation.Rotate(particles.X(k)),0,0,particles.mass(k),particles.volume(k));
                break;
            }
            for(int k=0;k<m;k++) particles.Add_To_Deletion_List(k);
            particles.Delete_Elements_On_Deletion_List(true);
            Add_Fixed_Corotated(150,0.4);
            Add_Gravity(TV(0,-9.8,0));
        } break;
        case 10:{ // torus into a box one by one
            grid.Initialize(TV_INT(resolution,resolution*2,resolution)+1,RANGE<TV>(TV(),TV(1,2,1)),true);
            T thickness=.2;
            Add_Collision_Object(RANGE<TV>(TV(-1,-1,-1),TV(2,.1,2)),COLLISION_TYPE::separate,.5);
            Add_Collision_Object(RANGE<TV>(TV(-1,-1,0.2-thickness),TV(2,.6,.2)),COLLISION_TYPE::separate,.5);
            Add_Collision_Object(RANGE<TV>(TV(-1,-1,0.8),TV(2,.6,0.8+thickness)),COLLISION_TYPE::separate,.5);
            Add_Collision_Object(RANGE<TV>(TV(0.2-thickness,-1,-1),TV(.2,.6,2)),COLLISION_TYPE::separate,.5);
            Add_Collision_Object(RANGE<TV>(TV(.8,-1,-1),TV(.8+thickness,.6,2)),COLLISION_TYPE::separate,.5);
            T density=5*scale_mass;
            ANALYTIC_IMPLICIT_OBJECT<TORUS<T> > torus(TORUS<T>(TV(),TV(1,0,0),0.02*2,0.03*2));
            Seed_Particles(torus,0,0,density,particles_per_cell);
            for(int i=0;i<particles.number;i++) particles.valid(i)=false;
            case10_m=particles.number;
            Add_Gravity(TV(0,-9.8,0));
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
    // static int i=0;
    switch(test_number)
    {
        case 10:{
            if(frame%8==0 && frame<200){
                TV center=random.Get_Uniform_Vector(TV(.4,1,.4),TV(.6,1,.6));
                T angle=random.Get_Uniform_Number((T)0,(T)pi*2);
                ROTATION<TV> rotation(angle,TV(0,1,0));
                int old_m=particles.number;
                for(int k=0;k<case10_m;k++) Add_Particle(center+rotation.Rotate(particles.X(k)),0,0,particles.mass(k),particles.volume(k));
                ARRAY<int> mpm_particles;
                for(int k=old_m;k<old_m+case10_m;k++) mpm_particles.Append(k);
                Add_Fixed_Corotated(150,0.3,&mpm_particles);}
        } break;
    }
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

//#####################################################################
// Function Initialize_Implicit_Surface
//
// This was copied from DEFORMABLES_STANDARD_TESTS.cpp
// TODO: put this function somewhere more convenient (maybe as a constructor of LEVELSET_IMPLICIT_OBJECT?)
//#####################################################################
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* STANDARD_TESTS<VECTOR<T,3> >::
Initialize_Implicit_Surface(TRIANGULATED_SURFACE<T>& surface,int max_res)
{
    typedef VECTOR<int,TV::m> TV_INT;
    LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >& undeformed_levelset=*LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >::Create();
    surface.Update_Bounding_Box();
    RANGE<VECTOR<T,3> > box=*surface.bounding_box;
    GRID<VECTOR<T,3> >& ls_grid=undeformed_levelset.levelset.grid;
    ARRAY<T,TV_INT>& phi=undeformed_levelset.levelset.phi;
    ls_grid=GRID<VECTOR<T,3> >::Create_Grid_Given_Cell_Size(box,box.Edge_Lengths().Max()/max_res,false,5);
    phi.Resize(ls_grid.Domain_Indices());
    LEVELSET_MAKER_UNIFORM<VECTOR<T,3> >::Compute_Level_Set(surface,ls_grid,0,phi);
    undeformed_levelset.Update_Box();
    return &undeformed_levelset;
}

template class STANDARD_TESTS<VECTOR<float,3> >;
template class STANDARD_TESTS<VECTOR<double,3> >;
}
