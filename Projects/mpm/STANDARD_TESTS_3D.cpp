//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/KD_TREE.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TORUS.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/OPENSUBDIV_SURFACE_CURVATURE_FORCE.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE_3D.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Forces/VOLUME_PRESERVING_OB_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include "STANDARD_TESTS_3D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type_input,parse_args),Nsurface(0),
    foo_int1(0),foo_T1(0),foo_T2(0),foo_T3(0),foo_T4(0),foo_surface1(0),foo_surface2(0),
    foo_levelset1(0),foo_cylinder(0) 
{
    parse_args.Add("-fooint1",&foo_int1,"int1","a interger");
    parse_args.Add("-fooT1",&foo_T1,"T1","a scalar");
    parse_args.Add("-fooT2",&foo_T2,"T2","a scalar");
    parse_args.Add("-fooT3",&foo_T3,"T3","a scalar");
    parse_args.Add("-fooT4",&foo_T4,"T4","a scalar");
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
~STANDARD_TESTS()
{
    if(foo_surface1) delete foo_surface1;
    if(foo_surface2) delete foo_surface2;
    if(foo_levelset1) delete foo_levelset1;
    if(foo_cylinder) delete foo_cylinder;
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
            grid.Initialize(TV_INT(resolution*3,resolution*2,resolution*3)+1,RANGE<TV>(TV(-1,0,-1),TV(2,2,2)),true);
            RANGE<TV> ym(TV(0.2,0.0,0.2),TV(0.8,0.2,0.8));
            RANGE<TV> xm(TV(0.2,0.2,0.2),TV(0.3,0.5,0.8));
            RANGE<TV> xM(TV(0.7,0.2,0.2),TV(0.8,0.5,0.8));
            RANGE<TV> zm(TV(0.3,0.2,0.2),TV(0.7,0.5,0.3));
            RANGE<TV> zM(TV(0.3,0.2,0.7),TV(0.7,0.5,0.8));
            Add_Penalty_Collision_Object(ym);
            Add_Penalty_Collision_Object(xm);
            Add_Penalty_Collision_Object(xM);
            Add_Penalty_Collision_Object(zm);
            Add_Penalty_Collision_Object(zM);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,true);
            T density=5*scale_mass;
            ANALYTIC_IMPLICIT_OBJECT<TORUS<T> > torus(TORUS<T>(TV(),TV(1,0,0),0.02*1.5,0.04*1.5));
            Seed_Particles(torus,0,0,density,particles_per_cell);
            for(int i=0;i<particles.number;i++) particles.valid(i)=false;
            foo_int1=particles.number;
            LOG::cout<<"one torus particle #: "<<foo_int1<<std::endl;
            Add_Gravity(TV(0,-9.8,0));
        } break;
        case 11:{ // skew impact of two elastic spheres with initial angular velocity
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(),TV(30,30,30)),true);
            T density=5*scale_mass;
            SPHERE<TV> sphere1(TV(10,13,15),2);
            VECTOR<T,3> angular_velocity1(TV(0,0,foo_T1));
            Seed_Particles_Helper(sphere1,[=](const TV& X){return angular_velocity1.Cross(X-sphere1.center)+TV(0.75,0,0);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity1);},density,particles_per_cell);
            SPHERE<TV> sphere2(TV(20,15,15),2);
            VECTOR<T,3> angular_velocity2(TV(0,0,foo_T2));
            Seed_Particles_Helper(sphere2,[=](const TV& X){return angular_velocity2.Cross(X-sphere2.center)+TV(-0.75,0,0);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity2);},density,particles_per_cell);
            Add_Neo_Hookean(31.685*scale_E,0.44022); //solve({E/(2*(1+r))=11,E*r/((1+r)*(1-2*r))=81},{E,r});
        } break;
        case 12:{ // surface tension test: fixed topology circle shell
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-1.5,-1.5,-1.5),TV(1.5,1.5,1.5)),true);
            T density=1*scale_mass;
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(.0f),data_directory+"/Rigid_Bodies/sphere.tri.gz",*surface);
            LOG::cout<<"Read mesh "<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh "<<surface->particles.number<<std::endl;
            TRIANGULATED_SURFACE<T>& new_sc=Seed_Lagrangian_Particles(*surface,[=](const TV& X){return TV();},0,density,true);
            SURFACE_TENSION_FORCE_3D<TV>* stf=new SURFACE_TENSION_FORCE_3D<TV>(new_sc,(T).1);
            Add_Force(*stf);
            this->deformable_body_collection.Test_Forces(0);
            Add_Neo_Hookean(31.685*scale_E,0.44022); //solve({E/(2*(1+r))=11,E*r/((1+r)*(1-2*r))=81},{E,r});
        } break;
        case 13:{ // drip drop
            grid.Initialize(TV_INT(1,3,1)*resolution,RANGE<TV>(TV(0,-2,0),TV(1,1,1)),true);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,true);
            Add_Gravity(TV(0,-1,0));
            T density=2*scale_mass;
            RANGE<TV> box(TV(.4,.5,.4),TV(.6,.85,.6));
            Seed_Particles_Helper(box,[=](const TV& X){return TV();},0,density,particles_per_cell);
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            bool no_mu=true;
            Add_Fixed_Corotated(scale_E*20,0.3,&mpm_particles,no_mu);
            Add_Force(*new MPM_VISCOSITY<TV>(force_helper,gather_scatter,0,0.000001));
            use_surface_tension=true;
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness,penalty_damping_stiffness);
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i)(1)>0.8){TV x=particles.X(i);pinning_force->Add_Target(i,[=](T time){return x;});}
            Add_Force(*pinning_force);
        } break;
        case 14:{ // drop an oldroyd-b to a ground SCA energy
            grid.Initialize(TV_INT(resolution*2,resolution,resolution*2),RANGE<TV>(TV(-1,0,-1),TV(1,1,1)),true);
            RANGE<TV> ym(TV(-5,0,-5),TV(5,.1,5));
            Add_Penalty_Collision_Object(ym);
            SPHERE<TV> sphere(TV(0,.5,0),.2);
            T density=2*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)10;
            particles.Store_S(use_oldroyd);            
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV();},0,density,particles_per_cell);
            particles.F.Fill(MATRIX<T,3>()+1);particles.S.Fill(SYMMETRIC_MATRIX<T,3>()+sqr(1));
            VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV> *neo=new VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>;
            neo->mu=38.462; // E=100, nu=0.3
            neo->lambda=57.692;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
            Add_Gravity(TV(0,-9.8,0));
        } break;
        case 15:{ // rotating cylinder oldroyd-b SCA energy
            //NEWTONIAN ./mpm -3d 15 -affine -max_dt 5e-4 -resolution 15 -fooint1 2 -fooT1 0 -fooT2 100 -fooT3 1e30 -fooT4 2e-4 -scale_mass 20 -last_frame 30
            grid.Initialize(TV_INT(resolution*2,resolution*3,resolution*2),RANGE<TV>(TV(-0.06,-0.06,-0.06),TV(0.06,0.12,0.06)),true);
            LOG::cout<<"GRID DX: "<<grid.dX<<std::endl;
            Add_Gravity(TV(0,-9.8,0));
            // Container glass
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/glass.tri.gz",*surface);
            LOG::cout<<"Read mesh elements "<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh particles "<<surface->particles.number<<std::endl;
            surface->mesh.Initialize_Adjacent_Elements();    
            surface->mesh.Initialize_Neighbor_Nodes();
            surface->mesh.Initialize_Incident_Elements();
            surface->Update_Bounding_Box();
            surface->Initialize_Hierarchy();
            surface->Update_Triangle_List();
            LOG::cout<<"Building levelset for the collision object..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Initialize_Implicit_Surface(*surface,100);
            LOG::cout<<"...done!"<<std::endl;
            Add_Penalty_Collision_Object(levelset);
            // Rotating cylinder
            if(0){
                foo_surface1=TRIANGULATED_SURFACE<T>::Create();
                FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/cylinder.tri.gz",*foo_surface1);
                foo_surface2=TRIANGULATED_SURFACE<T>::Create();
                FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/cylinder.tri.gz",*foo_surface2);
                LOG::cout<<"Read mesh elements "<<foo_surface1->mesh.elements.m<<std::endl;
                LOG::cout<<"Read mesh particles "<<foo_surface1->particles.number<<std::endl;
                foo_surface1->mesh.Initialize_Adjacent_Elements();    
                foo_surface1->mesh.Initialize_Neighbor_Nodes();
                foo_surface1->mesh.Initialize_Incident_Elements();
                foo_surface1->Update_Bounding_Box();
                foo_surface1->Initialize_Hierarchy();
                foo_surface1->Update_Triangle_List();
                LOG::cout<<"Building levelset for the cylinder stir..."<<std::endl;
                foo_levelset1=Initialize_Implicit_Surface(*foo_surface1,100);
                LOG::cout<<"...done!"<<std::endl;
                Add_Penalty_Collision_Object(foo_levelset1);}
            if(1){
                foo_cylinder=new CYLINDER<T>(TV(-0.025,-0.034,0),TV(0.025,-0.034,0),0.007);
                Add_Penalty_Collision_Object(*foo_cylinder);}
            // Polyethylene glycol
            SPHERE<TV> seeder1(TV(0,0.015,0),.03);
            T radius=0.03*1.3;TV P1(0,0.014-0.035-0.018,0);TV P2(0,0.014+0.05,0); CYLINDER<T> seeder2(P1,P2,radius);
            T density=1*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)1./foo_T3;
            particles.Store_S(use_oldroyd);            
            if(foo_int1==1) Seed_Particles_Helper(seeder1,[=](const TV& X){return TV();},0,density,particles_per_cell);
            else if(foo_int1==2) Seed_Particles_Helper(seeder2,[=](const TV& X){return TV();},0,density,particles_per_cell);
            else PHYSBAM_FATAL_ERROR();
            for(int k=0;k<particles.number;k++){
                TV X=particles.X(k);
                if(sqr(X(0))+sqr(X(2))>=sqr(0.046) || levelset->Extended_Phi(X)<=0 
                    || (foo_levelset1 && foo_levelset1->Extended_Phi(X)<=0) 
                    || (foo_cylinder && foo_cylinder->Lazy_Inside(X)))
                    particles.Add_To_Deletion_List(k);}
            particles.Delete_Elements_On_Deletion_List();
            particles.F.Fill(MATRIX<T,3>()+1);particles.S.Fill(SYMMETRIC_MATRIX<T,3>()+sqr(1));
            VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV> *neo=new VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>;
            neo->mu=foo_T1;
            neo->lambda=foo_T2;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
            Add_Force(*new MPM_VISCOSITY<TV>(force_helper,gather_scatter,0,foo_T4));
            LOG::cout<<"Polyethylene glycol added. mu="<<neo->mu<<", lambda="<<neo->lambda<<", Weissenbergi="<<foo_T3<<std::endl;
            LOG::cout<<"Particle count: "<<particles.number<<std::endl;
        } break;
        case 16:{ // rotating cylinder oldroyd-b SCA energy with pinned particles as the cylinder
            // NEWTONIAN ./mpm -3d 16 -affine -max_dt 5e-4 -resolution 15 -fooint1 2 -fooT1 0 -fooT2 100 -fooT3 1e30 -fooT4 1e-4 -scale_mass 20 -last_frame 200 -penalty_stiffness 10 -penalty_damping 0.001 -framerate 72 -cfl .7 -threads 16
            grid.Initialize(TV_INT(resolution*2,resolution*3,resolution*2),RANGE<TV>(TV(-0.06,-0.06,-0.06),TV(0.06,0.12,0.06)),true);
            LOG::cout<<"GRID DX: "<<grid.dX<<std::endl;
            Add_Gravity(TV(0,-9.8,0));
            // Container glass
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/glass_tall.tri.gz",*surface);
            LOG::cout<<"Read mesh elements "<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh particles "<<surface->particles.number<<std::endl;
            surface->mesh.Initialize_Adjacent_Elements();    
            surface->mesh.Initialize_Neighbor_Nodes();
            surface->mesh.Initialize_Incident_Elements();
            surface->Update_Bounding_Box();
            surface->Initialize_Hierarchy();
            surface->Update_Triangle_List();
            LOG::cout<<"Building levelset for the collision object..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Initialize_Implicit_Surface(*surface,100);
            LOG::cout<<"...done!"<<std::endl;
            Add_Penalty_Collision_Object(levelset);
            // Polyethylene glycol
            SPHERE<TV> seeder1(TV(0,0.015,0),.03);
            T radius=0.03*1.3;TV P1(0,0.014-0.035-0.018,0);TV P2(0,0.014+0.04,0); CYLINDER<T> seeder2(P1,P2,radius);
            T density=1*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)1./foo_T3;
            particles.Store_S(use_oldroyd);            
            if(foo_int1==1) Seed_Particles_Helper(seeder1,[=](const TV& X){return TV();},0,density,particles_per_cell);
            else if(foo_int1==2) Seed_Particles_Helper(seeder2,[=](const TV& X){return TV();},0,density,particles_per_cell);
            else PHYSBAM_FATAL_ERROR();
            for(int k=0;k<particles.number;k++){
                TV X=particles.X(k);
                if(sqr(X(0))+sqr(X(2))>=sqr(0.046) || levelset->Extended_Phi(X)<=0)
                    particles.Add_To_Deletion_List(k);}
            particles.Delete_Elements_On_Deletion_List();
            particles.F.Fill(MATRIX<T,3>()+1);particles.S.Fill(SYMMETRIC_MATRIX<T,3>()+sqr(1));
            VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV> *neo=new VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>;
            neo->mu=foo_T1;
            neo->lambda=foo_T2;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
            Add_Force(*new MPM_VISCOSITY<TV>(force_helper,gather_scatter,0,foo_T4));
            LOG::cout<<"Polyethylene glycol added. mu="<<neo->mu<<", lambda="<<neo->lambda<<", Weissenbergi="<<foo_T3<<std::endl;
            LOG::cout<<"Particle count: "<<particles.number<<std::endl;
            // Give particles in cylinder pinning forces
            foo_cylinder=new CYLINDER<T>(TV(-0.025,-0.028,0),TV(0.025,-0.028,0),0.007);
            VECTOR<T,3> angular_velocity(0,(T)31.41592653*scale_speed,0);
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness,
                penalty_damping_stiffness);
            for(int i=0;i<particles.X.m;i++){
                if(foo_cylinder->Lazy_Inside(particles.X(i))){
                    particles.mass(i)*=(T)1.1;
                    TV dx=particles.X(i)-foo_cylinder->Bounding_Box().Center();
                    pinning_force->Add_Target(i,
                        [=](T time){
                            ROTATION<TV> rot=ROTATION<TV>::From_Rotation_Vector(angular_velocity*time);
                            return rot.Rotate(dx)+foo_cylinder->Bounding_Box().Center();});}}
            Add_Force(*pinning_force);
        } break;
        case 17:{ // sand box drop
            use_plasticity=true;
            use_variable_coefficients=true;
            particles.Store_Fp(true);
            particles.Store_Mu(true);
            particles.Store_Lambda(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            Add_Collision_Object(RANGE<TV>(TV(-0.5,-1,-0.5),TV(1.5,.1,1.5)),COLLISION_TYPE::separate,10);

            T density=(T)1281*scale_mass;
            T E=5000*scale_E,nu=.4;
            this->mu0=E/(2*(1+nu));
            this->lambda0=E*nu/((1+nu)*(1-2*nu));
            if(theta_c==0) theta_c=0.01;
            if(theta_s==0) theta_s=.00001;
            if(hardening_factor==0) hardening_factor=80;
            if(max_hardening) max_hardening=5;
            Add_Fixed_Corotated(E,nu);
            RANGE<TV> box(TV(.4,.15,.4),TV(.6,.35,.6));
            Seed_Particles_Helper(box,0,0,density,particles_per_cell);
            for(int p=0;p<particles.number;++p){
                particles.mu(p)=this->mu0;
                particles.lambda(p)=this->lambda0;}
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
                for(int k=0;k<foo_int1;k++) Add_Particle(center+rotation.Rotate(particles.X(k)),0,0,particles.mass(k),particles.volume(k));
                ARRAY<int> mpm_particles;
                for(int k=old_m;k<old_m+foo_int1;k++) mpm_particles.Append(k);
                Add_Fixed_Corotated(150*scale_E,0.3,&mpm_particles);}
        } break;
        case 15:{
            if(foo_levelset1){
                PHYSBAM_ASSERT(!foo_cylinder);
                delete lagrangian_forces(lagrangian_forces.m-1);
                lagrangian_forces.Remove_End();
                PHYSBAM_ASSERT(this->deformable_body_collection.structures.m==0);
                LOG::cout<<"Building levelset for the cylinder stir..."<<std::endl;
                if(foo_levelset1) delete foo_levelset1;
                LOG::cout<<"Deleted old levelset."<<std::endl;
                ROTATION<TV> rotator((T)0.02*frame,TV(0,1,0));
                for(int k=0;k<foo_surface1->particles.number;k++)
                foo_surface1->particles.X(k)=rotator.Rotate(foo_surface2->particles.X(k));
                foo_surface1->Update_Bounding_Box();
                foo_surface1->Initialize_Hierarchy();
                foo_levelset1=Initialize_Implicit_Surface(*foo_surface1,50);
                LOG::cout<<"...done!"<<std::endl;
                Add_Penalty_Collision_Object(foo_levelset1);
                Dump_Levelset(foo_levelset1->levelset.grid,*foo_levelset1,VECTOR<T,3>(0,1,0));}
            else if(foo_cylinder){
                LOG::cout<<"Dumping cylinder stirer level set.."<<std::endl;
                GRID<TV> ghost_grid(TV_INT(50,50,10),foo_cylinder->Bounding_Box(),true);
                Dump_Levelset(ghost_grid,ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(*foo_cylinder),VECTOR<T,3>(1,1,0));
                LOG::cout<<"...done!"<<std::endl;}
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
    switch(test_number)
    {
        case 14:{
            if(time>=10/24.0){
                lagrangian_forces.Delete_Pointers_And_Clean_Memory();
                this->deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();
                RANGE<TV> ym(TV(-5,0,-5),TV(5,.1+(time-10/24.0)*.5,5));
                Add_Penalty_Collision_Object(ym);
                Add_Gravity(TV(0,-9.8,0));}
        } break;
        case 15:{
            if(foo_cylinder){
                PHYSBAM_ASSERT(!foo_levelset1);
                delete lagrangian_forces(lagrangian_forces.m-1);
                lagrangian_forces.Remove_End();
                PHYSBAM_ASSERT(this->deformable_body_collection.structures.m==0);
                LOG::cout<<"Adding new analytic cylinder stirer..."<<std::endl;
                ROTATION<TV> rotator((T)3.1415*10*time,TV(0,1,0));
                foo_cylinder->Set_Endpoints(rotator.Rotate(TV(-0.025,-0.034,0)),rotator.Rotate(TV(0.025,-0.034,0)));
                Add_Penalty_Collision_Object(*foo_cylinder);
                LOG::cout<<"...done!"<<std::endl;}
        } break;
        default: break;
    }

    if(use_surface_tension){

        bool use_bruteforce=false;
        bool use_kdtree=true;

        // Remove old surface particles
        int N_non_surface=particles.number-Nsurface;
        for(int k=N_non_surface;k<particles.number;k++){
            particles.Add_To_Deletion_List(k);
            int m=steal(k-N_non_surface);
            TV old_momentum=particles.mass(m)*particles.V(m)+particles.mass(k)*particles.V(k);
            particles.mass(m)+=particles.mass(k);
            particles.volume(m)+=particles.volume(k);
            particles.V(m)=old_momentum/particles.mass(m);
            if(use_affine) particles.B(m)=(particles.B(m)+particles.B(k))*0.5;}
        LOG::cout<<"deleting "<<Nsurface<<" particles..."<<std::endl;
        particles.Delete_Elements_On_Deletion_List();
        lagrangian_forces.Delete_Pointers_And_Clean_Memory();
        this->deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();

        // Dirty hack of rasterizing mass
        this->simulated_particles.Remove_All();
        for(int p=0;p<this->particles.number;p++)
            if(this->particles.valid(p))
                this->simulated_particles.Append(p);
        this->particle_is_simulated.Remove_All();
        this->particle_is_simulated.Resize(this->particles.X.m);
        this->particle_is_simulated.Subset(this->simulated_particles).Fill(true);
        this->weights->Update(this->particles.X);
        this->gather_scatter.Prepare_Scatter(this->particles);
        MPM_PARTICLES<TV>& my_particles=this->particles;
#pragma omp parallel for
        for(int i=0;i<this->mass.array.m;i++)
            this->mass.array(i)=0;
        this->gather_scatter.template Scatter<int>(
            [this,&my_particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
            {
                T w=it.Weight();
                TV_INT index=it.Index();
                this->mass(index)+=w*my_particles.mass(p);
            },false);

        // Marching cube
        TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
        MARCHING_CUBES<TV>::Create_Surface(*surface,grid,mass,particles.mass(0)*1.5);

        // Improve surface quality
        T min_edge_length=FLT_MAX;
        for(int k=0;k<surface->mesh.elements.m;k++){
            int a=surface->mesh.elements(k)(0),b=surface->mesh.elements(k)(1);
            TV A=surface->particles.X(a),B=surface->particles.X(b);
            T l2=(A-B).Magnitude_Squared();
            if(l2<min_edge_length) min_edge_length=l2;}
        min_edge_length=sqrt(min_edge_length);
        LOG::cout<<"Marching cube min edge length: "<<min_edge_length<<std::endl;

        // Seed surface particles
        int Nold=particles.number;
        Nsurface=surface->particles.number;
        LOG::cout<<"adding "<<Nsurface<<" particles..."<<std::endl;
        TRIANGULATED_SURFACE<T>& new_sc=Seed_Lagrangian_Particles(*surface,0,0,(T)0.001,true);
        ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
        for(int i=Nold;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(0,1,0);

        // Build K-d tree for non-surface particles
        LOG::cout<<"building kdtree..."<<std::endl;
        KD_TREE<TV> kdtree;
        ARRAY<TV> nodes(Nold);
        if(use_kdtree){
            for(int p=0;p<Nold;p++) nodes(p)=particles.X(p);
            kdtree.Create_Left_Balanced_KD_Tree(nodes);}

        // Assign physical quantities
        steal.Clean_Memory();
        for(int k=Nold;k<particles.number;k++){
            T dist2=FLT_MAX;
            int m=-1;

            // Find closest interior particle using brute force
            if(use_bruteforce){
                for(int q=0;q<Nold;q++){
                    T dd=(particles.X(q)-particles.X(k)).Magnitude_Squared();
                    if(dd<dist2){dist2=dd;m=q;}}}

            // Find closest interior particle using kdtree
            int number_of_points_in_estimate=1;
            ARRAY<int> points_found(number_of_points_in_estimate);
            ARRAY<T> squared_distance_to_points_found(number_of_points_in_estimate);
            if(use_kdtree){
                int number_of_points_found;T max_squared_distance_to_points_found;
                kdtree.Locate_Nearest_Neighbors(particles.X(k),FLT_MAX,points_found,
                    squared_distance_to_points_found,number_of_points_found,max_squared_distance_to_points_found,nodes);
                PHYSBAM_ASSERT(number_of_points_found==number_of_points_in_estimate);}

            // Debug kdtree
            if(use_bruteforce && use_kdtree && m!=points_found(0)){
                LOG::cout<<"Disagree!"<<std::endl; PHYSBAM_FATAL_ERROR();}

            if(use_kdtree) m=points_found(0);

            T split_mass=particles.mass(m)*.5;
            T split_volume=particles.volume(m)*.5;
            TV com=particles.X(m);
            particles.mass(m)=split_mass;
            particles.volume(m)=split_volume;
            particles.mass(k)=split_mass;
            particles.volume(k)=split_volume;
            particles.V(k)=particles.V(m);
            if(use_affine) particles.B(k)=particles.B(m);
            steal.Append(m);}

        // Add surface tension force
        SURFACE_TENSION_FORCE_3D<TV>* stf=new SURFACE_TENSION_FORCE_3D<TV>(new_sc,(T)2e-2);
        Add_Force(*stf);
    }
    
    if(test_number==13){
        Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,true);
        Add_Gravity(TV(0,-1,0));
        PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness,penalty_damping_stiffness);
        for(int i=0;i<particles.X.m;i++)
            if(particles.X(i)(1)>0.8){TV x=particles.X(i);pinning_force->Add_Target(i,[=](T time){return x;});}
        Add_Force(*pinning_force);
    }

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
    phi.Resize(ls_grid.Domain_Indices(3));
    LEVELSET_MAKER_UNIFORM<VECTOR<T,3> >::Compute_Level_Set(surface,ls_grid,3,phi);
    undeformed_levelset.Update_Box();
    return &undeformed_levelset;
}

template class STANDARD_TESTS<VECTOR<float,3> >;
template class STANDARD_TESTS<VECTOR<double,3> >;
}
