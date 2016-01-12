//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/KD_TREE.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Matrices/FRAME.h>
#include <Geometry/Basic_Geometry/CONE.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/HOURGLASS.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TORUS.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INTERSECTION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INVERT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UNION.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
#include <Deformables/Constitutive_Models/MPM_DRUCKER_PRAGER.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/OPENSUBDIV_SURFACE_CURVATURE_FORCE.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE_3D.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Forces/VOLUME_PRESERVING_OB_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include "POUR_SOURCE.h"
#include "STANDARD_TESTS_3D.h"
#include <fstream>
namespace PhysBAM{
//#####################################################################
// Function Initialize_Implicit_Surface
//
// This was copied from DEFORMABLES_STANDARD_TESTS.cpp
// TODO: put this function somewhere more convenient (maybe as a constructor of LEVELSET_IMPLICIT_OBJECT?)
//#####################################################################
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >*
Initialize_Implicit_Surface(TRIANGULATED_SURFACE<T>& surface,int max_res)
{
    typedef VECTOR<int,3> TV_INT;
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
//#####################################################################
// Function Levelset_From_File
//
// TODO: put this function where it belongs
//#####################################################################
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* 
Levelset_From_File(const std::string& filename,int max_resolution=200)
{
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),filename,*surface);
    LOG::printf("Read mesh: %d triangle, %d particles\n",surface->mesh.elements.m,surface->particles.number);
    surface->mesh.Initialize_Adjacent_Elements();
    surface->mesh.Initialize_Neighbor_Nodes();
    surface->mesh.Initialize_Incident_Elements();
    surface->Update_Bounding_Box();
    surface->Initialize_Hierarchy();
    surface->Update_Triangle_List();
    return Initialize_Implicit_Surface(*surface,max_resolution);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type_input,parse_args),Nsurface(0),
    foo_int1(0),foo_T1(0),foo_T2(0),foo_T3(0),foo_T4(0),foo_T5(0),
    use_foo_T1(false),use_foo_T2(false),use_foo_T3(false),use_foo_T4(false),use_foo_T5(false),
    foo_surface1(0),foo_surface2(0),foo_levelset1(0),foo_cylinder(0)
{
    parse_args.Add("-fooint1",&foo_int1,"int1","a interger");
    parse_args.Add("-fooT1",&foo_T1,&use_foo_T1,"T1","a scalar");
    parse_args.Add("-fooT2",&foo_T2,&use_foo_T2,"T2","a scalar");
    parse_args.Add("-fooT3",&foo_T3,&use_foo_T3,"T3","a scalar");
    parse_args.Add("-fooT4",&foo_T4,&use_foo_T4,"T4","a scalar");
    parse_args.Add("-fooT5",&foo_T5,&use_foo_T5,"T5","a scalar");
    parse_args.Parse();
    if(!this->override_output_directory) output_directory=LOG::sprintf("Test_3d_%i",test_number);
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
    if(destroy) destroy();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Write_Output_Files(const int frame)
{
    if(write_output_files) write_output_files(frame);
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Read_Output_Files(const int frame)
{
    if(read_output_files) read_output_files(frame);
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
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            VECTOR<T,3> angular_velocity(TV(0.4,0,0)/s);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity);}
                ,density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
        } break;
        case 2:{ // Oscillating sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,3>()+1.5);
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
        } break;
        case 3:{ // Freefall sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 4:{ // subdivision surface - 40x40 strip
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-3,-3,-3),TV(3,3,3))*m,true);
            
            
            std::string filename=data_directory+"/OpenSubdiv_Surfaces/strip_40.dat.gz";
            OPENSUBDIV_SURFACE<TV>* surf=OPENSUBDIV_SURFACE<TV>::Create();
            T thickness=1e-3*m;
            surf->Initialize(filename,thickness);

            T density=1000*unit_rho*scale_mass;
//            surf->Set_Masses(density,thickness);

            T c1=3.5e6; // TODO: units!!!
            T c2=1.3e6;
            T stiffness_multiplier=scale_E;
            T thickness_multiplier=1;

            OPENSUBDIV_SURFACE<TV>& new_surf=Seed_Lagrangian_Particles(*surf,0,[=](const TV&){return MATRIX<T,3>();},density,false);
            MOONEY_RIVLIN_CURVATURE<T> model=MOONEY_RIVLIN_CURVATURE<T>(c1*stiffness_multiplier,c2*stiffness_multiplier,thickness*thickness_multiplier);
            Add_Force(*new OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,3>(particles,new_surf,model));

            Add_Walls(8,COLLISION_TYPE::separate,.3,.1*m,false);

            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 5:{ // subdivision surface - duck
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-3,-3,-3),TV(3,3,3))*m,true);
            
            std::string filename=data_directory+"/OpenSubdiv_Surfaces/duck_1073f.dat.gz";
            OPENSUBDIV_SURFACE<TV>* surf=OPENSUBDIV_SURFACE<TV>::Create();
            T thickness=1e-3;
            surf->Initialize(filename,thickness);
            
            T density=1000*unit_rho*scale_mass;

            T c1=3.5e6;
            T c2=1.3e6;
            T stiffness_multiplier=unit_p*scale_E; // TODO: units!!!
            T thickness_multiplier=1;

            OPENSUBDIV_SURFACE<TV>& new_surf=Seed_Lagrangian_Particles(*surf,0,[=](const TV&){return MATRIX<T,3>();},density,false);
            MOONEY_RIVLIN_CURVATURE<T> model=MOONEY_RIVLIN_CURVATURE<T>(c1*stiffness_multiplier,c2*stiffness_multiplier,thickness*thickness_multiplier);
            Add_Force(*new OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,3>(particles,new_surf,model));

//            for(int p=0;p<particles.X.m;p++) // squish the duck.
//                particles.X(p)(0)*=0.5;

            Add_Walls(8,COLLISION_TYPE::separate,.3,.1*m,false);

            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 6:{ // subdivision surface - drop several ducks.
            int num_duckies=3;
            grid.Initialize(resolution*TV_INT(2,num_duckies,2),RANGE<TV>(TV(-6,-3,-6),TV(6,3+6*(num_duckies-1),6))*m,true);
            
            std::string filename=data_directory+"/OpenSubdiv_Surfaces/duck_1073f.dat.gz";

            T c1=3.5e6;
            T c2=1.3e6;
            T stiffness_multiplier=unit_p*scale_E; // TODO: units!!!
            T thickness_multiplier=1;
            T density=1000;
            T thickness=1e-3;

            MOONEY_RIVLIN_CURVATURE<T> model=MOONEY_RIVLIN_CURVATURE<T>(c1*stiffness_multiplier,c2*stiffness_multiplier,thickness*thickness_multiplier);

            OPENSUBDIV_SURFACE<TV>* surf=OPENSUBDIV_SURFACE<TV>::Create();
            surf->Initialize(filename,thickness);
            for(int i=0;i<num_duckies;i++){
                for(int p=0;p<particles.X.m;p++)
                    particles.X(p)(1)+=6;
                OPENSUBDIV_SURFACE<TV>& new_surf=Seed_Lagrangian_Particles(*surf,0,[=](const TV&){return MATRIX<T,3>();},density,false,false);
                Add_Force(*new OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,3>(particles,new_surf,model));
            }

            Add_Walls(8,COLLISION_TYPE::separate,.3,.1*m,false);

            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
            delete surf;
        } break;
        case 7:{ // skew impact of two elastic spheres
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(),TV(30,30,30))*m,true);
            T density=5*unit_rho*scale_mass;
            SPHERE<TV> sphere1(TV(10,13,15)*m,2*m);
            Seed_Particles(sphere1,[=](const TV& X){return TV(0.75,0,0)*(m/s);},[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            SPHERE<TV> sphere2(TV(20,15,15)*m,2*m);
            Seed_Particles(sphere2,[=](const TV& X){return TV(-0.75,0,0)*(m/s);},[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            Add_Neo_Hookean(31.685*unit_p*scale_E,0.44022); //solve({E/(2*(1+r))=11,E*r/((1+r)*(1-2*r))=81},{E,r});
        } break;
        case 8:{ // torus into a box
            // TODO: fix crash ./mpm -3d 8 -affine -last_frame 500 -midpoint -resolution 40 -newton_tolerance 1e-3 
            grid.Initialize(TV_INT(resolution,resolution*2,resolution)+1,RANGE<TV>(TV(),TV(1,2,1))*m,true);

            // Add_Walls(8,COLLISION_TYPE::separate,.3,.1*m,false);

            T thickness=.1;
            Add_Collision_Object(RANGE<TV>(TV(.2,0,.2),TV(.8,.1,.8))*m,COLLISION_TYPE::separate,.3);
            Add_Collision_Object(RANGE<TV>(TV(.2-thickness,0,0.2-thickness),TV(.8+thickness,.4,.2))*m,COLLISION_TYPE::separate,.3);
            Add_Collision_Object(RANGE<TV>(TV(.2-thickness,0,0.8),TV(.8+thickness,.4,0.8+thickness))*m,COLLISION_TYPE::separate,.3);
            Add_Collision_Object(RANGE<TV>(TV(0.2-thickness,0,0.2),TV(.2,.4,.8))*m,COLLISION_TYPE::separate,.3);
            Add_Collision_Object(RANGE<TV>(TV(.8,0,0.2),TV(.8+thickness,.4,.8))*m,COLLISION_TYPE::separate,.3);

            T density=5*unit_rho*scale_mass;
            ANALYTIC_IMPLICIT_OBJECT<TORUS<T> > torus(TORUS<T>(TV(),TV(1,0,0)*m,0.02*2*m,0.03*2*m));
            Seed_Particles(torus,0,0,density,particles_per_cell);
            int n=particles.number;
            GRID<TV> torus_grid(TV_INT(2,3,2),RANGE<TV>(TV(.4,1,.4),TV(.6,1.5,.6))*m);
            for(NODE_ITERATOR<TV> iterator(torus_grid);iterator.Valid();iterator.Next()){
                TV center=iterator.Location();
                T angle=random.Get_Uniform_Number((T)0,(T)pi*2);
                ROTATION<TV> rotation(angle,TV(0,1,0));
                for(int k=0;k<n;k++)
                    Add_Particle(center+rotation.Rotate(particles.X(k)),0,0,particles.mass(k),particles.volume(k));
                break;
            }
            for(int k=0;k<n;k++) particles.Add_To_Deletion_List(k);
            particles.Delete_Elements_On_Deletion_List(true);
            Add_Fixed_Corotated(150*unit_p*scale_E,0.4);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;

        case 9:{ // torus into a bowl
            // TODO: fix crash in LEVELSET_MAKER_UNIFORM ./mpm -3d 9 -resolution 80
            grid.Initialize(TV_INT(resolution,resolution,resolution)+1,RANGE<TV>(TV(-2,-2,-2),TV(2,2,2))*m,true);

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

            T density=5*unit_rho*scale_mass;
            ANALYTIC_IMPLICIT_OBJECT<TORUS<T> > torus(TORUS<T>(TV(),TV(1,0,0)*m,0.02*2*m,0.03*2*m));
            Seed_Particles(torus,0,0,density,particles_per_cell);
            int n=particles.number;
            GRID<TV> torus_grid(TV_INT(2,3,2),RANGE<TV>(TV(0,1.5,0),TV(.6,2.5,.6))*m);
            for(NODE_ITERATOR<TV> iterator(torus_grid);iterator.Valid();iterator.Next()){
                TV center=iterator.Location();
                T angle=random.Get_Uniform_Number((T)0,(T)pi*2);
                ROTATION<TV> rotation(angle,TV(0,1,0));
                for(int k=0;k<n;k++)
                    Add_Particle(center+rotation.Rotate(particles.X(k)),0,0,particles.mass(k),particles.volume(k));
                break;
            }
            for(int k=0;k<n;k++) particles.Add_To_Deletion_List(k);
            particles.Delete_Elements_On_Deletion_List(true);
            Add_Fixed_Corotated(150*unit_p*scale_E,0.4);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 10:{ // torus into a box one by one
            grid.Initialize(TV_INT(resolution*3,resolution*2,resolution*3)+1,RANGE<TV>(TV(-1,0,-1),TV(2,2,2))*m,true);
            RANGE<TV> ym(TV(0.2,0.0,0.2)*m,TV(0.8,0.2,0.8)*m);
            RANGE<TV> xm(TV(0.2,0.2,0.2)*m,TV(0.3,0.5,0.8)*m);
            RANGE<TV> xM(TV(0.7,0.2,0.2)*m,TV(0.8,0.5,0.8)*m);
            RANGE<TV> zm(TV(0.3,0.2,0.2)*m,TV(0.7,0.5,0.3)*m);
            RANGE<TV> zM(TV(0.3,0.2,0.7)*m,TV(0.7,0.5,0.8)*m);
            Add_Penalty_Collision_Object(ym);
            Add_Penalty_Collision_Object(xm);
            Add_Penalty_Collision_Object(xM);
            Add_Penalty_Collision_Object(zm);
            Add_Penalty_Collision_Object(zM);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,true);
            T density=5*unit_rho*scale_mass;
            ANALYTIC_IMPLICIT_OBJECT<TORUS<T> > torus(TORUS<T>(TV(),TV(1,0,0)*m,0.02*1.5*m,0.04*1.5*m));
            Seed_Particles(torus,0,0,density,particles_per_cell);
            for(int i=0;i<particles.number;i++) particles.valid(i)=false;
            foo_int1=particles.number;
            LOG::cout<<"one torus particle #: "<<foo_int1<<std::endl;
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
            begin_frame=[this](int frame)
                {
                    if(frame%8==0 && frame<200){
                        TV center=random.Get_Uniform_Vector(TV(.4,1,.4),TV(.6,1,.6))*m;
                        T angle=random.Get_Uniform_Number((T)0,(T)pi*2);
                        ROTATION<TV> rotation(angle,TV(0,1,0));
                        int old_m=particles.number;
                        for(int k=0;k<foo_int1;k++) Add_Particle(center+rotation.Rotate(particles.X(k)),0,0,particles.mass(k),particles.volume(k));
                        ARRAY<int> mpm_particles;
                        for(int k=old_m;k<old_m+foo_int1;k++) mpm_particles.Append(k);
                        Add_Fixed_Corotated(150*unit_p*scale_E,0.3,&mpm_particles);}
                };
        } break;
        case 11:{ // skew impact of two elastic spheres with initial angular velocity
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(),TV(30,30,30))*m,true);
            T density=5*unit_rho*scale_mass;
            SPHERE<TV> sphere1(TV(10,13,15)*m,2*m);
            VECTOR<T,3> angular_velocity1(TV(0,0,foo_T1));
            Seed_Particles(sphere1,[=](const TV& X){return angular_velocity1.Cross(X-sphere1.center)+TV(0.75,0,0)*(m/s);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity1);},density,particles_per_cell);
            SPHERE<TV> sphere2(TV(20,15,15)*m,2*m);
            VECTOR<T,3> angular_velocity2(TV(0,0,foo_T2));
            Seed_Particles(sphere2,[=](const TV& X){return angular_velocity2.Cross(X-sphere2.center)+TV(-0.75,0,0)*(m/s);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity2);},density,particles_per_cell);
            Add_Neo_Hookean(31.685*unit_p*scale_E,0.44022); //solve({E/(2*(1+r))=11,E*r/((1+r)*(1-2*r))=81},{E,r});
        } break;
        case 12:{ // surface tension test: fixed topology circle shell
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-1.5,-1.5,-1.5),TV(1.5,1.5,1.5))*m,true);
            T density=1*unit_rho*scale_mass;
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(.0f),data_directory+"/Rigid_Bodies/sphere.tri.gz",*surface);
            LOG::cout<<"Read mesh "<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh "<<surface->particles.number<<std::endl;
            TRIANGULATED_SURFACE<T>& new_sc=Seed_Lagrangian_Particles(*surface,0,0,density,true);
            SURFACE_TENSION_FORCE_3D<TV>* stf=new SURFACE_TENSION_FORCE_3D<TV>(new_sc,(T).1);
            Add_Force(*stf);
            this->deformable_body_collection.Test_Forces(0);
            Add_Neo_Hookean(31.685*unit_p*scale_E,0.44022); //solve({E/(2*(1+r))=11,E*r/((1+r)*(1-2*r))=81},{E,r});
        } break;
        case 13:{ // drip drop
            grid.Initialize(TV_INT(1,3,1)*resolution,RANGE<TV>(TV(0,-2,0),TV(1,1,1))*m,true);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,true);
            Add_Gravity(m/(s*s)*TV(0,-1,0));
            T density=2*unit_rho*scale_mass;
            RANGE<TV> box(TV(.4,.5,.4)*m,TV(.6,.85,.6)*m);
            Seed_Particles(box,0,0,density,particles_per_cell);
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            bool no_mu=true;
            Add_Fixed_Corotated(unit_p*scale_E*20,0.3,&mpm_particles,no_mu);
            Add_Force(*new MPM_VISCOSITY<TV>(force_helper,gather_scatter,0,0.000001));
            use_surface_tension=true;
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness*kg/(m*s*s),penalty_damping_stiffness*kg/s);
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i)(1)>0.8*m){TV x=particles.X(i);pinning_force->Add_Target(i,[=](T time){return x;});}
            Add_Force(*pinning_force);
            begin_time_step=[this](T time)
                {
                    Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,true);
                    Add_Gravity(m/(s*s)*TV(0,-1,0));
                    PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness*kg/(m*s*s),penalty_damping_stiffness*kg/s);
                    for(int i=0;i<particles.X.m;i++)
                        if(particles.X(i)(1)>0.8*m){TV x=particles.X(i);pinning_force->Add_Target(i,[=](T time){return x;});}
                    Add_Force(*pinning_force);
                };
        } break;
        case 14:{ // drop an oldroyd-b to a ground SCA energy
            grid.Initialize(TV_INT(resolution*2,resolution,resolution*2),RANGE<TV>(TV(-1,0,-1),TV(1,1,1))*m,true);
            RANGE<TV> ym(TV(-5,0,-5)*m,TV(5,.1,5)*m);
            Add_Penalty_Collision_Object(ym);
            SPHERE<TV> sphere(TV(0,.5,0)*m,.2*m);
            T density=2*unit_rho*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)10;
            particles.Store_S(use_oldroyd);            
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            particles.F.Fill(MATRIX<T,3>()+1);particles.S.Fill(SYMMETRIC_MATRIX<T,3>()+sqr(1));
            VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV> *neo=new VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>;
            neo->mu=38.462*unit_p*scale_E; // E=100, nu=0.3
            neo->lambda=57.692*unit_p*scale_E;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
            begin_time_step=[this](T time)
                {
                    if(time>=10/24.0*s){
                        lagrangian_forces.Delete_Pointers_And_Clean_Memory();
                        this->deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();
                        RANGE<TV> ym(TV(-5,0,-5)*m,TV(5,.1+(time/s-10/24.0)*.5,5)*m);
                        Add_Penalty_Collision_Object(ym);
                        Add_Gravity(m/(s*s)*TV(0,-9.8,0));}
                };
        } break;
        case 15:{ // rotating cylinder oldroyd-b SCA energy
            //NEWTONIAN ./mpm -3d 15 -affine -max_dt 5e-4 -resolution 15 -fooint1 2 -fooT1 0 -fooT2 100 -fooT3 1e30 -fooT4 2e-4 -scale_mass 20 -last_frame 30
            grid.Initialize(TV_INT(resolution*2,resolution*3,resolution*2),RANGE<TV>(TV(-0.06,-0.06,-0.06),TV(0.06,0.12,0.06))*m,true);
            LOG::cout<<"GRID DX: "<<grid.dX<<std::endl;
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
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
                foo_cylinder=new CYLINDER<T>(TV(-0.025,-0.034,0)*m,TV(0.025,-0.034,0)*m,0.007*m);
                Add_Penalty_Collision_Object(*foo_cylinder);}
            // Polyethylene glycol
            SPHERE<TV> seeder1(TV(0,0.015,0)*m,.03*m);
            T radius=0.03*1.3*m;TV P1(0,(0.014-0.035-0.018)*m,0);TV P2(0,(0.014+0.05)*m,0); CYLINDER<T> seeder2(P1,P2,radius);
            T density=1*unit_rho*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)1./foo_T3;
            particles.Store_S(use_oldroyd);            
            if(foo_int1==1) Seed_Particles(seeder1,0,0,density,particles_per_cell);
            else if(foo_int1==2) Seed_Particles(seeder2,0,0,density,particles_per_cell);
            else PHYSBAM_FATAL_ERROR();
            for(int k=0;k<particles.number;k++){
                TV X=particles.X(k);
                if(sqr(X(0))+sqr(X(2))>=sqr(0.046*m) || levelset->Extended_Phi(X)<=0 
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
            begin_frame=[this](int frame)
                {
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
                };
            begin_time_step=[this](T time)
                {
                    if(foo_cylinder){
                        PHYSBAM_ASSERT(!foo_levelset1);
                        delete lagrangian_forces(lagrangian_forces.m-1);
                        lagrangian_forces.Remove_End();
                        PHYSBAM_ASSERT(this->deformable_body_collection.structures.m==0);
                        LOG::cout<<"Adding new analytic cylinder stirer..."<<std::endl;
                        ROTATION<TV> rotator((T)3.1415*10*time/s,TV(0,1,0));
                        foo_cylinder->Set_Endpoints(rotator.Rotate(TV(-0.025,-0.034,0)*m),rotator.Rotate(TV(0.025,-0.034,0)*m));
                        Add_Penalty_Collision_Object(*foo_cylinder);
                        LOG::cout<<"...done!"<<std::endl;}
                };
        } break;
        case 16:{ // rotating cylinder oldroyd-b SCA energy with pinned particles as the cylinder
            // NEWTONIAN ./mpm -3d 16 -affine -max_dt 5e-4 -resolution 15 -fooint1 2 -fooT1 0 -fooT2 100 -fooT3 1e30 -fooT4 1e-4 -scale_mass 20 -last_frame 200 -penalty_stiffness 10 -penalty_damping 0.001 -framerate 72 -cfl .7 -threads 16
            grid.Initialize(TV_INT(resolution*2,resolution*3,resolution*2),RANGE<TV>(TV(-0.06,-0.06,-0.06),TV(0.06,0.12,0.06))*m,true);
            LOG::cout<<"GRID DX: "<<grid.dX<<std::endl;
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
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
            SPHERE<TV> seeder1(TV(0,0.015,0)*m,.03*m);
            T radius=0.03*1.3*m;TV P1(0,(0.014-0.035-0.018)*m,0);TV P2(0,(0.014+0.04)*m,0); CYLINDER<T> seeder2(P1,P2,radius);
            T density=1*unit_rho*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)1./foo_T3;
            particles.Store_S(use_oldroyd);            
            if(foo_int1==1) Seed_Particles(seeder1,0,0,density,particles_per_cell);
            else if(foo_int1==2) Seed_Particles(seeder2,0,0,density,particles_per_cell);
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
            foo_cylinder=new CYLINDER<T>(TV(-0.025,-0.028,0)*m,TV(0.025,-0.028,0)*m,0.007*m);
            VECTOR<T,3> angular_velocity(0,(T)31.41592653*scale_speed/s,0);
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness*kg/(m*s*s),
                penalty_damping_stiffness*kg/s);
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
            particles.Store_Fp(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            Add_Collision_Object(RANGE<TV>(TV(-0.5,-1,-0.5),TV(1.5,.1,1.5))*m,COLLISION_TYPE::separate,10);

            T density=(T)1281*unit_rho*scale_mass;
            T E=5000*unit_p*scale_E,nu=.4;
            if(!use_theta_c) theta_c=0.01;
            if(!use_theta_s) theta_s=.00001;
            if(!use_hardening_factor) hardening_factor=80;
            if(!use_max_hardening) max_hardening=5;
            Add_Fixed_Corotated(E,nu);
            RANGE<TV> box(TV(.4,.15,.4)*m,TV(.6,.35,.6)*m);
            Seed_Particles(box,0,0,density,particles_per_cell);
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 18: // sand box drop, better paramaters, with Hencky, usage: mpm 18 -3d -resolution 100 -friction_angle 0.65 -cohesion 15
        case 19:{
            particles.Store_Fp(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            if(use_penalty_collisions)
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1,-0.5),TV(1.5,.1,1.5))*m,0.9);
            else
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1,-0.5),TV(1.5,.1,1.5))*m,COLLISION_TYPE::stick,0);

            T density=(T)1281*unit_rho*scale_mass; // source: Sand, dry http://www.engineeringtoolbox.com/density-materials-d_1652.html
            T E=35.37e6*unit_p*scale_E,nu=.3;
            Add_St_Venant_Kirchhoff_Hencky_Strain(E,nu);
            T gap=grid.dX(1)*0.1;
            RANGE<TV> box(TV(.4*m,.1*m+gap,.4*m),TV(.45*m,.1*m+gap+0.5*m,.45*m));
            Seed_Particles(box,0,0,density,particles_per_cell);
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(m/(s*s)*TV(0,-9.81,0));
        } break;
        case 20:
        case 21:
        case 22:
        case 23:
        case 24:
        case 25:
        case 26:
        case 27:
        case 28:
        case 29:{ // Mast paper
            particles.Store_Fp(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1,-0.5),TV(1.5,.1,1.5))*m,0.9); // ground 
                Add_Penalty_Collision_Object(RANGE<TV>(TV(0.06-0.3,0.6-0.5,0.5-0.5),TV(0.06+0.3,0.6+0.5,0.5+0.5))*m,0); // xmin
                Add_Penalty_Collision_Object(RANGE<TV>(TV(0.94-0.3,0.6-0.5,0.5-0.5),TV(0.94+0.3,0.6+0.5,0.5+0.5))*m,0); // xmax
                Add_Penalty_Collision_Object(RANGE<TV>(TV(0.5-0.14,0.6-0.5,0.13-0.3),TV(0.5+0.14,0.6+0.5,0.13+0.3))*m,0); // zmin
                Add_Penalty_Collision_Object(RANGE<TV>(TV(0.5-0.14,0.6-0.5,0.86-0.3),TV(0.5+0.14,0.6+0.5,0.86+0.3))*m,0);} // zmax
            else
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1,-0.5),TV(1.5,.1,1.5))*m,COLLISION_TYPE::stick,0);

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Drucker_Prager_Case(E,nu,test_number-20);
            T gap=grid.dX(1)*0.1;
            T l0=0.05*m;
            T h0=l0*8;
            CYLINDER<T> cylinder(TV(.5*m,.1+gap,.5*m),TV(.5*m,.1*m+gap+h0,.5*m),l0);
            Seed_Particles(cylinder,0,0,density,particles_per_cell);
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(m/(s*s)*TV(0,-9.81,0));
        } break;
        case 30:{ // (fluid test) pool of water 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            RANGE<TV> box(TV(0,0,0)*m,TV(1,0.25,1)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 31:{ // (fluid test) circle drop 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.75,.5)*m,.2*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 32:
        case 37:{ // sand on wedge
            particles.Store_Fp(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            ORIENTED_BOX<TV> wedge(RANGE<TV>(TV(0,0,-0.1),TV(0.2,0.2,1.1))*m,ROTATION<TV>::From_Euler_Angles(TV(0,0,0.25*M_PI)),TV(0.5,0.4-sqrt(2.0)*0.1,0)*m);
            RANGE<TV> ground(TV(-0.1,0,-0.1)*m,TV(1.1,0.1,1.1)*m);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;

            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(wedge);
                Add_Penalty_Collision_Object(ground);}
            else{
                Add_Collision_Object(ground,COLLISION_TYPE::stick,0);
                Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> >(wedge),COLLISION_TYPE::separate,1);}

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e5*unit_p*scale_E,nu=.3;
            if(test_number==32){
                int case_num=use_hardening_mast_case?hardening_mast_case:2;
                Add_Drucker_Prager_Case(E,nu,case_num);}
            //add sand particles
            RANGE<TV> box(TV(.3,.7,.3),TV(.7,.9,.7));
            Seed_Particles(box,0,0,density,particles_per_cell);
            Set_Lame_On_Particles(E,nu);

            if(test_number==37){
                if(!use_foo_T5) foo_T5=1;
                ARRAY<int> sand_particles(particles.X.m);
                for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
                Add_Drucker_Prager(E,nu,(T)35,&sand_particles);
                Add_Lambda_Particles(&sand_particles,E*foo_T5,nu,foo_T2*unit_rho,true,foo_T3,foo_T4);}

            Add_Gravity(m/(s*s)*TV(0,-9.81,0));
       } break;
        case 33:{// dry sand notch dam break, wedging friction angle
            // ./mpm 33 -3d -threads 8 -resolution 30 -last_frame 40 -framerate 48 -fooT1 4 -fooT2 30 -max_dt 1e-4 -scale_E 0.01 -symplectic_euler -no_implicit_plasticity -o notch30
            particles.Store_Fp(true);
            grid.Initialize(TV_INT(5,2,5)*resolution,RANGE<TV>(TV(0.08,0.09,0.08)*m,TV(0.38,0.21,0.38)*m),true);
            LOG::cout<<"GRID dx: "<<grid.dX<<std::endl;
            RANGE<TV> ground(TV(-10,-5,-5)*m,TV(10,0.1,10)*m);
            RANGE<TV> left_wall(TV(-5,-5,-5)*m,TV(0.1,10,10)*m);
            RANGE<TV> back_wall(TV(-5,-1,-10)*m,TV(10,10,0.1)*m);
            IMPLICIT_OBJECT_UNION<TV>* bounds=new IMPLICIT_OBJECT_UNION<TV>(
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(ground),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(left_wall),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(back_wall));
            Add_Collision_Object(bounds,COLLISION_TYPE::slip,foo_T1);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/notch.tri.gz",*surface);
            LOG::cout<<"Read mesh of notch triangle #"<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh of notch particle # "<<surface->particles.number<<std::endl;
            for(int i=0;i<surface->particles.number;i++){
                surface->particles.X(i)=(surface->particles.X(i)-TV(0.1,0.1,0.1))*0.1+TV(0.1,0.1,0.1);}
            surface->mesh.Initialize_Adjacent_Elements();
            surface->mesh.Initialize_Neighbor_Nodes();
            surface->mesh.Initialize_Incident_Elements();
            surface->Update_Bounding_Box();
            surface->Initialize_Hierarchy();
            surface->Update_Triangle_List();
            LOG::cout<<"Converting the mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Initialize_Implicit_Surface(*surface,200);
            LOG::cout<<"Seeding particles..."<<std::endl;
            Seed_Particles(*levelset,0,0,density,particles_per_cell);
            LOG::cout<<"Particle count: "<<this->particles.number<<std::endl;
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(m/(s*s)*TV(0,-9.80665,0));
            //foo_T2 is the friction angle
            Add_Drucker_Prager(E,nu,foo_T2);
        }break;
        case 34:{ // sand dam break
            // usage:./mpm 34 -3d -use_exp_F -max_dt 1e-3 -unit_p*scale_E 10 -fooT1 10 -fooT2 1000 -fooT3 3 -last_frame 20 
            particles.Store_Fp(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            RANGE<TV> ground(TV(-0.1,0,-0.1)*m,TV(1.1,0.1,1.1)*m);
            RANGE<TV> left_wall(TV(-0.1,0,-0.1)*m,TV(0.1,1.1,1.1)*m);
            RANGE<TV> back_wall(TV(-0.1,0,-0.1)*m,TV(1.1,1.1,0.1)*m);
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(ground);
                Add_Penalty_Collision_Object(left_wall);
                Add_Penalty_Collision_Object(back_wall);}
            else{
                Add_Collision_Object(ground,COLLISION_TYPE::stick,0);
                Add_Collision_Object(left_wall,COLLISION_TYPE::stick,0);
                Add_Collision_Object(back_wall,COLLISION_TYPE::stick,0);}

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e5*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            T gap=grid.dX(1)*0.1;
            RANGE<TV> box(TV(.1*m+gap,.1*m+gap,.1*m+gap),TV(.3,.75,.3)*m);
            Seed_Particles(box,0,0,density,particles_per_cell);
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles);
            //int case_num=use_hardening_mast_case?hardening_mast_case:2;
            //Add_Drucker_Prager_Case(E,nu,case_num);
            Set_Lame_On_Particles(E,nu);
            
            if(test_number==34){
                T El=500*unit_p*foo_T1,nul=0.1*foo_T3;
                Add_Lambda_Particles(&sand_particles,El,nul,foo_T2,true);}

            Add_Gravity(m/(s*s)*TV(0,-9.81,0));
        } break;
        case 35:{ // cup
            particles.Store_Fp(true);

            TV_INT res(resolution,resolution*0.3,resolution);
            TV extent(res*m);
            extent /= 2*resolution;
            grid.Initialize(res,RANGE<TV>(TV(0.25,0,0.25)*m,TV(0.25,0,0.25)*m+extent),true);
            RANGE<TV> ground(TV(-0.1,0,-0.1)*m,TV(1.1*m,grid.dX(1),1.1*m));
            if(use_penalty_collisions)
                Add_Penalty_Collision_Object(ground);
            else
                Add_Collision_Object(ground,COLLISION_TYPE::stick,0);

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            int case_num=use_hardening_mast_case?hardening_mast_case:2;
            Add_Drucker_Prager_Case(E,nu,case_num);
            CONE<T> cone(TV(0.5,grid.dX(1)/m*1.1,0.5)*m,TV(0.5,0.235+grid.dX(1)/m*1.1,0.5)*m,0.042*m);
            RANGE<TV> box(TV(0.4,0,0.4)*m,TV(0.6*m,0.095*m+grid.dX(1)*1.1,0.6*m));
            IMPLICIT_OBJECT_INTERSECTION<TV> cup(new ANALYTIC_IMPLICIT_OBJECT<CONE<T>>(cone),new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV>>(box));
            Seed_Particles(cup,0,0,density,particles_per_cell);
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(m/(s*s)*TV(0,-9.81,0));
        } break;
        case 36:{ // hourglass
            particles.Store_Fp(true);
            grid.Initialize(TV_INT(4,9,4)*resolution,RANGE<TV>(TV(-0.2,-0.45,-0.2)*m,TV(0.2,0.45,0.2)*m),true);
            LOG::cout<<"GRID DX: " <<grid.dX<<std::endl;
            IMPLICIT_OBJECT<TV>* hg=new ANALYTIC_IMPLICIT_OBJECT<HOURGLASS<TV> >(HOURGLASS<TV>(TV::Axis_Vector(1),TV(),(T).16,(T).0225,(T).8,(T).0225));
            IMPLICIT_OBJECT<TV>* inv=new IMPLICIT_OBJECT_INVERT<TV>(hg);
            if(use_penalty_collisions) Add_Penalty_Collision_Object(inv);
            else Add_Collision_Object(inv,COLLISION_TYPE::separate,.3);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            RANGE<TV> fill_part=grid.domain;
            fill_part.min_corner.y=0;
            fill_part.max_corner.y=.1;
            ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> > io_fill_part(fill_part);
            IMPLICIT_OBJECT_INTERSECTION<TV> ioi(&io_fill_part,hg);
            ioi.owns_io.Fill(false);
            Seed_Particles(ioi,0,0,density,particles_per_cell);
            LOG::printf("added %i particles.",particles.X.m);
            Set_Lame_On_Particles(E,nu);
            Add_Drucker_Prager_Case(E,nu,2);
            Add_Gravity(m/(s*s)*TV(0,-9.81,0));
        } break;
        case 38:
        case 39:{//cup flip
            particles.Store_Fp(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            grid.Initialize(TV_INT(resolution*2,resolution,2*resolution),RANGE<TV>(TV(-1,-0.5,-1),TV(1,0.5,1)),true);
            RANGE<TV> cupbottom(TV(-0.3,-0.40,-0.3),TV(0.3,-0.35,0.3));
            RANGE<TV> cupleft(TV(-0.3,-0.40,-0.3),TV(-0.25,0.25,0.3));
            RANGE<TV> cupright(TV(0.25,-0.40,-0.3),TV(0.3,0.25,0.3));
            RANGE<TV> cupback(TV(-0.3,-0.40,-0.3),TV(0.3,0.25,-0.25));
            RANGE<TV> cupfront(TV(-0.3,-0.40,0.25),TV(0.3,0.25,0.3));
            IMPLICIT_OBJECT_UNION<TV>* cup=new IMPLICIT_OBJECT_UNION<TV>(
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupbottom),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupleft),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupright),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupback),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupfront));
            Add_Collision_Object(cup,COLLISION_TYPE::separate,0);
            Add_Walls(-1,COLLISION_TYPE::slip,0.3,0.04,false);
            //Add sands
            T density;
            if(use_cohesion && sigma_Y!=0)
                density=(T)1808.89*unit_rho*scale_mass;
            else
                density=(T)1582.22*unit_rho*scale_mass;
            T E=35.37e5*unit_p*scale_E,nu=.3;
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/voronoi_fracture_sphere.tri.gz",*surface);
            LOG::cout<<"Read mesh of voronoi sphere #"<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh of voronoi sphere particle # "<<surface->particles.number<<std::endl;
            surface->mesh.Initialize_Adjacent_Elements();    
            surface->mesh.Initialize_Neighbor_Nodes();
            surface->mesh.Initialize_Incident_Elements();
            surface->Update_Bounding_Box();
            surface->Initialize_Hierarchy();
            surface->Update_Triangle_List();
            LOG::cout<<"Converting the mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Initialize_Implicit_Surface(*surface,200);
            LOG::cout<<"Seeding particles..."<<std::endl;
            Seed_Particles(*levelset,0,0,density,particles_per_cell);
            LOG::cout<<"Particle count: "<<this->particles.number<<std::endl;
            //Sand properties
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles,false,test_number==38?sigma_Y:0);
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(m/(s*s)*TV(0,-9.81,0));
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
            for(int i=0;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(.8,.7,.7);
            //Add water particles for case 39
            if(test_number==39){
                if(!use_foo_T5) foo_T5=(T)1e-2;
                Add_Lambda_Particles(&sand_particles,E*foo_T5,nu,(T)1000*unit_rho,true,(T)0.3,(T)1);}
            int add_gravity_frame=restart?restart:10;
            begin_frame=[this,add_gravity_frame](int frame){if(frame==add_gravity_frame) Add_Gravity(TV(0,20,0));};
        } break;
        case 48:{ // sand ball
            particles.Store_Fp(true); if(!no_implicit_plasticity) use_implicit_plasticity=true;
            grid.Initialize(TV_INT(resolution,resolution,resolution),RANGE<TV>(TV(-.05,-.07,-.05),TV(.05,.03,.05)),true);
            //Add sands
            T density;
            if(use_cohesion && sigma_Y!=0)
                density=(T)1808.89*unit_rho*scale_mass;
            else
                density=(T)1582.22*unit_rho*scale_mass;
            T E=35.37e4*unit_p*scale_E,nu=.3;

            RANGE<TV> ground(TV(-10,-10,-10)*m,TV(10,-0.05,10)*m);
            Add_Collision_Object(ground,COLLISION_TYPE::stick,0);
            //strong levelset
            TRIANGULATED_SURFACE<T>* strong=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/voronoi_strong_50.tri.gz",*strong);
            LOG::cout<<"Read mesh of strong #"<<strong->mesh.elements.m<<std::endl;
            LOG::cout<<"Read particles of strong #"<<strong->particles.number<<std::endl;
            for(int i=0;i<strong->particles.number;i++) strong->particles.X(i)/=20;
            strong->mesh.Initialize_Adjacent_Elements();    
            strong->mesh.Initialize_Neighbor_Nodes();
            strong->mesh.Initialize_Incident_Elements();
            strong->Update_Bounding_Box();
            strong->Initialize_Hierarchy();
            strong->Update_Triangle_List();
            LOG::cout<<"Converting strong mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* strong_levelset=Initialize_Implicit_Surface(*strong,200);
            //full sphere
            TRIANGULATED_SURFACE<T>* full=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/voronoi_full_50.tri.gz",*full);
            LOG::cout<<"Read mesh of full #"<<full->mesh.elements.m<<std::endl;
            LOG::cout<<"Read particles of full #"<<full->particles.number<<std::endl;
            for(int i=0;i<full->particles.number;i++) full->particles.X(i)/=20;
            full->mesh.Initialize_Adjacent_Elements();    
            full->mesh.Initialize_Neighbor_Nodes();
            full->mesh.Initialize_Incident_Elements();
            full->Update_Bounding_Box();
            full->Initialize_Hierarchy();
            full->Update_Triangle_List();
            LOG::cout<<"Converting the mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* full_levelset=Initialize_Implicit_Surface(*full,200);
            LOG::cout<<"Seeding particles..."<<std::endl;
            Seed_Particles(*full_levelset,0,0,density,particles_per_cell);
            LOG::cout<<"Sand particle count: "<<this->particles.number<<std::endl;
            LOG::printf("Sand particle count=%P\n",particles.X.m);

            //Sand properties
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles,false,0);
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));

            T porosity=0.3;
            T saturation_level=1;
            T water_density=(T)1000*unit_rho;
            T water_E=E*foo_T5;
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
            ARRAY<int> strong_lambda_particles;
            ARRAY<int> weak_lambda_particles;
            T volume_lambda=particles.volume(0)*porosity*saturation_level;
            T mass_lambda=water_density*volume_lambda;
            T lambda=water_E*nu/((1+nu)*(1-2*nu));
            for(int k=0;k<sand_particles.m;k++){
                int p=particles.Add_Element();
                if(strong_levelset->Extended_Phi(particles.X(k))<=0){
                    strong_lambda_particles.Append(p);
                    (*color_attribute)(p)=VECTOR<T,3>(1,0,0);
                    (*color_attribute)(k)=VECTOR<T,3>(1,0,0);}
                else{
                    weak_lambda_particles.Append(p);
                    (*color_attribute)(p)=VECTOR<T,3>(0,1,0);
                    (*color_attribute)(k)=VECTOR<T,3>(0,1,0);}
                //lambda_particles(k)=p;
                int i=sand_particles(k);
                particles.valid(p)=true;
                particles.X(p)=particles.X(i);
                particles.V(p)=particles.V(i);
                particles.F(p)=particles.F(i);
                if(particles.store_Fp) particles.Fp(p)=particles.Fp(i); 
                if(particles.store_B) particles.B(p)=particles.B(i);
                if(particles.store_C) particles.C(p)=particles.C(i);
                if(particles.store_S) particles.S(p)=particles.S(i);
                particles.mass(p)=mass_lambda;
                particles.volume(p)=volume_lambda;
                particles.mu(p)=(T)0;
                particles.mu0(p)=(T)0;
                particles.lambda(p)=lambda;
                particles.lambda0(p)=lambda;}
            Add_Fixed_Corotated(water_E,nu,&strong_lambda_particles,true);
            Add_Fixed_Corotated(water_E*0.1,nu,&weak_lambda_particles,true);
            //if strong level set extended phi<0 then it means that we are in the strong level set
        } break;

        case 40:{ // dry sand siggraph letters drop
            particles.Store_Fp(true);
            grid.Initialize(TV_INT(6,1,2)*resolution,RANGE<TV>(TV(-0.75,0,-0.25)*m,TV(0.75,0.25,0.25)*m),true);
            LOG::cout<<"GRID dx: "<<grid.dX<<std::endl;
            RANGE<TV> ground(TV(-10,-10,-10)*m,TV(10,0.05,10)*m);
            if(use_penalty_collisions) Add_Penalty_Collision_Object(ground);
            else Add_Collision_Object(ground,COLLISION_TYPE::slip,10);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            int case_num=use_hardening_mast_case?hardening_mast_case:2;
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/siggraph_letters.tri.gz",*surface);
            LOG::cout<<"Read mesh of siggraph letters triangle #"<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh of siggraph letters particle # "<<surface->particles.number<<std::endl;
            surface->mesh.Initialize_Adjacent_Elements();    
            surface->mesh.Initialize_Neighbor_Nodes();
            surface->mesh.Initialize_Incident_Elements();
            surface->Update_Bounding_Box();
            surface->Initialize_Hierarchy();
            surface->Update_Triangle_List();
            LOG::cout<<"Converting the mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Initialize_Implicit_Surface(*surface,200);
            LOG::cout<<"Seeding particles..."<<std::endl;
            Seed_Particles(*levelset,0,0,density,particles_per_cell);
            LOG::cout<<"Particle count: "<<this->particles.number<<std::endl;
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(m/(s*s)*TV(0,-9.80665,0));
            Add_Drucker_Prager_Case(E,nu,case_num);
        } break;
        case 41:{ // Draw in sand
            particles.Store_Fp(true);
            grid.Initialize(TV_INT(4,1,4)*resolution,RANGE<TV>(TV(),TV(1,0.25,1)*m),true);
            LOG::printf("REAL GRID: %P\n",grid);

            if(!friction_is_set)friction=0.5;
            IMPLICIT_OBJECT<TV> *sandbox=Levelset_From_File<T>(data_directory+"/../Private_Data/sandbox_only.tri.gz");
            Add_Collision_Object(sandbox,COLLISION_TYPE::stick,friction);


            const T stylus_r=0.01*m;
            const T settle_wait(0.1);
            const T final_t(10-settle_wait);
            Add_Collision_Object(
                    CYLINDER<T>(TV(0,0.056*m,0),TV(0,0.3*m,0),stylus_r),COLLISION_TYPE::separate,friction,
                    [=](T time){
                        T t=std::max((T)0,time-settle_wait)*2*pi/final_t;
                        T x=-(4*cos(t)*cos(16*t)+17*cos(t)*sin(7*t)-20*cos(t)*sin(5*t)+20*cos(t)*cos(4*t)-25*cos(t)*sin(3*t)-30*cos(t)*cos(2*t)+5*cos(t)*sin(t)-90*cos(t)-150)/300.0*m;
                        T y=-(4*sin(t)*cos(16*t)+17*sin(t)*sin(7*t)-20*sin(t)*sin(5*t)+20*sin(t)*cos(4*t)-25*sin(t)*sin(3*t)-30*sin(t)*cos(2*t)+5*sqr(sin(t))-90*sin(t)-150)/300.0*m;
                        return FRAME<TV>(TV(x,0,y));},
                    [=](T time){
                        if(time<settle_wait) return TWIST<TV>();
                        T t=(time-settle_wait)*2*pi/final_t;
                        T dx=-(-64*cos(t)*sin(16*t)-4*sin(t)*cos(16*t)-17*sin(t)*sin(7*t)+119*cos(t)*cos(7*t)+20*sin(t)*sin(5*t)-100*cos(t)*cos(5*t)-80*cos(t)*sin(4*t)-20*sin(t)*cos(4*t)+25*sin(t)*sin(3*t)-75*cos(t)*cos(3*t)+60*cos(t)*sin(2*t)+30*sin(t)*cos(2*t)-5*sqr(sin(t))+90*sin(t)+5*sqr(cos(t)))/300.0;
                        T dy=-(-64*sin(t)*sin(16*t)+4*cos(t)*cos(16*t)+17*cos(t)*sin(7*t)+119*sin(t)*cos(7*t)-20*cos(t)*sin(5*t)-100*sin(t)*cos(5*t)-80*sin(t)*sin(4*t)+20*cos(t)*cos(4*t)-25*cos(t)*sin(3*t)-75*sin(t)*cos(3*t)+60*sin(t)*sin(2*t)-30*cos(t)*cos(2*t)+10*cos(t)*sin(t)-90*cos(t))/300.0;
                        return TWIST<TV>(TV(dx,0,dy),typename TV::SPIN());});

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;

            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Levelset_From_File<T>(data_directory+"/../Private_Data/sanddune_square.tri.gz");
            Seed_Particles(*levelset,0,0,density,particles_per_cell);
            LOG::printf("Particle count: %d\n",particles.number);
            Set_Lame_On_Particles(E,nu);

            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            int case_num=use_hardening_mast_case?hardening_mast_case:2;
            Add_Drucker_Prager_Case(E,nu,case_num);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 42:{ // Raking
            particles.Store_Fp(true);
            grid.Initialize(TV_INT(4,1,4)*resolution,RANGE<TV>(TV(),TV(1,0.25,1)*m),true);
            LOG::printf("REAL GRID: %P\n",grid);

            if(!friction_is_set)friction=0.5;
            const TV start_pos(0.25*m,0,0);

            IMPLICIT_OBJECT<TV> *sandbox=Levelset_From_File<T>(data_directory+"/../Private_Data/sandbox.tri.gz");
            Add_Collision_Object(sandbox,COLLISION_TYPE::stick,friction);

            LEVELSET_IMPLICIT_OBJECT<TV>* rake=Levelset_From_File<T>(data_directory+"/../Private_Data/rake.tri.gz");

            const T settle_wait(0.1);
            const T final_t(10);
            Add_Collision_Object(rake,COLLISION_TYPE::separate,friction,
                    [=](T time){
                        time=max((T)0,time-settle_wait);
                        ROTATION<TV> R=ROTATION<TV>::From_Euler_Angles(0,2*pi*time/final_t,0);
                        return FRAME<TV>(TV(0.5,0.056,0.5)*m+R.Rotate(start_pos),R);},
                    [=](T time){return TWIST<TV>(TV(),typename TV::SPIN(0,time<settle_wait?0:2*pi/final_t,0));});

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;

            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Levelset_From_File<T>(data_directory+"/../Private_Data/sanddune_square.tri.gz");
            Seed_Particles(*levelset,0,0,density,particles_per_cell);
            LOG::printf("Particle count: %d\n",particles.number);
            Set_Lame_On_Particles(E,nu);

            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            int case_num=use_hardening_mast_case?hardening_mast_case:2;
            Add_Drucker_Prager_Case(E,nu,case_num);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 43:{ // Rotating table
            particles.Store_Fp(true);
            grid.Initialize(TV_INT(4,1,4)*resolution,RANGE<TV>(TV(),TV(1,0.25,1)*m),true);
            LOG::printf("REAL GRID: %P\n",grid);

            if(!friction_is_set)friction=0.5;
            const T omega=0.5*pi/s;
            Add_Collision_Object(
                    CYLINDER<T>(TV(0,-0.1,0),TV(0,-0.1,0),0.5),
                    COLLISION_TYPE::separate,friction,
                    [=](T time){return FRAME<TV>(TV(),ROTATION<TV>::From_Euler_Angles(0,omega*time,0));},
                    [=](T time){return TWIST<TV>(TV(),typename TV::SPIN(0,omega,0));});
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;

        case 44:{ // sand falling into a pile.
            particles.Store_Fp(true);
            grid.Initialize(TV_INT(5,2,5)*resolution,RANGE<TV>(TV(-.1,-0.02,-0.1),TV(0.1,0.06,0.1))*m,true);
            RANGE<TV> ground(TV(-10,-10,-10)*m,TV(10,0,10)*m);
            if(use_penalty_collisions) Add_Penalty_Collision_Object(ground);
            else Add_Collision_Object(ground,COLLISION_TYPE::stick,0);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e4*unit_p*scale_E,nu=.3;
            T spout_width=8.334e-3*m;
            T spout_height=.01*m;
            T seed_buffer=grid.dX.y*5;
            T pour_speed=.16319*m/s;
            TV gravity=TV(0,-9.8*m/(s*s),0);
            CYLINDER<T> seed_range(TV(0,0.06*m-spout_height,0),TV(0,0.06*m+seed_buffer,0),spout_width/2);

            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            POUR_SOURCE<TV>* source=new POUR_SOURCE<TV>(*this,
                *new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(seed_range),TV(0,-1,0),grid.domain.max_corner,
                TV(0,-pour_speed,0),gravity,max_dt*pour_speed+grid.dX.y,seed_buffer,mass,volume);
            destroy=[=](){delete source;};
            write_output_files=[=](int frame){source->Write_Output_Files(frame);};
            read_output_files=[=](int frame){source->Read_Output_Files(frame);};
            begin_time_step=[=](T time)
                {
                    if(time<0.08||time>=foo_T3) return;
                    ARRAY<int> affected_particles;
                    int n=particles.number;
                    source->Begin_Time_Step(time);
                    T mu=E/(2*(1+nu));
                    T lambda=E*nu/((1+nu)*(1-2*nu));
                    for(int i=n;i<particles.number;i++){
                        particles.mu(i)=mu;
                        particles.mu0(i)=mu;
                        particles.lambda(i)=lambda;
                        particles.lambda0(i)=lambda;
                        affected_particles.Append(n);}
                    for(int i=0;i<plasticity_models.m;i++)
                        if(MPM_DRUCKER_PRAGER<TV>* dp=dynamic_cast<MPM_DRUCKER_PRAGER<TV>*>(plasticity_models(i)))
                            dp->Initialize_Particles(&affected_particles);
                
                };
            end_time_step=[=](T time){
                if(time<=foo_T3)
                    source->End_Time_Step(time);
            };



            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            //int case_num=use_hardening_mast_case?hardening_mast_case:2;

            //ARRAY<int> sand_particles(particles.X.m);
            //for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;

            // Add_Drucker_Prager(E,nu,foo_T2,&sand_particles);
            Add_Drucker_Prager(E,nu,foo_T2);

            // Add_Drucker_Prager_Case(E,nu,case_num);
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(gravity);
        } break;

         case 45:{ // sand castle
             // ./mpm 45 -3d -threads 8 -resolution 20 -max_dt 5e-5 -scale_E 0.01 -mast_case 0 -last_frame 120 -fooT2 10 -framerate 72 -symplectic_euler -no_implicit_plasticity
            particles.Store_Fp(true);
            grid.Initialize(TV_INT(15,7,10)*resolution,RANGE<TV>(TV(-1.0,-0.1,-1.0)*m,TV(2.0,1.3,1.0)*m),true);
            LOG::cout<<"GRID dx: "<<grid.dX<<std::endl;
            RANGE<TV> boxymin(TV(-10,-10,-10)*m,TV(10,0,10)*m);
            RANGE<TV> boxymax(TV(-10,1.2,-10)*m,TV(10,10,10)*m);
            RANGE<TV> boxxmin(TV(-10,-10,-10)*m,TV(-0.9,10,10)*m);
            RANGE<TV> boxxmax(TV(1.9,-10,-10)*m,TV(10,10,10)*m);
            RANGE<TV> boxzmin(TV(-10,-10,-10)*m,TV(10,10,-0.9)*m);
            RANGE<TV> boxzmax(TV(-10,-10,0.9)*m,TV(10,10,10)*m);
            IMPLICIT_OBJECT_UNION<TV>* box=new IMPLICIT_OBJECT_UNION<TV>(
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(boxxmin),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(boxxmax),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(boxymin),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(boxymax),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(boxzmin),
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(boxzmax));
            if(use_penalty_collisions) PHYSBAM_FATAL_ERROR();
            else Add_Collision_Object(box,COLLISION_TYPE::stick,0);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            int case_num=use_hardening_mast_case?hardening_mast_case:2;
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/castle.tri.gz",*surface);
            LOG::cout<<"Read mesh of castle triangle #"<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh of castle particle # "<<surface->particles.number<<std::endl;
            surface->mesh.Initialize_Adjacent_Elements();    
            surface->mesh.Initialize_Neighbor_Nodes();
            surface->mesh.Initialize_Incident_Elements();
            surface->Update_Bounding_Box();
            surface->Initialize_Hierarchy();
            surface->Update_Triangle_List();
            LOG::cout<<"Converting the mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Initialize_Implicit_Surface(*surface,200);
            LOG::cout<<"Seeding particles..."<<std::endl;
            Seed_Particles(*levelset,0,0,density,particles_per_cell);
            LOG::cout<<"Particle count: "<<this->particles.number<<std::endl;
            Set_Lame_On_Particles(E,nu);
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager_Case(E,nu,case_num,&sand_particles);
            SPHERE<TV> sphere1(TV(-0.64,0.52,0)*m,.15*m);
            T density_sphere=1e4*unit_rho*scale_mass;
            VECTOR<T,3> angular_velocity1(TV(0,0,foo_T2));
            Seed_Particles(sphere1,[=](const TV& X){return angular_velocity1.Cross(X-sphere1.center)+TV(5,-1.5,0)*(m/s);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity1);},density_sphere,8);
            Add_Fixed_Corotated(40e5*unit_p,0.3);
            Add_Gravity(m/(s*s)*TV(0,-9.80665,0));
        } break;    

        case 46:{ // sandbox with simulated ball
            if(!use_foo_T1) foo_T1=1;
            if(!use_foo_T2) foo_T2=1;
            if(!use_foo_T3) foo_T3=1;
            particles.Store_Fp(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            T sand_depth=0.35*m;
            T air_height=0.65*m;
            T sand_length=1*m;
            T sand_width=1*m;
            T wall_thickness=0.2*m;
            RANGE<TV> sand(TV(-sand_width/2,-sand_depth,-sand_length/2),TV(sand_width/2,0,sand_length/2));
            RANGE<TV> sand_box(sand.Thickened(wall_thickness));
            RANGE<TV> sand_box_interior(sand);
            sand_box_interior.max_corner(1)=sand_box.max_corner(1);
            RANGE<TV> domain(sand);
            domain.max_corner(1)=air_height;
            domain=domain.Thickened(0.05*m);
            grid.Initialize(TV_INT(domain.Edge_Lengths()*resolution),domain,true);
            auto sand_box_walls=
                new IMPLICIT_OBJECT_INTERSECTION<TV>(
                        new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(sand_box),
                        new IMPLICIT_OBJECT_INVERT<TV>(
                            new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(sand_box_interior)));
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(sand_box_walls);}
            else{
                Add_Collision_Object(sand_box_walls,COLLISION_TYPE::stick,0);}

            T density=(T)1582*unit_rho;
            T E=35.37e6*unit_p*scale_E,nu=.3;

            Seed_Particles(sand,0,0,density,particles_per_cell); 
            // SAND PROPERTIES
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)21,&sand_particles);
            Set_Lame_On_Particles(E,nu);
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
            T g=9.81*m/(s*s);
            for(int i=0;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(.8,.7,.7);
            {
                int N_sand=particles.number;
                VECTOR<T,3> velocity(TV(5,-5,0)*foo_T1*m/s);
                VECTOR<T,3> angular_velocity(TV(0,0,10)*foo_T2/s);
                T y0=0.15*m;
                T t_impact=(velocity(1)+sqrt(sqr(velocity(1))+2*y0*g))/g;
                SPHERE<TV> sphere(TV(-velocity(0)*t_impact,y0,-velocity(2)*t_impact),0.1*foo_T3*m);
                T density=11340*unit_rho*scale_mass;
                Seed_Particles(sphere,[=](const TV& X){return velocity+angular_velocity.Cross(X-sphere.center);},
                        [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity);}
                        ,density,particles_per_cell);
                int N_box_particles=particles.number-N_sand;
                LOG::cout<<N_sand<<" sand particles "<<N_box_particles<<" ball particles\n";
                ARRAY<int> ball_particles(N_box_particles);
                for(int k=0;k<ball_particles.m;k++) ball_particles(k)=k+N_sand;
                Add_Fixed_Corotated(16e7*unit_p*scale_E,0.44,&ball_particles);
            }
            Add_Gravity(TV(0,-g,0));
        } break;

        case 47:{ // sandbox with collision object ball
            if(!use_foo_T1) foo_T1=1;
            if(!use_foo_T2) foo_T2=1;
            if(!use_foo_T3) foo_T3=1;
            particles.Store_Fp(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            T sand_depth=0.35*m;
            T air_height=0.65*m;
            T sand_length=1*m;
            T sand_width=1*m;
            T wall_thickness=0.2*m;
            RANGE<TV> sand(TV(-sand_width/2,-sand_depth,-sand_length/2),TV(sand_width/2,0,sand_length/2));
            RANGE<TV> domain(sand);
            domain.max_corner(1)=air_height;
            RANGE<TV> sand_box(domain.Thickened(wall_thickness));
            RANGE<TV> sand_box_interior(domain);
            domain=domain.Thickened(0.05*m);
            grid.Initialize(TV_INT(domain.Edge_Lengths()*resolution),domain,true);
            auto sand_box_walls=
                new IMPLICIT_OBJECT_INTERSECTION<TV>(
                        new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(sand_box),
                        new IMPLICIT_OBJECT_INVERT<TV>(
                            new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(sand_box_interior)));
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(sand_box_walls);}
            else{
                Add_Collision_Object(sand_box_walls,COLLISION_TYPE::slip,0);}

            T y0=0.15*m;
            T v0=-1*m/s*foo_T1;
            T g=9.81*m/(s*s);
            T radius=0.1*foo_T3*m;
            SPHERE<TV> sphere(TV(0,y0,0),radius);
            T stop_time=(v0+sqrt(sqr(v0)+2*(y0-radius+sand_depth-0.05*m)*g))/g;
            Add_Collision_Object(sphere,COLLISION_TYPE::separate,friction,
                    [=](T time){
                        if(time>stop_time) time=stop_time;
                        return FRAME<TV>(TV(0,v0*time-0.5*g*sqr(time),0));},
                    [=](T time){
                        if(time>=stop_time) return TWIST<TV>(TV(),typename TV::SPIN());
                        return TWIST<TV>(TV(0,v0-g*time,0),typename TV::SPIN());});

            T density=(T)1582*unit_rho;
            T E=35.37e6*unit_p*scale_E,nu=.2;

            Seed_Particles(sand,0,0,density,particles_per_cell); 
            // SAND PROPERTIES
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)21,&sand_particles);
            Set_Lame_On_Particles(E,nu);
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
            for(int i=0;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(.8,.7,.7);
            Add_Gravity(TV(0,-g,0));
        } break;

        case 69:{ // Magic balls
            particles.Store_Fp(true);
            particles.Store_Lame(true);
            grid.Initialize(TV_INT(1,1,1)*resolution,RANGE<TV>::Unit_Box()*m,true);
            LOG::cout<<"GRID dx: "<<grid.dX<<std::endl;
            RANGE<TV> boxymin(TV(-10,-10,-10)*m,TV(10,0.03,10)*m);
            if(use_penalty_collisions) PHYSBAM_FATAL_ERROR();
            else Add_Collision_Object(boxymin,COLLISION_TYPE::slip,0.7);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;

            const SPHERE<TV> ball(TV(0,0.75,0)*m,.1*m);
            const SPHERE<TV> padded_ball(ball.center,0.11*m);
            const T ball_volume=ball.Size();
            const int number_of_particles=particles_per_cell*ball_volume/grid.dX.Product();
            const T volume_per_particle=ball_volume/number_of_particles;

            ARRAY<SPHERE<TV> > spheres;
            for(int i=1;i<=500;i++){
                TV offset=random.template Get_Vector_In_Unit_Sphere<TV>()*padded_ball.radius;
                if(offset.Magnitude()<padded_ball.radius*0.5 || offset.Magnitude()>padded_ball.radius*0.9)
                {i--;continue;}
                T radius=random.Get_Uniform_Number(0,padded_ball.radius-offset.Magnitude());
                if(radius<padded_ball.radius*0.1)
                {i--;continue;}
                spheres.Append(SPHERE<TV>(offset+padded_ball.center,radius));}

            const T h=10;
            const T a=1.025;
            const T ea=exp((a-1)*h); 
            particles.Preallocate(number_of_particles);
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
            for(int p=0;p<number_of_particles;p++){
                particles.valid(p)=true;
                particles.X(p)=random.template Get_Vector_In_Unit_Sphere<TV>()*padded_ball.radius+padded_ball.center;
                particles.F(p)=MATRIX<T,TV::m>()+1;
                if(particles.store_Fp) particles.Fp(p).Set_Identity_Matrix();

                bool inside=ball.Inside(particles.X(p),0);
                for(int j=0;j<spheres.m;j++)
                    inside|=spheres(j).Inside(particles.X(p),0);
                if(!inside){p--;continue;}

                particles.V(p)=TV(0,-1.2,0);
                particles.mass(p)=density*volume_per_particle*Perturb(0.5);
                particles.volume(p)=volume_per_particle;

                (*color_attribute)(p)=VECTOR<T,3>(1,1,1);
                T E_local=E*Perturb(0.5);
                T nu_local=nu*Perturb(0.1);
                particles.mu(p)=E_local/(2*(1+nu_local));
                particles.lambda(p)=E_local*nu_local/((1+nu_local)*(1-2*nu_local));

                T mult=1;
                for(int i=0;i<10;i++)
                    if((particles.X(p)-ball.center).Magnitude()>(i+1)*ball.radius/10){
                        particles.mass(p)*=a;
                        particles.lambda(p)*=ea;
                        particles.mu(p)*=ea;}
                    else mult/=a;

                if(Uniform(0,1)>mult) p--;}

            for(int k=0;k<500;k++){
                const TV offset=random.template Get_Vector_In_Unit_Sphere<TV>()*padded_ball.radius*5;
                const T radius=random.Get_Uniform_Number(padded_ball.radius,padded_ball.radius*10);
                SPHERE<TV> sphere_big(offset+padded_ball.center,radius+grid.dX.Min()/2*2);
                SPHERE<TV> sphere_small(offset+padded_ball.center,radius-grid.dX.Min()/2*2);
                const TV offset_new=random.template Get_Vector_In_Unit_Sphere<TV>()*padded_ball.radius*5;
                const T radius_new=random.Get_Uniform_Number(padded_ball.radius,padded_ball.radius*10);
                SPHERE<TV> sphere_new(offset_new+padded_ball.center,radius_new);
                if(offset_new.Magnitude()+ball.radius<radius_new) continue;
                if((offset-offset_new).Magnitude()+radius<radius_new) continue;
                for(int i=0;i<particles.X.m;i++)
                    if(sphere_big.Inside(particles.X(i),0)&&!sphere_small.Inside(particles.X(i),0)&&sphere_new.Inside(particles.X(i),0)){
                        particles.lambda(i)*=0.5;
                        particles.mu(i)*=0.5;
                        (*color_attribute)(i)(2)*=0.5;}}
            

            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles,false,0);
            Add_Lambda_Particles(&sand_particles,E,nu,unit_rho,true,0.3,1.0);
            Add_Gravity(TV(0,-9.81,0));
        }

        case 948:{ // column collapse wedge firction angle
            particles.Store_Fp(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1,-0.5),TV(1.5,.1,1.5))*m,0.9); // ground 
                Add_Penalty_Collision_Object(RANGE<TV>(TV(0.06-0.3,0.6-0.5,0.5-0.5),TV(0.06+0.3,0.6+0.5,0.5+0.5))*m,0); // xmin
                Add_Penalty_Collision_Object(RANGE<TV>(TV(0.94-0.3,0.6-0.5,0.5-0.5),TV(0.94+0.3,0.6+0.5,0.5+0.5))*m,0); // xmax
                Add_Penalty_Collision_Object(RANGE<TV>(TV(0.5-0.14,0.6-0.5,0.13-0.3),TV(0.5+0.14,0.6+0.5,0.13+0.3))*m,0); // zmin
                Add_Penalty_Collision_Object(RANGE<TV>(TV(0.5-0.14,0.6-0.5,0.86-0.3),TV(0.5+0.14,0.6+0.5,0.86+0.3))*m,0);} // zmax
            else
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1,-0.5),TV(1.5,.1,1.5))*m,COLLISION_TYPE::stick,0);

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Drucker_Prager_Case(E,nu,test_number-20);
            T gap=grid.dX(1)*0.1;
            T l0=0.05*m;
            T h0=l0*8;
            CYLINDER<T> cylinder(TV(.5*m,.1+gap,.5*m),TV(.5*m,.1*m+gap+h0,.5*m),l0);
            Seed_Particles(cylinder,0,0,density,particles_per_cell);
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(m/(s*s)*TV(0,-9.81,0));
        } break;


        case 949:{ // lambda voronoi sand ball
            // ./mpm 949 -3d -resolution 60 -threads 1 -max_dt 1e-4 -scale_E 0.01 -framerate 120 -last_frame 20 -fooT3 1 -fooT5 5 -symplectic_euler -no_implicit_plasticity -o zz
            particles.Store_Fp(true);
            grid.Initialize(TV_INT(4,1,4)*resolution,RANGE<TV>(TV(-0.5,0,-0.5)*m,TV(0.5,0.25,0.5)*m),true);
            LOG::cout<<"GRID dx: "<<grid.dX<<std::endl;
            RANGE<TV> boxymin(TV(-10,-10,-10)*m,TV(10,0.03,10)*m);
            if(use_penalty_collisions) PHYSBAM_FATAL_ERROR();
            else Add_Collision_Object(boxymin,COLLISION_TYPE::slip,0.7);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;

            // SEEDING FULL
            TRIANGULATED_SURFACE<T>* surface_strong=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/voronoi_strong_50.tri.gz",*surface_strong);
            LOG::cout<<"Read mesh of strong voronoi triangle #"<<surface_strong->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh of sterong voronoi particle # "<<surface_strong->particles.number<<std::endl;
            for(int i=0;i<surface_strong->particles.number;i++){
                surface_strong->particles.X(i)/=3;
                surface_strong->particles.X(i)+=TV(0,0.15,0);}
            surface_strong->mesh.Initialize_Adjacent_Elements();    
            surface_strong->mesh.Initialize_Neighbor_Nodes();
            surface_strong->mesh.Initialize_Incident_Elements();
            surface_strong->Update_Bounding_Box();
            surface_strong->Initialize_Hierarchy();
            surface_strong->Update_Triangle_List();
            LOG::cout<<"Converting the mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* strong_levelset=Initialize_Implicit_Surface(*surface_strong,200);
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/voronoi_full_50.tri.gz",*surface);
            LOG::cout<<"Read mesh of full sandball triangle #"<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh of full sandbval particle # "<<surface->particles.number<<std::endl;
            for(int i=0;i<surface->particles.number;i++){
                surface->particles.X(i)/=3;
                surface->particles.X(i)+=TV(0,0.15,0);}
            surface->mesh.Initialize_Adjacent_Elements();    
            surface->mesh.Initialize_Neighbor_Nodes();
            surface->mesh.Initialize_Incident_Elements();
            surface->Update_Bounding_Box();
            surface->Initialize_Hierarchy();
            surface->Update_Triangle_List();
            LOG::cout<<"Converting the full mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Initialize_Implicit_Surface(*surface,200);
            LOG::cout<<"Seeding particles..."<<std::endl;
            Seed_Particles(*levelset,[=](const TV& X){return TV(0,-6,0);},0,density,particles_per_cell);
            LOG::cout<<"Particle count: "<<this->particles.number<<std::endl;
            
            // SAND
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles,false,0);
            Set_Lame_On_Particles(E,nu);

            // GRAVITY
            Add_Gravity(m/(s*s)*TV(0,-9.80665,0));

            // LAMBDA
            T porosity=0.3;
            T saturation_level=foo_T3; // 0 - 1
            T water_density=(T)1000*unit_rho;
            T water_E=E*foo_T5; // 0.5
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
            ARRAY<int> strong_lambda_particles;
            ARRAY<int> weak_lambda_particles;
            T volume_lambda=particles.volume(0)*porosity*saturation_level;
            T mass_lambda=water_density*volume_lambda;
            T lambda=water_E*nu/((1+nu)*(1-2*nu));
            for(int k=0;k<sand_particles.m;k++){
                int i=sand_particles(k);
                if(strong_levelset->Extended_Phi(particles.X(k))<=0){
                    int p=particles.Add_Element();
                    particles.mass(p)=mass_lambda;
                    strong_lambda_particles.Append(p);
                    particles.lambda(p)=lambda;
                    particles.lambda0(p)=lambda;
                    (*color_attribute)(p)=VECTOR<T,3>(1,0,0);
                    (*color_attribute)(k)=VECTOR<T,3>(1,0,0);
                    particles.valid(p)=true;
                    particles.X(p)=particles.X(i);
                    particles.V(p)=particles.V(i);
                    particles.F(p)=particles.F(i);
                    particles.mass(i)*=0.99; // fanfu hack
                    if(particles.store_Fp) particles.Fp(p)=particles.Fp(i); 
                    if(particles.store_B) particles.B(p)=particles.B(i);
                    if(particles.store_C) particles.C(p)=particles.C(i);
                    if(particles.store_S) particles.S(p)=particles.S(i);
                    particles.volume(p)=volume_lambda;
                    particles.mu(p)=(T)0;
                    particles.mu0(p)=(T)0;}
                else{
                    (*color_attribute)(k)=VECTOR<T,3>(0,1,0);}}
            Add_Fixed_Corotated(water_E,nu,&strong_lambda_particles,true);

        } break;

        case 950:{ // kdtree wet sand ball (filled with weak)
            // ./mpm 950 -3d -resolution 50 -threads 10 -max_dt 1e-4 -framerate 120 -last_frame 120 -fooT1 0.000001 -fooT2 1000 -fooT4 0.2 -symplectic_euler -no_implicit_plasticity -o bbb
            particles.Store_Fp(true);
            particles.Store_Lame(true);
            grid.Initialize(TV_INT(4,1,4)*resolution,RANGE<TV>(TV(-0.5,0,-0.5)*m,TV(0.5,0.25,0.5)*m),true);
            LOG::cout<<"GRID dx: "<<grid.dX<<std::endl;

            RANGE<TV> boxymin(TV(-10,-10,-10)*m,TV(10,0.03,10)*m);
            // RANGE<TV> boxzmin(TV(-10,-10,-10)*m,TV(10,10,-0.47)*m);
            // IMPLICIT_OBJECT_UNION<TV>* bounds=new IMPLICIT_OBJECT_UNION<TV>(
            //     new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(boxymin),
            //     new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(boxzmin));

            if(use_penalty_collisions) PHYSBAM_FATAL_ERROR();
            else Add_Collision_Object(boxymin,COLLISION_TYPE::slip,0.3);

            T density_sand=(T)2200*unit_rho*scale_mass;
            // T density_water=(T)1000*unit_rho;
            T nu=0.3;
            T E_strong_sand=35.37e4*unit_p*scale_E;
            // T E_strong_water=E_strong_sand*foo_T5; // foo_T5 = 5 is the water youngs moudulus over sand youngs modulus
            T E_weak_sand=E_strong_sand*foo_T1; // foo_T1 is the softening ratio
            // T E_weak_water=E_strong_water*foo_T1;
            // T porosity=0.3;
            // T saturation_level=foo_T3; // foo_T3 = 0 - 1 is saturation
            T ratio_of_max_distance_considered_to_be_weak=foo_T4; // foo_T4 is the thickness of weak layer (0 - 1)
            T cohesion=foo_T2; // foo_T2 is cohesion

            if(!no_implicit_plasticity) use_implicit_plasticity=true;

            // SEEDING ALL SAND
            TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
            FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.f),data_directory+"/../Private_Data/voronoi_full_100.tri.gz",*surface);
            LOG::cout<<"Read mesh of full sandball triangle #"<<surface->mesh.elements.m<<std::endl;
            LOG::cout<<"Read mesh of full sandbval particle # "<<surface->particles.number<<std::endl;
            for(int i=0;i<surface->particles.number;i++){surface->particles.X(i)/=3;surface->particles.X(i)+=TV(0,0.15,0);}
            surface->mesh.Initialize_Adjacent_Elements();    
            surface->mesh.Initialize_Neighbor_Nodes();
            surface->mesh.Initialize_Incident_Elements();
            surface->Update_Bounding_Box();
            surface->Initialize_Hierarchy();
            surface->Update_Triangle_List();
            LOG::cout<<"Converting the full mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<TV>* levelset=Initialize_Implicit_Surface(*surface,100);
            LOG::cout<<"Seeding particles..."<<std::endl;
            Seed_Particles(*levelset,[=](const TV& X){return TV(0,-3,0);},0,density_sand,particles_per_cell);  // initial velocity
            LOG::cout<<"Particle count: "<<this->particles.number<<std::endl;
            
            // CLASSIFY STRONG AND WEAK SAND
            ARRAY<int> strong_sand;
            ARRAY<int> weak_sand;
            LOG::cout<<"Building kdtree from grain boundary point cloud..."<<std::endl;
            KD_TREE<TV> kdtree;
            ARRAY<TV> point_cloud;
            std::ifstream fs;
            std::string filename=data_directory+"/../Private_Data/voronoi_100seeds_grain_pc.obj";
            fs.open(filename.c_str());
            std::string line;
            while(std::getline(fs,line)){
                std::stringstream ss(line);
                if(line[0]!='v') continue;
                else{ss.ignore();TV x;
                    for(int i=0;i<3;++i) ss>>x(i);
                    x/=3;x+=TV(0,0.15,0);
                    point_cloud.Append(x);}}
            kdtree.Create_Left_Balanced_KD_Tree(point_cloud);
            ARRAY<T> distance2;
            int number_of_points_in_estimate=1;
            for(int i=0;i<particles.number;i++){
                ARRAY<int> points_found(number_of_points_in_estimate);
                ARRAY<T> squared_distance_to_points_found(number_of_points_in_estimate);
                int number_of_points_found;T max_squared_distance_to_points_found;
                kdtree.Locate_Nearest_Neighbors(particles.X(i),FLT_MAX,points_found,
                    squared_distance_to_points_found,number_of_points_found,max_squared_distance_to_points_found,point_cloud);
                distance2.Append(max_squared_distance_to_points_found);
                PHYSBAM_ASSERT(number_of_points_found==number_of_points_in_estimate);}
            T max_distance=sqrt(distance2.Max());
            T min_distance=sqrt(distance2.Min());
            LOG::cout<<"sand to grain boudnary distance measure: "<<std::endl;
            LOG::cout<<"max distance: "<<max_distance<<std::endl;
            LOG::cout<<"min distance: "<<min_distance<<std::endl;
            T critical_distance2=sqr(ratio_of_max_distance_considered_to_be_weak*max_distance);
            for(int i=0;i<particles.number;i++){
                if(distance2(i)<critical_distance2) weak_sand.Append(i);
                else strong_sand.Append(i);}
            LOG::cout<<"# strong sand: "<<strong_sand.m<<std::endl;
            LOG::cout<<"# weak sand: "<<weak_sand.m<<std::endl;
            
            // DELETE WEAK PARTICLES (HOLLOOWWW)
            for(int k=0;k<weak_sand.m;k++){
                int i=weak_sand(k);
                particles.Add_To_Deletion_List(i);}
            particles.Delete_Elements_On_Deletion_List();
            strong_sand.Clean_Memory();
            weak_sand.Clean_Memory();
            for(int k=0;k<particles.number;k++) strong_sand.Append(k);

            ARRAY<int> all_sand;
            for(int i=0;i<particles.number;i++) all_sand.Append(i);
            Add_Drucker_Prager(0,0,(T)35,&all_sand,false,cohesion);

            // GRAVITY
            Add_Gravity(m/(s*s)*TV(0,-9.80665,0));
            // Add_Gravity(m/(s*s)*TV(0,0,-9.80665));


            // // WATER MASS AND VOLUME
            // T volume_water=particles.volume(0)*porosity*saturation_level;
            // T mass_water=density_water*volume_water;

            // // STRONG WATER
            // ARRAY<int> strong_water;
            // T lambda_strong_water=E_strong_water*nu/((1+nu)*(1-2*nu));
            // for(int k=0;k<strong_sand.m;k++){
            //     int i=strong_sand(k);
            //     int p=particles.Add_Element();
            //     particles.mass(p)=mass_water;
            //     strong_water.Append(p);
            //     particles.lambda(p)=lambda_strong_water;
            //     particles.lambda0(p)=lambda_strong_water;
            //     particles.valid(p)=true;
            //     particles.X(p)=particles.X(i);
            //     particles.V(p)=particles.V(i);
            //     particles.F(p)=particles.F(i);
            //     if(particles.store_Fp) particles.Fp(p)=particles.Fp(i); 
            //     if(particles.store_B) particles.B(p)=particles.B(i);
            //     if(particles.store_C) particles.C(p)=particles.C(i);
            //     if(particles.store_S) particles.S(p)=particles.S(i);
            //     particles.volume(p)=volume_water;
            //     particles.mu(p)=(T)0;
            //     particles.mu0(p)=(T)0;}
            // Add_Fixed_Corotated(E_strong_water,nu,&strong_water,true);

            // // WEAK WATER
            // ARRAY<int> weak_water;
            // T lambda_weak_water=E_weak_water*nu/((1+nu)*(1-2*nu));
            // for(int k=0;k<weak_sand.m;k++){
            //     int i=weak_sand(k);
            //     int p=particles.Add_Element();
            //     particles.mass(p)=mass_water;
            //     weak_water.Append(p);
            //     particles.lambda(p)=lambda_weak_water;
            //     particles.lambda0(p)=lambda_weak_water;
            //     particles.valid(p)=true;
            //     particles.X(p)=particles.X(i);
            //     particles.V(p)=particles.V(i);
            //     particles.F(p)=particles.F(i);
            //     if(particles.store_Fp) particles.Fp(p)=particles.Fp(i); 
            //     if(particles.store_B) particles.B(p)=particles.B(i);
            //     if(particles.store_C) particles.C(p)=particles.C(i);
            //     if(particles.store_S) particles.S(p)=particles.S(i);
            //     particles.volume(p)=volume_water;
            //     particles.mu(p)=(T)0;
            //     particles.mu0(p)=(T)0;}
            // Add_Fixed_Corotated(E_weak_water,nu,&weak_water,true);

            // ASSIGN MU AND LAMBDA FOR ALL SAND
            T mu_strong_sand=E_strong_sand/(2*(1+nu));
            T lambda_strong_sand=(E_strong_sand*nu)/((1+nu)*(1-2*nu));
            for(int k=0;k<strong_sand.m;k++){
                int i=strong_sand(k);
                particles.mu(i)=mu_strong_sand;
                particles.mu0(i)=mu_strong_sand;
                particles.lambda(i)=lambda_strong_sand;
                particles.lambda0(i)=lambda_strong_sand;}
            T mu_weak_sand=E_weak_sand/(2*(1+nu));
            T lambda_weak_sand=(E_weak_sand*nu)/((1+nu)*(1-2*nu));
            for(int k=0;k<weak_sand.m;k++){
                int i=weak_sand(k);
                particles.mass(i)*=0.99; // fanfu hack
                particles.mu(i)=mu_weak_sand;
                particles.mu0(i)=mu_weak_sand;
                particles.lambda(i)=lambda_weak_sand;
                particles.lambda0(i)=lambda_weak_sand;}
            this->Update_Variable_Lame_Parameters_On_Constitutive_Models();

        } break;

        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
    if(forced_collision_type!=-1)
        for(int i=0;i<collision_objects.m;i++)
            collision_objects(i)->type=(COLLISION_TYPE)forced_collision_type;
    Update_Variable_Lame_Parameters_On_Constitutive_Models();
}
//#####################################################################
// Function Begin_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Begin_Frame(const int frame)
{
    if(begin_frame) begin_frame(frame);
}
//#####################################################################
// Function End_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
End_Frame(const int frame)
{
    if(end_frame) end_frame(frame);
}
//#####################################################################

// Function Begin_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Begin_Time_Step(const T time)
{
    if(begin_time_step) begin_time_step(time);

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
        this->gather_scatter.template Scatter<int>(false,0,
            [this,&my_particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
            {
                T w=it.Weight();
                TV_INT index=it.Index();
                this->mass(index)+=w*my_particles.mass(p);
            });

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
        TRIANGULATED_SURFACE<T>& new_sc=Seed_Lagrangian_Particles(*surface,0,0,(T)0.001*unit_rho,true);
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
}
//#####################################################################
// Function End_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
End_Time_Step(const T time)
{
    if(end_time_step) end_time_step(time);
}

template class STANDARD_TESTS<VECTOR<float,3> >;
template class STANDARD_TESTS<VECTOR<double,3> >;
}
