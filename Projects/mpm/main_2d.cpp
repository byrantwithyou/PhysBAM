//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/TYPED_STREAM.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include "DELAUNAY_TRIANGULATION_2D.h"
#include "MPM_SIMULATION.h"
#include "MPM_PROJECTION.h"
#include "SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL.h"
#include "SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON.h"
#include "TIMING.h"
#include "VORONOI_2D.h"
#include <omp.h>
using namespace PhysBAM;

template<class TV>
void Run_Simulation(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,2> TV_INT;
    typedef float RW;
    MPM_SIMULATION<TV> sim;
    VORONOI_2D<T> voronoi;
    int test_number=-1;
    std::string output_directory="";
    bool use_output_directory=false;
    sim.dt=(T)1e-3;
    int frame_jump=20;
    int grid_res=16,particle_res=32,particle_count=100;
    T particle_exclude_radius=0;
    T density_scale=(T)1;
    T ym=1;
    T pr=(T).3;
    bool use_voronoi=false;
    bool use_voronoi_boundary=false;
    bool use_turk=false;
    bool use_bridson=false;
    bool use_delaunay=false;
    bool use_projection=false;
    T delaunay_maximum_edge_length=(T)99999;
    T delaunay_minimum_angle=(T)0;
    sim.xi=(T)0;
    sim.visco_kappa=(T)0;
    sim.PROFILING=false;
    parse_args.Add("-test",&test_number,"test","test number");
    parse_args.Add("-o",&output_directory,&use_output_directory,"o","output directory");
    parse_args.Add("-dt",&sim.dt,"dt","dt");
    parse_args.Add("-fj",&frame_jump,"fj","frame jump");
    parse_args.Add("-dump_matrix",&sim.dump_matrix,"dump linear system");
    parse_args.Add("-test_system",&sim.test_system,"test linear system");
    parse_args.Add("-voronoi",&use_voronoi,"use voronoi mesh reconstruction on the fly");
    parse_args.Add("-voronoib",&use_voronoi_boundary,"use voronoi boundary mesh reconstruction on the fly");
    parse_args.Add("-turk",&use_turk,"use turk mesh reconstruction on the fly");
    parse_args.Add("-bridson",&use_bridson,"use bridson mesh reconstruction on the fly");
    parse_args.Add("-delaunay",&use_delaunay,"use delaunay to generate voronoi");
    parse_args.Add("-delaunay_maxl",&delaunay_maximum_edge_length,"value","triangles with edge longer than this will get deleted");
    parse_args.Add("-delaunay_mina",&delaunay_minimum_angle,"value","triangles with angle smaller than this will get deleted");
    parse_args.Add("-stiffness",&ym,"value","scale stiffness");
    parse_args.Add("-poisson_ratio",&pr,"value","poisson's ratio");
    parse_args.Add("-hardening",&sim.xi,"value","harderning coefficient for normal plasticity");
    parse_args.Add("-visco_hardening",&sim.visco_kappa,"value","harderning coefficient for visco plasticity");
    parse_args.Add("-flip",&sim.FLIP_alpha,"value","flip fraction");
    parse_args.Add("-gres",&grid_res,"value","grid resolution");
    parse_args.Add("-pres",&particle_res,"value","particle resolution");
    parse_args.Add("-pn",&particle_count,"value","particle number");
    parse_args.Add("-exclude",&particle_exclude_radius,"value","particle exclude radius when using pn");
    parse_args.Add("-rho",&density_scale,"value","scale object density");
    parse_args.Add("-profile",&sim.PROFILING,"print out timing statements");
    parse_args.Add("-projection",&use_projection,"use poisson projection to enforce incompressibility");
    parse_args.Parse(true);
    parse_args.Parse();

    if(!use_output_directory) output_directory=STRING_UTILITIES::string_sprintf("MPM_%dD_test_%d",TV::m,test_number);

    T object_density=1;
    T object_mass=1;


    // {
    //     HASHTABLE<int,int> HT;
    //     HT.Get_Or_Insert(2)++;
    //     LOG::cout<<*HT.Get_Pointer(2)<<std::endl;

    //     GRID<TV> g(TV_INT(4,3),RANGE<TV>(TV(-3,-2),TV(3,2)));
        
    //     for(int axis=0;axis<2;axis++){
    //         for(int face=0;face<2;face++){
    //             TV_INT node(2,1);
    //             TV_INT face_index=g.Node_Face_Index(axis,node,face);
    //             LOG::cout<<"node "<<node<<std::endl;
    //             LOG::cout<<"axis "<<axis<<" face "<<face<<" face_index "<<face_index<<std::endl;}}

    //     TV_INT nodes[4];
    //     g.Nodes_In_Cell_From_Minimum_Corner_Node(TV_INT(1,0),nodes);
    //     for(int k=0;k<4;k++) LOG::cout<<nodes[k]<<std::endl;
        
    //     for(int axis=0;axis<2;axis++){
    //         for(int face=0;face<2;face++){
    //             TV_INT node(3,1);
    //             TV_INT face_index=g.Node_Face_Index(axis,node,face);
    //             LOG::cout<<"node "<<node<<std::endl;
    //             LOG::cout<<"axis "<<axis<<" face "<<face<<" face_index "<<face_index<<std::endl;}}
    // }

    // {
    //     GRID<TV> g(TV_INT(4,3),RANGE<TV>(TV(-3,-2),TV(3,2)),true);
        
    //     for(int axis=0;axis<2;axis++){
    //         for(int face=0;face<2;face++){
    //             TV_INT node(2,1);
    //             TV_INT face_index=g.Node_Face_Index(axis,node,face);
    //             LOG::cout<<"node "<<node<<std::endl;
    //             LOG::cout<<"axis "<<axis<<" face "<<face<<" face_index "<<face_index<<std::endl;}}

    //     TV_INT nodes[4];
    //     g.Nodes_In_Cell_From_Minimum_Corner_Node(TV_INT(1,0),nodes);
    //     for(int k=0;k<4;k++) LOG::cout<<nodes[k]<<std::endl;
        
    //     for(int axis=0;axis<2;axis++){
    //         for(int face=0;face<2;face++){
    //             TV_INT node(3,1);
    //             TV_INT face_index=g.Node_Face_Index(axis,node,face);
    //             LOG::cout<<"node "<<node<<std::endl;
    //             LOG::cout<<"axis "<<axis<<" face "<<face<<" face_index "<<face_index<<std::endl;}}
    // }


    // {
    //     typedef VECTOR<T,3> TV3;
    //     typedef VECTOR<int,3> TV_INT3;
    
    //     GRID<TV3> g(TV_INT3(3,3,3),RANGE<TV3>(TV3(-3,-3,-3),TV3(3,3,3)));
        
    //     for(int axis=0;axis<3;axis++){
    //         for(int face=0;face<4;face++){
    //             TV_INT3 node(2,1,1);
    //             TV_INT3 face_index=g.Node_Face_Index(axis,node,face);
    //             LOG::cout<<"node "<<node<<std::endl;
    //             LOG::cout<<"axis "<<axis<<" face "<<face<<" face_index "<<face_index<<std::endl;}}
    // }



    // geometry setting
    switch(test_number){
        case 1: // stretching beam
            sim.grid.Initialize(TV_INT(2*grid_res+1,0.5*grid_res+1),RANGE<TV>(TV(-1,-0.25),TV(1,0.25)));
            sim.particles.Initialize_X_As_A_Randomly_Sampled_Box(particle_count,RANGE<TV>(TV(-0.2,-0.12),TV(0.2,0.12)));
            sim.dirichlet_box.Append(RANGE<TV>(TV(0.16,-0.15),TV(0.25,0.15)));
            sim.dirichlet_velocity.Append(TV(0.2,0));
            sim.dirichlet_box.Append(RANGE<TV>(TV(-0.25,-0.15),TV(-0.16,0.15)));
            sim.dirichlet_velocity.Append(TV(-0.2,0));
            break;
        case 2: // shit fall
            sim.grid.Initialize(TV_INT(1.2*grid_res+1,1.2*grid_res+1),RANGE<TV>(TV(-0.6,-0.6),TV(0.6,0.6)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.2*particle_res+1,0.6*particle_res+1),RANGE<TV>(TV(-0.1,-0.5),TV(0.1,0.1)));
            sim.particles.Add_X_As_A_Grid(TV_INT(0.2*particle_res+1,0.2*particle_res+1),RANGE<TV>(TV(-0.1,0.2),TV(0.1,0.4)));
            sim.rigid_ball.Append(SPHERE<TV>(TV(0,-0.2),0.03));
            sim.rigid_ball_velocity.Append(TV());
            sim.rigid_ball.Append(SPHERE<TV>(TV(0.09,-0.3),0.03));
            sim.rigid_ball_velocity.Append(TV());
            break;
        case 3: // snow 
            sim.grid.Initialize(TV_INT(1.0*grid_res+1,0.8*grid_res+1),RANGE<TV>(TV(-0.5,-0.4),TV(0.5,0.4)));
            sim.particles.Initialize_X_As_A_Randomly_Sampled_Box(particle_count,RANGE<TV>(TV(-0.3,-0.3),TV(-0.1,0.2)));
            // sim.particles.Add_X_As_A_Randomly_Sampled_Box(particle_count,RANGE<TV>(TV(0.1,-0.1),TV(0.3,0.1)));
            // for(int p=0;p<sim.particles.number;p++){
            //     if(sim.particles.X(p).x<0) sim.particles.V(p)=TV(2,0);
            //     if(sim.particles.X(p).x>0) sim.particles.V(p)=TV(-2,0);}
            // sim.rigid_ball.Append(SPHERE<TV>(TV(0,-0.2),0.05));
            // sim.rigid_ball_velocity.Append(TV());
            break;
        case 4: // notch
            sim.grid.Initialize(TV_INT(2*grid_res+1,2*grid_res+1),RANGE<TV>(TV(-1,-1),TV(1,1)));
            sim.particles.Initialize_X_As_A_Randomly_Sampled_Box(particle_count,RANGE<TV>(TV(-0.4,-0.2),TV(0.4,0.2)));
            break;
        default: PHYSBAM_FATAL_ERROR("Missing test");};

    // VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),GRID<TV>(sim.grid.numbers_of_cells+1,RANGE<TV>(sim.grid.domain.min_corner-sim.grid.dX*0.5,sim.grid.domain.max_corner+sim.grid.dX*0.5),true),output_directory);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,output_directory);

    // Delaunay Triangulation
    TRIANGULATED_AREA<T> ta;
    if(use_delaunay){
        DELAUNAY_TRIANGULATION_2D<T>::Triangulate(sim.particles.X,ta,delaunay_maximum_edge_length,delaunay_minimum_angle);
        for(int s=0;s<ta.mesh.elements.m;s++){
            Add_Debug_Object(VECTOR<TV,TV::m>(sim.particles.X(ta.mesh.elements(s)(0)),sim.particles.X(ta.mesh.elements(s)(1))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));
            Add_Debug_Object(VECTOR<TV,TV::m>(sim.particles.X(ta.mesh.elements(s)(1)),sim.particles.X(ta.mesh.elements(s)(2))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));
            Add_Debug_Object(VECTOR<TV,TV::m>(sim.particles.X(ta.mesh.elements(s)(2)),sim.particles.X(ta.mesh.elements(s)(0))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("delaunay triangulation");
        SEGMENTED_CURVE_2D<T>& boundary_curve=ta.Get_Boundary_Object();
        for(int s=0;s<boundary_curve.mesh.elements.m;s++) Add_Debug_Object(VECTOR<TV,TV::m>(ta.particles.X.Subset(boundary_curve.mesh.elements(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("delaunay triangulation boundary");}
    
    // voronoi reconstruction
    if(use_voronoi){
        // voronoi.Initialize_With_A_Regular_Grid_Of_Particles(GRID<TV>(TV_INT(0.4*particle_res+1,0.24*particle_res+1),RANGE<TV>(TV(-0.2,-0.12),TV(0.2,0.12))));
        // voronoi.Initialize_With_A_Triangulated_Area(ta);
        voronoi.Initialize_With_And_As_A_Triangulated_Area_And_Relocate_Particles_To_Tri_Centers(ta,sim.particles);
        voronoi.Build_Segments();
        for(int i=0;i<voronoi.X.m;i++) {
            if(voronoi.type(i)==1) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,0,0));
            else if(voronoi.type(i)==10) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,1,0));
            else if(voronoi.type(i)==100) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(0,0,1));}
        for(int s=0;s<voronoi.segments.m;s++){
            Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("voronoi cells");}
    if(use_voronoi_boundary){
        // voronoi.Initialize_With_A_Regular_Grid_Of_Particles(GRID<TV>(TV_INT(0.4*particle_res+1,0.24*particle_res+1),RANGE<TV>(TV(-0.2,-0.12),TV(0.2,0.12))));
        // voronoi.Initialize_With_A_Triangulated_Area(ta);
        voronoi.Initialize_With_And_As_A_Triangulated_Area_And_Relocate_Particles_To_Tri_Centers(ta,sim.particles);
        voronoi.Build_Boundary_Segments();
        for(int i=0;i<voronoi.X.m;i++) {
            if(voronoi.type(i)==1) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,0,0));
            else if(voronoi.type(i)==10) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,1,0));
            else if(voronoi.type(i)==100) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(0,0,1));}
        for(int s=0;s<voronoi.boundary_segments.m;s++){
            Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.boundary_segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("voronoi cell boundary");}

    // material setting
    sim.Initialize();
    if(abs(sim.grid.dX(0)-sim.grid.dX(1))>(T)1e-10){
        LOG::cout<<"grid dx: "<<sim.grid.dX<<std::endl;
        exit(0);}
    switch(test_number){
        case 1:
            object_density=(T)1200*density_scale;
            object_mass=object_density*RANGE<TV>(TV(-0.2,-0.12),TV(0.2,0.12)).Size();
            ym*=(T)5e3;
            for(int p=0;p<sim.particles.number;p++){
                // T this_ym=0.7*ym/0.04*sqr(sim.particles.X(p)(0))+0.3*ym;
                T this_ym=ym;
                sim.mu(p)=(this_ym/((T)2*((T)1+pr)));
                sim.lambda(p)=(this_ym*pr/(((T)1+pr)*((T)1-2*pr)));}
            sim.use_gravity=false;
            sim.use_plasticity_yield=false;
            sim.yield_min=-100;
            sim.yield_max=1.3;
            break;
        case 2:
            object_density=(T)1200*density_scale;
            object_mass=object_density*(RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1))).Size()*2;
            ym*=(T)1e5;
            sim.mu.Fill(ym/((T)2*((T)1+pr)));
            sim.lambda.Fill(ym*pr/(((T)1+pr)*((T)1-2*pr)));
            sim.use_visco_plasticity=true;
            sim.visco_nu=1e3;
            sim.visco_tau.Fill(1000);
            break;
        case 3:
            object_density=(T)400*density_scale;
            object_mass=object_density*(RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1))).Size();
            ym*=1e5;
            for(int p=0;p<sim.particles.number;p++){
                T this_ym=0;
                sim.mu(p)=(this_ym/((T)2*((T)1+pr)));
                sim.lambda(p)=(this_ym*pr/(((T)1+pr)*((T)1-2*pr)));}
            sim.use_gravity=true;
            sim.use_plasticity_yield=false;
            sim.yield_min=(T)1-(T)0.025;
            sim.yield_max=(T)1+(T)0.0075;
            break;
        case 4:
            object_density=(T)1200*density_scale;
            object_mass=object_density*(RANGE<TV>(TV(-0.1,-0.2),TV(0.1,0.2)).Size()-SPHERE<TV>(TV(0.1,0),0.04).Size());
            ym*=(T)5e3;
            for(int p=0;p<sim.particles.number;p++){
                T alpha=1;
                T this_ym=((T)1-alpha)/(T)0.04*ym*sqr(sim.particles.X(p)(1))+alpha*ym;
                sim.mu(p)=(this_ym/((T)2*((T)1+pr)));
                sim.lambda(p)=(this_ym*pr/(((T)1+pr)*((T)1-2*pr)));}
            sim.use_gravity=false;
            break;
        default: PHYSBAM_FATAL_ERROR("Missing test");};
    for(int p=0;p<sim.particles.number;p++){
        sim.particles.mass(p)=object_mass/sim.particles.number;
        sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
        sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
    sim.mu0=sim.mu;
    sim.lambda0=sim.lambda;

    // projection init
    MPM_PROJECTION<TV> projection(sim);

    // use voronoi polygon to initialize particle mass and volume and particle domain
    if(use_voronoi || use_voronoi_boundary){
        sim.assigned_volume_externally=true;
        for(int p=0;p<sim.particles.number;p++){
            POLYGON<TV> poly(voronoi.elements(p).m);
            for(int i=0;i<voronoi.elements(p).m;i++)
                poly.X(i)=voronoi.Xm(voronoi.elements(p)(i));
            T area=poly.Area();
            sim.particles.volume(p)=area;
            sim.particles.mass(p)=object_density*area;
            if(poly.X.m==TV::m+1){ // triangle mesh
                for(int i=0;i<poly.X.m;i++) sim.particles.particle_domain(p)(i)=poly.X(i);
                if(!((poly.X(1)-poly.X(0)).x*(poly.X(2)-poly.X(0)).y-(poly.X(1)-poly.X(0)).y*(poly.X(2)-poly.X(0)).x>0)){
                    TV temp=sim.particles.particle_domain(p)(1);
                    sim.particles.particle_domain(p)(1)=sim.particles.particle_domain(p)(2);
                    sim.particles.particle_domain(p)(2)=temp;}}}}

    // Greg Turk
    if(use_turk){
        SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV> recons;
        ARRAY<TV> xbar;
        ARRAY<MATRIX<T,TV::m> > G;
        ARRAY<T> density;
        GRID<TV> recons_grid(TV_INT(0.6*300+1,2*300+1),RANGE<TV>(TV(-0.3,-1),TV(0.3,1)));
        ARRAY<T,TV_INT> phi;
        recons.Compute_Kernal_Centers_And_Transformation_And_Density(sim.particles.X,sim.particles.mass,0.02,0.04,0.9,10,4,1400,0.5,xbar,G,density);
        recons.Build_Scalar_Field(xbar,sim.particles.mass,density,G,recons_grid,phi);
        T k;std::cin>>k;
        for(int i=0;i<phi.array.m;i++) phi.array(i)-=k;
        // Dump_Levelset(recons_grid,phi,VECTOR<T,3>(0,0,1),VECTOR<T,3>(0,0,1));
    }

    // Zhu and Bridson
    SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<TV> recons_bridson;
    GRID<TV> recons_bridson_grid(TV_INT()+(grid_res*8+1),RANGE<TV>(TV(-1,-1),TV(1,1)));
    ARRAY<T,TV_INT> phi_bridson;
    if(use_bridson){
        recons_bridson.Initialize(particle_exclude_radius*2);
        recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,1);
        phi_bridson.array+=1;}

    // draw wall
    for(T x=sim.grid.domain.min_corner(0)+3.5*sim.grid.dX.Min();x<sim.grid.domain.max_corner(0)-3.5*sim.grid.dX.Min();x+=sim.grid.dX.Min()*0.2){
        Add_Debug_Particle(TV(x,sim.grid.domain.min_corner(1)+3.5*sim.grid.dX.Min()),VECTOR<T,3>(0,0,1));
        Add_Debug_Particle(TV(x,sim.grid.domain.max_corner(1)-3.5*sim.grid.dX.Min()),VECTOR<T,3>(0,0,1));}
    for(T y=sim.grid.domain.min_corner(1)+3.5*sim.grid.dX.Min();y<sim.grid.domain.max_corner(1)-3.5*sim.grid.dX.Min();y+=sim.grid.dX.Min()*0.2){
        Add_Debug_Particle(TV(sim.grid.domain.min_corner(0)+3.5*sim.grid.dX.Min(),y),VECTOR<T,3>(0,0,1));
        Add_Debug_Particle(TV(sim.grid.domain.max_corner(0)-3.5*sim.grid.dX.Min(),y),VECTOR<T,3>(0,0,1));}

    // draw collision objects
    for(int b=0;b<sim.rigid_ball.m;b++){
        for(int k=0;k<50;k++){
            T theta=k*2.0*3.14/50.0;
            TV disp;disp(0)=sim.rigid_ball(b).radius*cos(theta);disp(1)=sim.rigid_ball(b).radius*sin(theta);
            Add_Debug_Particle(sim.rigid_ball(b).center+disp,VECTOR<T,3>(0,0,1));}}

    // draw MPM particles
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,1,1));

    Flush_Frame<TV>("mpm");

    for(int f=1;f<2900977;f++){
        TIMING_START;
        LOG::cout<<"MPM TIMESTEP "<<f<<std::endl;

        sim.Build_Weights_And_Grad_Weights();
        sim.Build_Helper_Structures_For_Constitutive_Model();
        sim.Rasterize_Particle_Data_To_The_Grid();
        if(sim.frame==0 && !sim.assigned_volume_externally) sim.Compute_Particle_Volumes_And_Densities();
        sim.Compute_Grid_Forces();
        if(sim.use_gravity) sim.Apply_Gravity_To_Grid_Forces();
        sim.Update_Velocities_On_Grid();
        // sim.Grid_Based_Body_Collisions();
        sim.Solve_The_Linear_System(); // so far sim.node_V is achieved via MPM
        if(use_projection){
            projection.Reinitialize();
            projection.Identify_Dirichlet_Cells();
            projection.Identify_Neumann_Cells();
            projection.Velocities_Corners_To_faces();
            projection.Build_Velocity_Divergence();
            projection.Solve_For_Pressure(sim.dt,1);
            projection.Do_Projection(sim.dt,1);
            projection.Velocities_Faces_To_Corners(); // this step modifies sim.node_V
        }
        sim.Update_Deformation_Gradient();
        sim.Update_Particle_Velocities();
        sim.Particle_Based_Body_Collisions();
        sim.Update_Particle_Positions();
        sim.Update_Dirichlet_Box_Positions();
        sim.Update_Colliding_Object_Positions();
        sim.frame++;

        TIMING_END("Current MPM time step totally");

        if(f%frame_jump==0){
            // draw MPM particles
            for(int i=0;i<sim.particles.X.m;i++){
                if(!sim.failed(i) && sim.valid(i)) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,1,1));
                else if(sim.failed(i) && sim.valid(i)) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,0,1));
                else Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,0,0));}

            // projection: visualize MAC grid velocities
            // if(use_projection){
            //     for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),projection.mac_grid.counts));it.Valid();it.Next()){
            //         Add_Debug_Particle(projection.mac_grid.X(it.index),VECTOR<T,3>(0,0,0)); // cell centers: red
            //         Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.node_V(it.index));}
            //     for(FACE_ITERATOR<TV> iterator(projection.mac_grid);iterator.Valid();iterator.Next()){
            //         TV location=iterator.Location();
            //         int axis=iterator.Axis();
            //         if(axis==0){
            //             Add_Debug_Particle((location),VECTOR<T,3>(0,0,0)); // face centers with x component velocity: green
            //             Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(projection.face_velocities(iterator.Full_Index()),0));}
            //         else if(axis==1){
            //             Add_Debug_Particle((location),VECTOR<T,3>(0,0,0)); // face centers with y component velocity: blue
            //             Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(0,projection.face_velocities(iterator.Full_Index())));}}}

            // Zhu and Bridson
            if(use_bridson){
                recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,1);
                Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(1,0,0));
                LEVELSET<TV> ls(recons_bridson_grid,phi_bridson,0);
                ls.Fast_Marching_Method();
                phi_bridson.array+=particle_exclude_radius*3;
                ls.Fast_Marching_Method();
                phi_bridson.array-=particle_exclude_radius*2;
                ls.Fast_Marching_Method();
                Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(0,1,0));}

            // voronoi reconstruction
            if(use_voronoi){
                voronoi.Crack(sim.particles.X,sim.grid.dX.Min()*1.8);
                voronoi.Build_Association();
                voronoi.Build_Segments();
                voronoi.Deform_Mesh_Using_Particle_Deformation(sim.particles.Xm,sim.particles.X,sim.particles.Fe,sim.particles.Fp,true,3);
                for(int s=0;s<voronoi.segments.m;s++){
                    Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}}
            if(use_voronoi_boundary){
                voronoi.Crack(sim.particles.X,sim.grid.dX.Min()*1.3);
                voronoi.Build_Association();
                voronoi.Build_Boundary_Segments();
                voronoi.Deform_Mesh_Using_Particle_Deformation(sim.particles.Xm,sim.particles.X,sim.particles.Fe,sim.particles.Fp,true,1);
                for(int s=0;s<voronoi.boundary_segments.m;s++){
                    Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.boundary_segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}}

            // draw wall
            for(T x=sim.grid.domain.min_corner(0)+4.5*sim.grid.dX.Min();x<sim.grid.domain.max_corner(0)-4.5*sim.grid.dX.Min();x+=sim.grid.dX.Min()*0.2){
                Add_Debug_Particle(TV(x,sim.grid.domain.min_corner(1)+4.5*sim.grid.dX.Min()),VECTOR<T,3>(0,0,1));
                Add_Debug_Particle(TV(x,sim.grid.domain.max_corner(1)-4.5*sim.grid.dX.Min()),VECTOR<T,3>(0,0,1));}
            for(T y=sim.grid.domain.min_corner(1)+4.5*sim.grid.dX.Min();y<sim.grid.domain.max_corner(1)-4.5*sim.grid.dX.Min();y+=sim.grid.dX.Min()*0.2){
                Add_Debug_Particle(TV(sim.grid.domain.min_corner(0)+4.5*sim.grid.dX.Min(),y),VECTOR<T,3>(0,0,1));
                Add_Debug_Particle(TV(sim.grid.domain.max_corner(0)-4.5*sim.grid.dX.Min(),y),VECTOR<T,3>(0,0,1));}
            
            // draw collision objects
            for(int b=0;b<sim.rigid_ball.m;b++){
                for(int k=0;k<50;k++){
                    T theta=k*2.0*3.14/50.0;
                    TV disp;disp(0)=sim.rigid_ball(b).radius*cos(theta);disp(1)=sim.rigid_ball(b).radius*sin(theta);
                    Add_Debug_Particle(sim.rigid_ball(b).center+disp,VECTOR<T,3>(0,0,1));}}

            Flush_Frame<TV>("mpm");}
        LOG::cout<<std::endl;}
}

int main(int argc,char *argv[])
{
    // ./mpm_2d -test 3 -gres 160 -pn 5000 -dt 1e-4 -fj 80 -hardening 10

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse(true);
    Run_Simulation<VECTOR<double,2> >(parse_args);
    LOG::Finish_Logging();
    return 0;
}


