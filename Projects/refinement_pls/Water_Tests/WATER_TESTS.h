//#####################################################################
// Copyright 2009, Michael Lentine, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_TESTS__
#define __SMOKE_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Dynamics/PLS_REFINEMENT_EXAMPLE.h>
#include "REFINEMENT_THREADS.h"

namespace PhysBAM{

template<class TV>
class WATER_TESTS:public PLS_REFINEMENT_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef VECTOR<int,TV::dimension-1> TV_INT2;
    typedef PLS_REFINEMENT_EXAMPLE<TV> BASE;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;    

public:
    using BASE::stream_type;using BASE::coarse_mac_grid;using BASE::fine_mac_grid;using BASE::coarse_face_velocities;using BASE::fine_face_velocities;using BASE::fine_mpi_grid;using BASE::coarse_mpi_grid;
    using BASE::incompressible;using BASE::projection;using BASE::domain_boundary;using BASE::use_collidable_advection;using BASE::collision_bodies_affecting_fluid;
    using BASE::particle_levelset_evolution;using BASE::output_directory;using BASE::write_substeps_level;using BASE::last_frame;using BASE::rigid_geometry_collection;using BASE::split_dir;using BASE::write_debug_data;
    using BASE::coarse_phi;using BASE::restart;using BASE::boundary;using BASE::advection_scalar;using BASE::non_mpi_boundary;using BASE::phi_boundary;using BASE::gravity;using BASE::cfl;
    
    int sub_scale,scale;
    T sub_scale_face_inverse; //tmp data to make more efficient operations    
    int binary_refinement_levels;
    ARRAY<T,TV_INT2> weights;
    T alpha;
    ARRAY<T,FACE_INDEX<TV::dimension> > fine_face_velocities_save,coarse_face_velocities_save;
    ARRAY<bool,FACE_INDEX<TV::dimension> > fine_psi_N;
    bool use_direct_solver;
    int number_of_threads_per_mpi_proc;
    int buffer;
    bool use_cubic_interpolation;
    bool surface_solve;
    TV ratio;
    THREAD_QUEUE* thread_queue;

    GRID<TV> local_mac_grid;
    ARRAY<T,FACE_INDEX<TV::dimension> > local_face_velocities;
    PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > local_projection,local_projection_analytic;
    PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > levelset_projection;
    //INCOMPRESSIBLE_UNIFORM<GRID<TV> > levelset_incompressible;

    ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> >*> local_face_velocities_array,local_face_velocities_save_array;
    ARRAY<PROJECTION_UNIFORM<GRID<TV> >*> local_projection_array;

    SPARSE_MATRIX_FLAT_NXN<T> A;
    VECTOR_ND<T> b;
    ARRAY<int,TV_INT> cell_index_to_matrix_index;
    ARRAY<TV_INT,int> matrix_index_to_cell_index;

    ARRAY<T,TV_INT> local_phi;
    ARRAY<T,FACE_INDEX<TV::dimension> > beta_face;    
    
    INTERPOLATION_UNIFORM<GRID<VECTOR<T,TV::dimension-1> >,T>* interpolation;
    INTERPOLATION_UNIFORM<GRID<TV>,T>* phi_interpolation;
    VECTOR<GRID<VECTOR<T,TV::dimension-1> >,TV::dimension> lower_dim_grids;
    RANGE<TV_INT> lower_dim_domain;
    VECTOR<ARRAY<T,VECTOR<int,TV::dimension-1> >,TV::dimension> lower_dim_velocities;
            
    CYLINDER<T> source_cyl;
    RANGE<TV> source;
    
    WATER_TESTS(const STREAM_TYPE stream_type,const PARSE_ARGS& parse_args)
        :PLS_REFINEMENT_EXAMPLE<TV>(stream_type),binary_refinement_levels(0),use_direct_solver(false),buffer(1),use_cubic_interpolation(false),surface_solve(true),thread_queue(0),
        local_mac_grid(TV_INT::All_Ones_Vector(),RANGE<TV>::Centered_Box(),true),local_face_velocities(local_mac_grid),local_projection(local_mac_grid),local_projection_analytic(local_mac_grid),levelset_projection(fine_mac_grid),
        //levelset_incompressible(fine_mac_grid,levelset_projection),
        interpolation(0),phi_interpolation(0)
    {
        int test_number=2;last_frame=200;
        restart=parse_args.Get_Integer_Value("-restart");
        alpha=parse_args.Get_Double_Value("-alpha");
        write_substeps_level=parse_args.Get_Integer_Value("-substep");
        number_of_threads_per_mpi_proc=parse_args.Get_Integer_Value("-threads");
        scale=parse_args.Get_Integer_Value("-scale");
        cfl=parse_args.Get_Double_Value("-cfl");
        binary_refinement_levels=parse_args.Get_Integer_Value("-binary");
        split_dir=parse_args.Get_String_Value("-split");
        if(binary_refinement_levels){
            sub_scale=1;
            for(int i=0;i<binary_refinement_levels;i++) sub_scale<<=1;
            Initialize_Binary_Local_Solve();}
        else sub_scale=parse_args.Get_Integer_Value("-subscale");
        if(sub_scale==1) LOG::cout<<"WARNING: Please use one_phase_incompressible instead"<<std::endl;
        buffer=parse_args.Get_Integer_Value("-buffer");
        use_cubic_interpolation=parse_args.Get_Option_Value("-cubic");
        write_debug_data=parse_args.Get_Option_Value("-write_debug");
        surface_solve=!parse_args.Get_Option_Value("-nosurface");
        output_directory=STRING_UTILITIES::string_sprintf("Water_Tests/Test_%d_%d_%d%s_%1.2f_%d%s%s",test_number,scale,sub_scale,binary_refinement_levels?"_binary":"",alpha,buffer,use_cubic_interpolation?"_cubic":"",surface_solve?"_surface":"");
        ratio=TV::All_Ones_Vector();
        if(parse_args.Get_Double_Value("-x")!=1) ratio(1)*=parse_args.Get_Double_Value("-x");
        if(parse_args.Get_Double_Value("-y")!=1) ratio(2)*=parse_args.Get_Double_Value("-y");
        if(parse_args.Get_Double_Value("-z")!=1) ratio(3)*=parse_args.Get_Double_Value("-z");
        TV_INT dimensions=TV_INT(ratio*scale);
        fine_mac_grid.Initialize(dimensions,RANGE<TV>(TV(),ratio),true);
        coarse_mac_grid.Initialize(dimensions/sub_scale,RANGE<TV>(TV(),ratio),true);
        local_mac_grid.Initialize(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>::Centered_Box(),true);
        local_projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
        local_projection.elliptic_solver->pcg.Use_Modified_Incomplete_Cholesky();
        //levelset_incompressible.Set_Custom_Advection(advection_scalar);
        local_face_velocities.Resize(local_mac_grid);
        local_projection.elliptic_solver->Set_Relative_Tolerance(1e-7);
        local_projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
        local_projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
        local_projection.elliptic_solver->pcg.cg_restart_iterations=0;
        local_projection.elliptic_solver->Solve_Neumann_Regions(true);
        local_projection.Initialize_Grid(local_mac_grid);
        local_projection_analytic.Initialize_Grid(local_mac_grid);
        local_projection.collidable_solver->Use_External_Level_Set(*new T_LEVELSET(local_mac_grid,local_phi));
        if(!use_direct_solver) Initialize_Fast_Unwrapped_PCG_Solve(local_projection);
        weights.Resize(RANGE<TV_INT2>(TV_INT2::All_Ones_Vector(),TV_INT2::All_Ones_Vector()*sub_scale));
        Initialize_Weights(weights,sub_scale);
        PHYSBAM_ASSERT(scale%sub_scale==0 && binary_refinement_levels>=0);
        sub_scale_face_inverse=(TV::dimension==2)?(1./sub_scale):(1./(sub_scale*sub_scale));
        if(use_cubic_interpolation) interpolation=new CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<T,TV::dimension-1> >,T>();
        else interpolation=new LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,TV::dimension-1> >,T>();
        phi_interpolation=new LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T>();
        if(TV::dimension==3){
            VECTOR<T,3> point1,point2;
            point1=VECTOR<T,3>::All_Ones_Vector()*(T).6;point1(1)=.4;point1(2)=.95;point2=VECTOR<T,3>::All_Ones_Vector()*(T).6;point2(1)=.4;point2(2)=1;
            source_cyl.Set_Endpoints(point1,point2);source_cyl.radius=.1;}
        else{
            TV point1,point2;
            point1=TV::All_Ones_Vector()*(T).5;point1(1)=.4;point1(2)=.95;point2=TV::All_Ones_Vector()*(T).65;point2(1)=.55;point2(2)=1;
            source.min_corner=point1;source.max_corner=point2;}
        //gravity=2;
    }

    ~WATER_TESTS()
    {
        if(thread_queue) delete thread_queue;
        if(interpolation) delete interpolation;
        if(phi_interpolation) delete phi_interpolation;
        for(int i=0;i<binary_refinement_levels;i++){
            delete local_face_velocities_array(i);
            delete local_face_velocities_save_array(i);
            delete local_projection_array(i);}
    }

    void Initialize_MPI()
    {
        coarse_face_velocities_save.Resize(coarse_mac_grid);
        fine_face_velocities_save.Resize(fine_mac_grid);
        beta_face.Resize(coarse_mac_grid);
        local_phi.Resize(local_mac_grid.Domain_Indices(1));
        fine_psi_N.Resize(fine_mac_grid);
        for(int i=0;i<TV::dimension;i++){lower_dim_grids(i)=coarse_mac_grid.Remove_Dimension(i);lower_dim_velocities(i).Resize(lower_dim_grids(i).Domain_Indices(3));}
        lower_dim_domain=coarse_mac_grid.Domain_Indices(3);
        
        //levelset_incompressible.mpi_grid=fine_mpi_grid;
        levelset_projection.elliptic_solver->mpi_grid=fine_mpi_grid;
        //levelset_incompressible.Initialize_Grids(fine_mac_grid);
        //levelset_incompressible.Set_Custom_Boundary(*boundary);
        levelset_projection.elliptic_solver->Set_Relative_Tolerance(1e-7);
        levelset_projection.elliptic_solver->pcg.Set_Maximum_Iterations(400);
        levelset_projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
        levelset_projection.elliptic_solver->pcg.cg_restart_iterations=40;
        levelset_projection.elliptic_solver->pcg.Show_Results();
        levelset_projection.elliptic_solver->Solve_Neumann_Regions(false);
        levelset_projection.Initialize_Grid(fine_mac_grid);
        levelset_projection.collidable_solver->Use_External_Level_Set(particle_levelset_evolution.particle_levelset.levelset);

        /*levelset_incompressible.projection.Use_Non_Zero_Divergence(false);
        levelset_incompressible.projection.elliptic_solver->Solve_Neumann_Regions(false);
        levelset_incompressible.projection.elliptic_solver->solve_single_cell_neumann_regions=false;
        levelset_incompressible.Use_Explicit_Part_Of_Implicit_Viscosity(false);
        levelset_incompressible.Set_Maximum_Implicit_Viscosity_Iterations(40);
        levelset_incompressible.Use_Variable_Vorticity_Confinement(false);
        levelset_incompressible.Set_Surface_Tension(0);
        levelset_incompressible.Set_Variable_Surface_Tension(false);
        levelset_incompressible.Set_Viscosity(0);
        levelset_incompressible.Set_Variable_Viscosity(false);
        levelset_incompressible.projection.Set_Density(1e3);*/
    }
        
    void Initialize_Weights(ARRAY<T,VECTOR<int,0> >& weights,const int sub_scale)
    {PHYSBAM_NOT_IMPLEMENTED();}
    
    void Initialize_Weights(ARRAY<T,VECTOR<int,1> >& weights,const int sub_scale)
    {for(int i=0;i<sub_scale;i++) weights(TV_INT2(i))=1/(T)sub_scale;}
    
    void Initialize_Weights(ARRAY<T,VECTOR<int,2> >& weights,const int sub_scale)
    {for(int i=0;i<sub_scale;i++) for(int j=0;j<sub_scale;j++) weights(TV_INT2(i,j))=1/(T)(sub_scale*sub_scale);}

    void Write_Output_Files(const int frame)
    {BASE::Write_Output_Files(frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",levelset_projection.p);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",levelset_projection.elliptic_solver->psi_N);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",levelset_projection.elliptic_solver->psi_D);}

    void Set_Boundary_Conditions_For_Coarse_Only(const T time)
    {projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*(axis-1)+axis_side;
        TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset=axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(domain_boundary(axis)(axis_side)){
            for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT face=iterator.Face_Index()+boundary_face_offset;
                TV local_offset=TV::All_Ones_Vector()*.5*coarse_mac_grid.DX();
                GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>(iterator.Location()-local_offset,iterator.Location()+local_offset),true);
                bool inside=false;
                for(int i=0;i<sub_scale;i++){TV_INT offset=TV_INT::All_Ones_Vector()*(i-1);offset(iterator.Axis())=0;
                    TV_INT coarse_face=face+interior_cell_offset;
                    TV_INT fine_face=sub_scale*(coarse_face-TV_INT::All_Ones_Vector())+TV_INT::All_Ones_Vector()+offset;
                    if(particle_levelset_evolution.particle_levelset.levelset.phi(fine_face)<=0){inside=true;break;}}
                if(inside){
                    if(coarse_face_velocities.Component(axis).Valid_Index(face)){
                        projection.elliptic_solver->psi_N.Component(axis)(face)=true;coarse_face_velocities.Component(axis)(face)=0;}}
                else{TV_INT cell=face+exterior_cell_offset;projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}}
        else for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
            projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
        TV local_offset=TV::All_Ones_Vector()*.5*coarse_mac_grid.DX();
        GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>(iterator.Location()-local_offset,iterator.Location()+local_offset),true);
        bool inside=false;
        for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
            TV_INT fine_face=iterator.Cell_Index()*sub_scale+local_iterator.Cell_Index()-sub_scale*TV_INT::All_Ones_Vector();
            if(particle_levelset_evolution.particle_levelset.levelset.phi(fine_face)<=0){inside=true;break;}}
        if(!inside){
            projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=0;}}}

    void Set_Coarse_Boundary_Conditions(const T time)
    {for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()) beta_face.Component(iterator.Axis())(iterator.Face_Index())=0;
    TV local_offset=TV::All_Ones_Vector()*.5*coarse_mac_grid.DX();
    GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>(-local_offset,local_offset),true);
    FACE_INDEX<TV::dimension> fine_index;beta_face.Fill(0);
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
            fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(iterator.Cell_Index()-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
            if(!fine_psi_N(fine_index)){
                if(local_iterator.First_Boundary()) beta_face(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.First_Face_Index(local_iterator.Axis())))+=T(1);
                else beta_face(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Second_Face_Index(local_iterator.Axis())))+=T(1);}}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
        if(beta_face.Component(iterator.Axis())(iterator.Face_Index())==0) projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
        else{
            int factor=0;
            if(iterator.First_Cell_Index()(iterator.Axis())>0){factor++;}
            if(iterator.Second_Cell_Index()(iterator.Axis())<=coarse_mac_grid.Counts()(iterator.Axis())){factor++;}
            beta_face.Component(iterator.Axis())(iterator.Face_Index())*=(factor==1)?sub_scale_face_inverse:0.5*sub_scale_face_inverse;}}
    projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    int ghost=max(3,sub_scale);
    ARRAY<T,TV_INT> phi_ghost(fine_mac_grid.Domain_Indices(ghost));
    phi_boundary->Fill_Ghost_Cells(fine_mac_grid,particle_levelset_evolution.phi,phi_ghost,0,time,ghost);
    for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*(axis-1)+axis_side;
        TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset=axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(domain_boundary(axis)(axis_side)){
            for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT face=iterator.Face_Index()+boundary_face_offset;
                //if(coarse_phi(face+interior_cell_offset)<=0){
                if(!Contains_Outside(face+interior_cell_offset,phi_ghost,0)){
                    if(coarse_face_velocities.Component(axis).Valid_Index(face)){
                            projection.elliptic_solver->psi_N.Component(axis)(face)=true;coarse_face_velocities.Component(axis)(face)=0;}}
                else{TV_INT cell=face+exterior_cell_offset;projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}}
        else for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
            projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()) if(Contains_Outside(iterator.Cell_Index(),particle_levelset_evolution.phi,0)) {
        projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=0;}
    if(coarse_mpi_grid){
        coarse_mpi_grid->Exchange_Boundary_Cell_Data(projection.elliptic_solver->psi_D,1,false);
        coarse_mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
    for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
        if(time<=3 && Lazy_Inside_Source(iterator.Location())){
            projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==2) coarse_face_velocities(iterator.Full_Index())=-1;
            else coarse_face_velocities(iterator.Full_Index())=0;}}}

    void Set_Coarse_Boundary_Conditions_Only()
    {
        ARRAY<bool,TV_INT>& psi_D=projection.elliptic_solver->psi_D;ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N=projection.elliptic_solver->psi_N;
        for(int axis=0;axis<GRID<TV>::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){
            int side=2*(axis-1)+axis_side;
            TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
            TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
            TV_INT boundary_face_offset=axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
            if(domain_boundary(axis)(axis_side)){
                for(typename GRID<TV>::FACE_ITERATOR iterator(projection.elliptic_solver->grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                    TV_INT face=iterator.Face_Index()+boundary_face_offset;
                    if(particle_levelset_evolution.particle_levelset.levelset.phi(face+interior_cell_offset)<=0){
                        if(coarse_face_velocities.Component(axis).Valid_Index(face)){
                            psi_N.Component(axis)(face)=true;coarse_face_velocities.Component(axis)(face)=0;}}
                    else{TV_INT cell=face+exterior_cell_offset;
                        psi_D(cell)=true;projection.elliptic_solver->u(cell)=0;}}}
            else
                for(typename GRID<TV>::FACE_ITERATOR iterator(projection.elliptic_solver->grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                    psi_D(cell)=true;projection.elliptic_solver->u(cell)=0;}}
    }

    void Set_Fine_Boundary_Conditions(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,const T time)
    {for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*(axis-1)+axis_side;
        TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset=axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(domain_boundary(axis)(axis_side)){
            for(typename GRID<TV>::FACE_ITERATOR iterator(local_mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT face=iterator.Face_Index()+boundary_face_offset;
                if(particle_levelset_evolution.particle_levelset.levelset.phi(face+interior_cell_offset)<=0){
                    if(local_face_velocities.Component(axis).Valid_Index(face)){
                        psi_N.Component(axis)(face)=true;local_face_velocities.Component(axis)(face)=0;}}}}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(local_mac_grid);iterator.Valid();iterator.Next()){
        if(time<=3 && Lazy_Inside_Source(iterator.Location())){
            psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==2) local_face_velocities(iterator.Full_Index())=-1;
            else local_face_velocities(iterator.Full_Index())=0;}}}

    bool Set_Local_Boundary_Conditions(GRID<TV>& local_mac_grid,PROJECTION_UNIFORM<GRID<TV> >& projection,TV_INT coarse_index)
    {projection.elliptic_solver->psi_D.Fill(false);bool contains_surface=false;
    for(typename GRID<TV>::CELL_ITERATOR iterator(local_mac_grid);iterator.Valid();iterator.Next()){
        TV_INT fine_cell=coarse_index*sub_scale+iterator.Cell_Index()-sub_scale*TV_INT::All_Ones_Vector();
        if(particle_levelset_evolution.particle_levelset.levelset.phi(fine_cell)>0){contains_surface=true;projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=0;}}
    return contains_surface;}
        
    void Set_Boundary_Conditions(const T time)
    {fine_psi_N.Fill(false);projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    Set_Fine_Boundary_Conditions(fine_mac_grid,fine_face_velocities,fine_psi_N,time);
    Set_Coarse_Boundary_Conditions(time);
    //projection.poisson->beta_face.Fill(1);
    }//projection.poisson->beta_given_on_faces=true;}

    bool All_Cells_Inside(TV_INT cell_index,int ring_number,const ARRAY<T,TV_INT>& levelset_phi)
    {
        RANGE<TV_INT> ring_domain;ring_domain.min_corner=cell_index-TV_INT::All_Ones_Vector()*ring_number;ring_domain.max_corner=cell_index+TV_INT::All_Ones_Vector()*ring_number;
        coarse_mac_grid.Clamp(ring_domain.min_corner);coarse_mac_grid.Clamp(ring_domain.max_corner);
        for(typename GRID<TV>::CELL_ITERATOR ring_iterator(coarse_mac_grid,ring_domain);ring_iterator.Valid();ring_iterator.Next()) if(Contains_Outside(ring_iterator.Cell_Index(),levelset_phi)) return false;
        //for(typename GRID<TV>::CELL_ITERATOR ring_iterator(coarse_mac_grid,ring_domain);ring_iterator.Valid();ring_iterator.Next()) if(ring_iterator.Location()(2)>.35) return false;
        return true;
    }

    bool Contains_Outside(TV_INT cell_index,const ARRAY<T,TV_INT>& levelset_phi,int buffer)
    { 
        RANGE<TV_INT> domain=local_mac_grid.Domain_Indices();domain.max_corner+=buffer*TV_INT::All_Ones_Vector();domain.min_corner-=buffer*TV_INT::All_Ones_Vector();
        for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid,domain);local_iterator.Valid();local_iterator.Next()){
            TV_INT fine_face=(cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Cell_Index(); //Assume fine_grid is levelset_grid for now
            if(levelset_phi(fine_face)>0) return true;}
        return false;
    }

    bool Contains_Inside(TV_INT cell_index,const ARRAY<T,TV_INT>& levelset_phi,int buffer)
    { 
        RANGE<TV_INT> domain=local_mac_grid.Domain_Indices();domain.max_corner+=buffer*TV_INT::All_Ones_Vector();domain.min_corner-=buffer*TV_INT::All_Ones_Vector();
        for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid,domain);local_iterator.Valid();local_iterator.Next()){
            TV_INT fine_face=(cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Cell_Index(); //Assume fine_grid is levelset_grid for now
            if(levelset_phi(fine_face)<=0) return true;}
        return false;
    }

    void Set_Levelset_Boundary_Conditions(const GRID<TV>& levelset_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& levelset_velocities,const ARRAY<T,TV_INT>& levelset_phi,const T time)
    {levelset_projection.elliptic_solver->psi_D.Fill(false);levelset_projection.elliptic_solver->psi_N.Fill(false);
    for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*(axis-1)+axis_side;
        TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset=axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(domain_boundary(axis)(axis_side)){
            for(typename GRID<TV>::FACE_ITERATOR iterator(levelset_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT face=iterator.Face_Index()+boundary_face_offset;
                if(levelset_phi(face+interior_cell_offset)<=0){
                    if(levelset_velocities.Component(axis).Valid_Index(face)){
                        levelset_projection.elliptic_solver->psi_N.Component(axis)(face)=true;levelset_velocities.Component(axis)(face)=0;}}
                else{TV_INT cell=face+exterior_cell_offset;levelset_projection.elliptic_solver->psi_D(cell)=true;levelset_projection.p(cell)=0;}}}
        else for(typename GRID<TV>::FACE_ITERATOR iterator(levelset_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
            levelset_projection.elliptic_solver->psi_D(cell)=true;levelset_projection.p(cell)=0;}}
    for(typename GRID<TV>::CELL_ITERATOR iterator(levelset_grid);iterator.Valid();iterator.Next()){
        if(levelset_phi(iterator.Cell_Index())>0){
            levelset_projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;levelset_projection.p(iterator.Cell_Index())=0;}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()){
        if(time<=3 && source.Lazy_Inside(iterator.Location())){
            levelset_projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==2) levelset_velocities(iterator.Full_Index())=-1;
            else levelset_velocities(iterator.Full_Index())=0;}}
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){ 
        if(!Contains_Outside(iterator.Cell_Index(),levelset_phi,buffer)) for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){
            TV_INT offset=(axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT::Axis_Vector(axis));
            TV_INT face=iterator.Cell_Index()+(axis_side==1?TV_INT():TV_INT::Axis_Vector(axis));
            if(non_mpi_boundary(axis)(axis_side) && ((iterator.Cell_Index()(axis)+offset(axis))<1 || (iterator.Cell_Index()(axis)+offset(axis))>coarse_mac_grid.Counts()(axis))) continue;
            TV_INT adjacent_cell=iterator.Cell_Index()+offset;
            bool adjacent_outside=false;
            if(Contains_Outside(adjacent_cell,levelset_phi,buffer)) adjacent_outside=true;
            RANGE<TV_INT> domain;domain.min_corner=TV_INT::All_Ones_Vector();domain.max_corner=TV_INT::All_Ones_Vector()*sub_scale;domain.max_corner(axis)=1;
            if(adjacent_outside) for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid,domain);local_iterator.Valid();local_iterator.Next()){
                TV_INT fine_face=sub_scale*(face-TV_INT::All_Ones_Vector())+local_iterator.Cell_Index();
                levelset_projection.elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(axis,fine_face))=true;}}}}
 
    bool Lazy_Inside_Source(const VECTOR<T,2> X)
    {
        return source.Lazy_Inside(X);
    }

    bool Lazy_Inside_Source(const VECTOR<T,3> X)
    {
        return source_cyl.Lazy_Inside(X);
    }

    void Set_Phi_Inside_Source(const VECTOR<int,2>& index,const VECTOR<T,2> X)
    {
        if(source.Lazy_Inside(X)) particle_levelset_evolution.phi(index)=min(particle_levelset_evolution.phi(index),source.Signed_Distance(X));
    }

    void Set_Phi_Inside_Source(const VECTOR<int,3>& index,const VECTOR<T,3> X)
    {
        if(source_cyl.Lazy_Inside(X)) particle_levelset_evolution.phi(index)=min(particle_levelset_evolution.phi(index),source_cyl.Signed_Distance(X));
    }

    void Adjust_Phi_With_Sources(const T time)
    {
        if(time>3) return;
        for(typename GRID<TV>::CELL_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()){Set_Phi_Inside_Source(iterator.Cell_Index(),iterator.Location());}
    }

    void Initialize_Phi()
    {
        ARRAY<T,TV_INT>& phi=particle_levelset_evolution.phi;
        for(typename GRID<TV>::CELL_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()){
            //TV center=ratio*.5;center(2)=.75*ratio(2);
            //static SPHERE<TV> circle(center,(T).2*ratio.Min());
            const TV &X=iterator.Location();
            //phi(iterator.Cell_Index())=min(circle.Signed_Distance(X),X.y-(T).412134*ratio(2));}
            //phi(iterator.Cell_Index())=X.y-(T).412134;}
            //phi(iterator.Cell_Index())=X.y-(T).442134;}
            phi(iterator.Cell_Index())=X.y-(T)fine_mac_grid.min_dX*5;}
    }

    //TODO: Speed these up if needed
    void Adjust_Phi_With_Objects(const T time)
    {
        if(use_collidable_advection) return;
        T tolerance=(T)9.8/24; // dt*gravity where dt=1/24 is based on the length of a frame
        for(int id=0;id<rigid_geometry_collection.particles.Size();id++){
            for(typename GRID<TV>::CELL_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Cell_Index();TV location=fine_mac_grid.X(index);
                if(particle_levelset_evolution.phi(index)<0 && rigid_geometry_collection.Rigid_Geometry(id).Implicit_Geometry_Extended_Value(location)<0){
                    TV V_fluid;
                    for(int i=0;i<TV::dimension;i++) V_fluid(i)=(fine_face_velocities(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i)))+fine_face_velocities(FACE_INDEX<TV::dimension>(i,iterator.Second_Face_Index(i))))/2.;
                    TV V_object=rigid_geometry_collection.Rigid_Geometry(id).Pointwise_Object_Velocity(location); // velocity object should be spatially varying
                    TV V_relative=V_fluid-V_object;
                    TV normal=rigid_geometry_collection.Rigid_Geometry(id).Implicit_Geometry_Normal(location);
                    T VN=TV::Dot_Product(V_relative,normal),magnitude=V_relative.Magnitude();
                    if(VN > max(tolerance,(T).1*magnitude)) particle_levelset_evolution.phi(index)=-rigid_geometry_collection.Rigid_Geometry(id).Implicit_Geometry_Extended_Value(location);}}}
    }

    void Extrapolate_Phi_Into_Objects(const T time)
    {
        if(use_collidable_advection) return;
        for(int id=0;id<rigid_geometry_collection.particles.Size();id++){
            ARRAY<T,TV_INT> phi_object(fine_mac_grid.Domain_Indices(3));
            for(typename GRID<TV>::CELL_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next())
                phi_object(iterator.Cell_Index())=-rigid_geometry_collection.Rigid_Geometry(id).Implicit_Geometry_Extended_Value(iterator.Location());
            EXTRAPOLATION_UNIFORM<GRID<TV>,T> extrapolate(fine_mac_grid,phi_object,particle_levelset_evolution.particle_levelset.levelset.phi,3);extrapolate.Set_Band_Width(3);extrapolate.Extrapolate();}
    }

    void Set_Band_Width(const int band_width)
    {
        int scale=sub_scale;
        int clamp_band_width=2*max(band_width/2,scale);
        particle_levelset_evolution.particle_levelset.Set_Band_Width(clamp_band_width);        
    }

    void Extrapolate_Velocity_Across_Interface(const int band_width,const T time)
    {
        assert(band_width<8);
        int scale=sub_scale;
        int clamp_band_width=max(band_width,scale);
        int ghost=min(max(8,clamp_band_width),fine_mac_grid.Counts().Min());
        ARRAY<T,TV_INT> exchanged_phi_ghost(fine_mac_grid.Domain_Indices(ghost));
        particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(fine_mac_grid,particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,ghost);
        if(use_collidable_advection) incompressible.Extrapolate_Velocity_Across_Interface(fine_face_velocities,exchanged_phi_ghost,false,clamp_band_width,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
        else incompressible.Extrapolate_Velocity_Across_Interface(fine_face_velocities,exchanged_phi_ghost,false,clamp_band_width,0,TV());
    }

    void Sync_Fine_To_Coarse()
    {
        for(typename GRID<TV>::FACE_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()) fine_face_velocities(iterator.Full_Index())=0;
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            TV local_offset=TV::All_Ones_Vector()*.5*coarse_mac_grid.DX();
            GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>(iterator.Location()-local_offset,iterator.Location()+local_offset),true);
            ARRAY<T,FACE_INDEX<TV::dimension> > local_face_velocities(local_mac_grid);
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()) local_face_velocities(local_iterator.Full_Index())=0;
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
                if(local_iterator.First_Boundary()) local_face_velocities(local_iterator.Full_Index())=coarse_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.First_Face_Index(local_iterator.Axis())));
                else local_face_velocities(local_iterator.Full_Index())=coarse_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Second_Face_Index(local_iterator.Axis())));}
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
                fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Cell_Index()*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector()))=local_face_velocities(local_iterator.Full_Index());}}
    }


    void Sync_Coarse_To_Fine()
    {
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            TV local_offset=TV::All_Ones_Vector()*.5*coarse_mac_grid.DX();
            GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>(iterator.Location()-local_offset,iterator.Location()+local_offset),true);
            ARRAY<T,FACE_INDEX<TV::dimension> > local_face_velocities(local_mac_grid);
            Map_Coarse_To_Local_For_Cell(local_mac_grid,local_face_velocities,coarse_face_velocities,iterator,false);
            Map_Local_To_Fine_For_Cell(local_mac_grid,local_face_velocities,iterator.Cell_Index());
        }
    }

    void Set_Coarse_Phi_From_Fine_Phi(ARRAY<T,TV_INT>& coarse_phi,const ARRAY<T,TV_INT>& fine_phi)
    {
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid,1);iterator.Valid();iterator.Next())
            coarse_phi(iterator.Cell_Index())=phi_interpolation->Clamped_To_Array(fine_mac_grid,fine_phi,iterator.Location());
    }

    void Set_Local_Phi_From_Fine_Phi(GRID<TV>& local_mac_grid,ARRAY<T,TV_INT>& local_phi,const ARRAY<T,TV_INT>& fine_phi,TV_INT cell_index)
    {
        for(typename GRID<TV>::CELL_ITERATOR iterator(local_mac_grid,1);iterator.Valid();iterator.Next()){
            TV_INT fine_index=(cell_index-TV_INT::All_Ones_Vector())*sub_scale+iterator.Cell_Index();
            local_phi(iterator.Cell_Index())=fine_phi(fine_index);}
    }

    bool Map_Fine_To_Local_Boundaries_For_Cell(GRID<TV>& local_mac_grid,ARRAY<bool,FACE_INDEX<TV::dimension> >& local_psi_N,TV_INT cell_index)
    {bool has_solids=false;
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
        local_psi_N(local_iterator.Full_Index())=fine_psi_N(FACE_INDEX<TV::dimension>(local_iterator.Axis(),(cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index()));
        if(local_psi_N(local_iterator.Full_Index())) has_solids=true;}
    return has_solids;}

    void Map_Coarse_To_Fine(GRID<TV>& local_mac_grid)
    {for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
        ARRAY<T,FACE_INDEX<TV::dimension> > sum(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()));
        FACE_INDEX<TV::dimension> fine_index,sum_index=FACE_INDEX<TV::dimension>(1,TV_INT::All_Ones_Vector());
        if(alpha!=1){
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){TV_INT offset;offset(local_iterator.Axis())++;
                fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Cell_Index()*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector());
                sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(local_iterator.First_Boundary())?(TV_INT::All_Ones_Vector()):(TV_INT::All_Ones_Vector()+offset));
                if(fine_psi_N(fine_index)) sum(sum_index)+=fine_face_velocities_save(fine_index);}}
        //update based on (1-a)*((1/A)*(Vco-avg Vs))+a*(Vfo+(Vcn-Vco)*(1/A))
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
            TV_INT2 index;int offset=0;for(int i=0;i<TV::dimension;i++){if(i==local_iterator.Axis()){offset++;continue;}index(i-offset)=local_iterator.Face_Index()(i);}
            FACE_INDEX<TV::dimension> coarse_index,fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Cell_Index()*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector());
            if(local_iterator.First_Boundary()) coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_mac_grid.First_Face_Index_In_Cell(local_iterator.Axis(),iterator.Cell_Index()));
            else coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_mac_grid.Second_Face_Index_In_Cell(local_iterator.Axis(),iterator.Cell_Index()));
            T product=1;for(int i=1;i<TV::dimension;i++) product*=sub_scale;
            T area=projection.poisson->beta_face(coarse_index),one_over_area=1./area;
            T beta=fine_psi_N(fine_index)?0:1;
            if(alpha!=1){TV_INT offset;offset(local_iterator.Axis())++;sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),local_iterator.First_Boundary()?TV_INT::All_Ones_Vector():TV_INT::All_Ones_Vector()+offset);}
            if(!fine_psi_N(fine_index))
                fine_face_velocities(fine_index)=(alpha!=1?((1-alpha)*one_over_area*(coarse_face_velocities(coarse_index)-sum(sum_index)/sub_scale)):0)+
                    (alpha?(alpha*(fine_face_velocities_save(fine_index)+weights(index)*product*one_over_area*beta*(coarse_face_velocities(coarse_index)-coarse_face_velocities_save(coarse_index)))):0);}}}

    void Map_Coarse_To_Local_Interior_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,
        TV_INT coarse_cell_index,bool zero_out)
    {for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
        FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector());
        if(zero_out) local_face_velocities(local_iterator.Full_Index())=0;
        else local_face_velocities(local_iterator.Full_Index())=fine_face_velocities_save(fine_index);}}

    void Map_Coarse_To_Local_Boundary_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,
        TV_INT coarse_cell_index)
    {for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
        FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector());
        local_face_velocities(local_iterator.Full_Index())=fine_face_velocities(fine_index);}}

    void Map_Coarse_To_Local_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,
        TV_INT coarse_cell_index,bool zero_out)
    {ARRAY<T,FACE_INDEX<TV::dimension> > sum(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()));
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
        FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(coarse_cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
        if(zero_out) local_face_velocities(local_iterator.Full_Index())=0;
        else local_face_velocities(local_iterator.Full_Index())=fine_face_velocities_save(fine_index);}
    FACE_INDEX<TV::dimension> fine_index,sum_index=FACE_INDEX<TV::dimension>(1,TV_INT::All_Ones_Vector());
    if(alpha!=1){
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){TV_INT offset;offset(local_iterator.Axis())++;
            sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(local_iterator.First_Boundary())?(TV_INT::All_Ones_Vector()):(TV_INT::All_Ones_Vector()+offset));
            sum(sum_index)=0;}
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){TV_INT offset;offset(local_iterator.Axis())++;
            fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(coarse_cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
            sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(local_iterator.First_Boundary())?(TV_INT::All_Ones_Vector()):(TV_INT::All_Ones_Vector()+offset));
            if(fine_psi_N(fine_index)) sum(sum_index)+=fine_face_velocities_save(fine_index);}}
    //update based on (1-a)*((1/A)*(Vco-avg Vs))+a*(Vfo+(Vcn-Vco)*(1/A))
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
        TV_INT2 index;int offset=0;for(int i=0;i<TV::dimension;i++){if(i==local_iterator.Axis()){offset++;continue;}index(i-offset)=local_iterator.Face_Index()(i);}
        FACE_INDEX<TV::dimension> coarse_index,fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(coarse_cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
        if(local_iterator.First_Boundary()) coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_mac_grid.First_Face_Index_In_Cell(local_iterator.Axis(),coarse_cell_index));
        else coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_mac_grid.Second_Face_Index_In_Cell(local_iterator.Axis(),coarse_cell_index));
        T product=1;for(int i=1;i<TV::dimension;i++) product*=sub_scale;
        T area=beta_face(coarse_index),one_over_area=1./area;
        T beta=fine_psi_N(fine_index)?0:1;
        if(alpha!=1){TV_INT offset;offset(local_iterator.Axis())++;sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),local_iterator.First_Boundary()?TV_INT::All_Ones_Vector():TV_INT::All_Ones_Vector()+offset);}
        if(!fine_psi_N(fine_index))
            local_face_velocities(local_iterator.Full_Index())=(alpha!=1?((1-alpha)*one_over_area*(coarse_face_velocities(coarse_index)-sum(sum_index)/sub_scale)):0)+
                (alpha?(alpha*(fine_face_velocities_save(fine_index)+weights(index)*product*one_over_area*beta*(coarse_face_velocities(coarse_index)-coarse_face_velocities_save(coarse_index)))):0);}}

    T Get_Interpolated_Velocity(const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,FACE_INDEX<TV::dimension>& coarse_index,TV X,bool use_interpolation)
    {
        GRID<VECTOR<T,TV::dimension-1> >& lower_dim_grid=lower_dim_grids(coarse_index.axis);
        ARRAY<T,VECTOR<int,TV::dimension-1> >& face_velocity=lower_dim_velocities(coarse_index.axis);
        RANGE<TV_INT>& domain=lower_dim_domain;
        int range=use_cubic_interpolation?3:1;
        domain.min_corner=coarse_index.index-range*TV_INT::All_Ones_Vector();domain.max_corner=coarse_index.index+range*TV_INT::All_Ones_Vector();
        domain.min_corner(coarse_index.axis)=coarse_index.index(coarse_index.axis);domain.max_corner(coarse_index.axis)=coarse_index.index(coarse_index.axis);
        VECTOR<int,TV::dimension-1> lower_dim_coarse_index;int offset=0;for(int j=0;j<TV::dimension;j++){if(j==coarse_index.axis) offset++;else{lower_dim_coarse_index(j-offset)=coarse_index.index(j);}}
        for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid,domain,coarse_index.axis);iterator.Valid();iterator.Next()){
            VECTOR<int,TV::dimension-1> lower_dim_index;int offset=0;for(int j=0;j<TV::dimension;j++){if(j==coarse_index.axis) offset++;else{lower_dim_index(j-offset)=iterator.Face_Index()(j);}}
            TV_INT face=iterator.Face_Index();coarse_mac_grid.Clamp(face);
            if(use_interpolation){
                TV_INT face2=face+TV_INT::Axis_Vector(iterator.Axis());coarse_mac_grid.Clamp(face2);
                face_velocity(lower_dim_index)=LINEAR_INTERPOLATION<T,T>::Linear(coarse_mac_grid.Face(iterator.Axis(),face)(iterator.Axis()),coarse_mac_grid.Face(iterator.Axis(),face2)(iterator.Axis()),
                    coarse_face_velocities(FACE_INDEX<TV::dimension>(iterator.Axis(),face)),coarse_face_velocities(FACE_INDEX<TV::dimension>(iterator.Axis(),face)),X(iterator.Axis()));}
            else face_velocity(lower_dim_index)=coarse_face_velocities(FACE_INDEX<TV::dimension>(iterator.Axis(),face));}
        VECTOR<T,TV::dimension-1> pos;offset=0;for(int j=0;j<TV::dimension;j++){if(j==coarse_index.axis) offset++;else{pos(j-offset)=X(j);}}
        return interpolation->Clamped_To_Array(lower_dim_grid,face_velocity,pos);
    }

    void Map_Coarse_To_Fine_For_Cell(GRID<TV>& local_mac_grid,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,TV_INT coarse_cell_index,bool use_interpolation)
    {ARRAY<T,FACE_INDEX<TV::dimension> > sum(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()));
    FACE_INDEX<TV::dimension> fine_index,sum_index=FACE_INDEX<TV::dimension>(1,TV_INT::All_Ones_Vector());
    if(alpha!=1){
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){TV_INT offset;offset(local_iterator.Axis())++;
            sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(local_iterator.First_Boundary())?(TV_INT::All_Ones_Vector()):(TV_INT::All_Ones_Vector()+offset));
            sum(sum_index)=0;}
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){TV_INT offset;offset(local_iterator.Axis())++;
            fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(coarse_cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
            sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(local_iterator.First_Boundary())?(TV_INT::All_Ones_Vector()):(TV_INT::All_Ones_Vector()+offset));
            if(fine_psi_N(fine_index)) sum(sum_index)+=fine_face_velocities_save(fine_index);}}
    //update based on (1-a)*((1/A)*(Vco-avg Vs))+a*(Vfo+(Vcn-Vco)*(1/A))
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
        TV_INT2 index;int offset=0;for(int i=0;i<TV::dimension;i++){if(i==local_iterator.Axis()){offset++;continue;}index(i-offset)=local_iterator.Face_Index()(i);}
        FACE_INDEX<TV::dimension> coarse_index,fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(coarse_cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
        if(local_iterator.First_Boundary()) coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_mac_grid.First_Face_Index_In_Cell(local_iterator.Axis(),coarse_cell_index));
        else coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_mac_grid.Second_Face_Index_In_Cell(local_iterator.Axis(),coarse_cell_index));
        T product=1;for(int i=1;i<TV::dimension;i++) product*=sub_scale;
        T area=beta_face(coarse_index),one_over_area=1./area;
        T beta=fine_psi_N(fine_index)?0:1;
        if(alpha!=1){TV_INT offset;offset(local_iterator.Axis())++;sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),local_iterator.First_Boundary()?TV_INT::All_Ones_Vector():TV_INT::All_Ones_Vector()+offset);}
        T coarse_face_velocity,coarse_face_velocity_old;
        if(use_interpolation){
            coarse_face_velocity=Get_Interpolated_Velocity(coarse_face_velocities,coarse_index,fine_mac_grid.Face(fine_index.axis,fine_index.index),false);
            coarse_face_velocity_old=Get_Interpolated_Velocity(coarse_face_velocities_save,coarse_index,fine_mac_grid.Face(fine_index.axis,fine_index.index),false);}
        else{coarse_face_velocity=coarse_face_velocities(coarse_index);coarse_face_velocity_old=coarse_face_velocities_save(coarse_index);}
        T delta_v_coarse=coarse_face_velocity-coarse_face_velocity_old;
        bool psi_D=levelset_projection.elliptic_solver->psi_D(fine_index.index) && levelset_projection.elliptic_solver->psi_D(fine_index.index+TV_INT::Axis_Vector(fine_index.axis));
        if(!fine_psi_N(fine_index) && !psi_D) fine_face_velocities(fine_index)=(alpha!=1?((1-alpha)*one_over_area*(coarse_face_velocity-sum(sum_index)/sub_scale)):0)+
                (alpha?(alpha*(fine_face_velocities_save(fine_index)+one_over_area*beta*delta_v_coarse)):0);}
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
        FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(coarse_cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
        bool psi_D=levelset_projection.elliptic_solver->psi_D(fine_index.index) && levelset_projection.elliptic_solver->psi_D(fine_index.index+TV_INT::Axis_Vector(fine_index.axis));
        TV_INT local_left=local_iterator.Face_Index(),local_right=local_iterator.Face_Index();local_left(local_iterator.Axis())=1;local_right(local_iterator.Axis())=sub_scale+1;
        TV_INT fine_left=(coarse_cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_left,fine_right=(coarse_cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_right;
        //std::cout<<local_iterator.Axis()<<","<<coarse_cell_index<<std::endl;
        //std::cout<<"Interpoating between "<<local_left<<","<<local_right<<" fine "<<fine_left<<","<<fine_right;
        //std::cout<<" pos "<<fine_mac_grid.Face(local_iterator.Axis(),fine_left)<<","<<fine_mac_grid.Face(local_iterator.Axis(),fine_right);
        //std::cout<<" v "<<fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),fine_left))<<","<<fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),fine_right))<<std::endl;
        if(use_interpolation && !psi_D) fine_face_velocities(fine_index)=LINEAR_INTERPOLATION<T,T>::Linear(
             fine_mac_grid.Face(local_iterator.Axis(),fine_left)(local_iterator.Axis()),fine_mac_grid.Face(local_iterator.Axis(),fine_right)(local_iterator.Axis()),
             fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),fine_left)),fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),fine_right)), 
             fine_mac_grid.Face(local_iterator.Axis(),fine_index.index)(local_iterator.Axis()));
        //std::cout<<"Res is "<<fine_face_velocities(fine_index)<<" at "<<fine_mac_grid.Face(local_iterator.Axis(),fine_index.index)(local_iterator.Axis())<<std::endl;
    }}

    void Map_Local_To_Fine_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,TV_INT cell_index)
    {for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
        fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),(cell_index-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index()))=local_face_velocities(local_iterator.Full_Index());}}

    void Map_Local_To_Fine_Interior_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,TV_INT cell_index)
    {for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
        fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector()))=local_face_velocities(local_iterator.Full_Index());}}

    void Save_Velocities(const GRID<TV>& mac_grid,const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities_save)
    {for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities_save(iterator.Full_Index())=face_velocities(iterator.Full_Index());}

    void Average_Velocities_From_Fine_To_Coarse(const GRID<TV>& coarse_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& fine_face_velocities,
        const int sub_scale,const TV_INT& fine_index_offset)
    {
        for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()) coarse_face_velocities(iterator.Full_Index())=0;
        TV local_offset=TV::All_Ones_Vector()*.5*coarse_mac_grid.DX();
        GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>(-local_offset,local_offset),true); //NOTE: Grid is incorrect but does not matter
        ARRAY<T,FACE_INDEX<TV::dimension> > local_face_velocities(local_mac_grid);
        FACE_INDEX<TV::dimension> fine_index;
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
                fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),fine_index_offset+(iterator.Cell_Index()-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
                if(local_iterator.First_Boundary()) coarse_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.First_Face_Index(local_iterator.Axis())))+=fine_face_velocities(fine_index);
                else coarse_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Second_Face_Index(local_iterator.Axis())))+=fine_face_velocities(fine_index);}}    
        for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            T factor=0;
            if(iterator.First_Cell_Index()(iterator.Axis())>0){factor++;}
            if(iterator.Second_Cell_Index()(iterator.Axis())<=coarse_mac_grid.Counts()(iterator.Axis())){factor++;}
            coarse_face_velocities.Component(iterator.Axis())(iterator.Face_Index())*=factor==1?sub_scale_face_inverse:0.5*sub_scale_face_inverse;}
    }

    void Average_Velocities_From_Fine_To_Coarse(const GRID<TV>& coarse_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& fine_face_velocities,const int sub_scale)
    {
        Average_Velocities_From_Fine_To_Coarse(coarse_mac_grid,coarse_face_velocities,fine_face_velocities,sub_scale,TV_INT::Constant_Vector(0));
    }

    void Preprocess_Projection(const T dt,const T time)
    {
        //for(typename GRID<TV>::FACE_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()){coarse_face_velocities(iterator.Full_Index())=fine_face_velocities(iterator.Full_Index());}
        Average_Velocities_From_Fine_To_Coarse(coarse_mac_grid,coarse_face_velocities,fine_face_velocities,sub_scale);
        Save_Velocities(coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_save);
        Save_Velocities(fine_mac_grid,fine_face_velocities,fine_face_velocities_save);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("preprocess example",0,1);
    }
    
    void Postprocess_Projection(const T dt,const T time)
    {
        //for(typename GRID<TV>::FACE_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()){fine_face_velocities(iterator.Full_Index())=coarse_face_velocities(iterator.Full_Index());}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("postprocess example",0,1);
        if(binary_refinement_levels){
            if(use_direct_solver) Local_Projection_Analytic(dt,time,binary_refinement_levels,coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_save,beta_face,TV_INT::All_Ones_Vector());
            else Local_Projection_Analytic2(dt,time,binary_refinement_levels,coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_save,beta_face,TV_INT::All_Ones_Vector());}
        else if(number_of_threads_per_mpi_proc) Threaded_Local_Projection_PCG(dt,time);
        else Local_Projection_PCG(dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after local projection example",0,1);
    }

    void Initialize_Fast_Unwrapped_PCG_Solve(PROJECTION_UNIFORM<GRID<TV> >& local_projection)
    {
        local_projection.elliptic_solver->Set_Neumann_Outer_Boundaries();
        int number_of_elements=TV::dimension==2?sub_scale*sub_scale:sub_scale*sub_scale*sub_scale;
        cell_index_to_matrix_index.Resize(local_mac_grid.Domain_Indices());
        matrix_index_to_cell_index.Resize(number_of_elements);
        b.Resize(number_of_elements);
        ARRAY<int> row_counts;
        row_counts.Resize(number_of_elements);
        int count=0;
        for(typename GRID<TV>::CELL_ITERATOR iterator(local_mac_grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            int matrix_index;
            count++;cell_index_to_matrix_index(cell_index)=matrix_index=count;
            matrix_index_to_cell_index(matrix_index)=cell_index;}
        for(int i=0;i<row_counts.m;i++){
            int boundary=0;
            for(int j=0;j<TV::dimension;j++) if(matrix_index_to_cell_index(i)(j)==1 || matrix_index_to_cell_index(i)(j)==local_mac_grid.Counts()(j)) boundary++;
            row_counts(i)=(2*TV::dimension+1)-boundary;}
        A.Set_Row_Lengths(row_counts);
        TV one_over_dx2=Inverse(local_mac_grid.dX*local_mac_grid.dX);
        T default_row_sum=-2*one_over_dx2.L1_Norm();
        TV_INT grid_counts=local_mac_grid.counts;
        for(typename GRID<TV>::CELL_ITERATOR iterator(local_mac_grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            T row_sum=default_row_sum;
            int matrix_index=cell_index_to_matrix_index(cell_index);
            for(int axis=0;axis<GRID<TV>::dimension;axis++){TV_INT offset;offset[axis]=1;
                if(local_projection.elliptic_solver->psi_N.Component(axis)(cell_index)) row_sum+=one_over_dx2[axis];
                else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),one_over_dx2[axis]);
                if(local_projection.elliptic_solver->psi_N.Component(axis)(cell_index+offset)) row_sum+=one_over_dx2[axis];
                else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),one_over_dx2[axis]);}
            A.Set_Element(matrix_index,matrix_index,row_sum);}
    }

    void Fast_Unwrapped_PCG_Solve(PROJECTION_UNIFORM<GRID<TV> >& local_projection,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,T dt,T time)
    {
        local_projection.Compute_Divergence(typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP(face_velocities),local_projection.elliptic_solver);
        for(typename GRID<TV>::CELL_ITERATOR iterator(local_mac_grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            int matrix_index=cell_index_to_matrix_index(cell_index);
            b(matrix_index)=local_projection.elliptic_solver->f(cell_index);}
        local_projection.elliptic_solver->Solve_Subregion(matrix_index_to_cell_index,A,b);
        local_projection.Apply_Pressure(face_velocities,dt,time);
        A.Negate();
    }

    void Fast_Direct_Solve(PROJECTION_UNIFORM<GRID<TV> >& local_projection,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,T dt,T time,int scale)
    {
        for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
            TV_INT cell_index=local_iterator.Cell_Index();
            local_projection.elliptic_solver->f(cell_index)=T(0);
            for(int axis=0;axis<TV::dimension;axis++){
                TV_INT first_face_index=local_iterator.First_Face_Index(axis),second_face_index=local_iterator.Second_Face_Index(axis);
                local_projection.elliptic_solver->f(cell_index)+=face_velocities.Component(axis)(second_face_index)-face_velocities.Component(axis)(first_face_index);}}        
        PHYSBAM_FATAL_ERROR();
        //Scaled_Direct_Solvers(local_projection.p,local_projection.elliptic_solver->f,scale);
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
            TV_INT first_cell_index=local_iterator.First_Cell_Index(),second_cell_index=local_iterator.Second_Cell_Index();
            face_velocities.Component(local_iterator.Axis())(local_iterator.Face_Index())-=(local_projection.p(second_cell_index)-local_projection.p(first_cell_index));}
    }

    void Threaded_Local_Projection_PCG(const T dt,const T time)
    {
        if(!thread_queue) thread_queue=new THREAD_QUEUE(number_of_threads_per_mpi_proc);
        Map_Coarse_To_Fine(local_mac_grid);
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            REFINEMENT_TASK<TV> *refinement_task=new REFINEMENT_TASK<TV>();
            refinement_task->smoke_tests=this;refinement_task->cell_index=iterator.Cell_Index();refinement_task->dt=dt;refinement_task->time=time;
            thread_queue->Queue(refinement_task);}
        thread_queue->Wait();
    }

    void Local_Projection_PCG(const T dt,const T time)
    {
        //bool ring_level=0;
        /*Set_Levelset_Boundary_Conditions(fine_mac_grid,fine_face_velocities,particle_levelset_evolution.phi);
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            if(!All_Cells_Inside(iterator.Cell_Index(),ring_level,particle_levelset_evolution.phi)){
                Map_Coarse_To_Fine_For_Cell(local_mac_grid,coarse_face_velocities,iterator.Cell_Index(),true);
            }}*/
        PHYSBAM_DEBUG_WRITE_SUBSTEP("done with mapping",0,1);
        LOG::Time("Fine Projection Small Cells");
        int ghost=max(buffer+sub_scale,3);
        ARRAY<T,TV_INT> phi_ghost(fine_mac_grid.Domain_Indices(ghost));
        phi_boundary->Fill_Ghost_Cells(fine_mac_grid,particle_levelset_evolution.phi,phi_ghost,dt,time,ghost);
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            //LOG::Time("Fine Projection Small Cell");
            if(surface_solve && Contains_Outside(iterator.Cell_Index(),phi_ghost,buffer)) continue;
            if(!surface_solve && !Contains_Inside(iterator.Cell_Index(),phi_ghost,buffer)) continue;
            Map_Coarse_To_Local_For_Cell(local_mac_grid,local_face_velocities,coarse_face_velocities,iterator.Cell_Index(),false);
            for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
                //local_projection.p(local_iterator.Cell_Index())=projection.p(iterator.Cell_Index());}
                local_projection.p(local_iterator.Cell_Index())=0;}
            bool contains_solids=Map_Fine_To_Local_Boundaries_For_Cell(local_mac_grid,local_projection.elliptic_solver->psi_N,iterator.Cell_Index());
            bool contains_surface=false;
            if(!surface_solve) contains_surface=Set_Local_Boundary_Conditions(local_mac_grid,local_projection,iterator.Cell_Index());
            local_projection.elliptic_solver->Set_Neumann_Outer_Boundaries();
            local_projection.p*=dt;
            if(!surface_solve){
                Set_Local_Phi_From_Fine_Phi(local_mac_grid,local_phi,particle_levelset_evolution.phi,iterator.Cell_Index());
                local_projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();}
            if(contains_solids || contains_surface) local_projection.Make_Divergence_Free(local_face_velocities,dt,time);
            else if(use_direct_solver) Fast_Direct_Solve(local_projection,local_face_velocities,dt,time,sub_scale);
            else Fast_Unwrapped_PCG_Solve(local_projection,local_face_velocities,dt,time);
            local_projection.p/=dt;
            Map_Local_To_Fine_For_Cell(local_mac_grid,local_face_velocities,iterator.Cell_Index());
            if(!surface_solve) local_projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);}
        if(surface_solve){
            LOG::Time("Fine Projection Surface Boundary");
            for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
                for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
                    TV_INT fine_index=(iterator.Cell_Index()-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Cell_Index();
                    levelset_projection.p(fine_index)=0;/*projection.p(iterator.Cell_Index())*/;}}
            LOG::Time("Fine Projection Boundary Conditions");
            Set_Levelset_Boundary_Conditions(fine_mac_grid,fine_face_velocities,phi_ghost,time);
            LOG::Time("Fine Projection Exchange");
            if(fine_mpi_grid){
                fine_mpi_grid->Union_Common_Face_Data(levelset_projection.elliptic_solver->psi_N);
                //fine_mpi_grid->Exchange_Boundary_Face_Data(levelset_projection.elliptic_solver->psi_N,1);
                fine_mpi_grid->Exchange_Boundary_Cell_Data(levelset_projection.elliptic_solver->psi_D,1,false);
                fine_mpi_grid->Exchange_Boundary_Cell_Data(levelset_projection.p,1,false);}
            //levelset_incompressible.Set_Dirichlet_Boundary_Conditions(&particle_levelset_evolution.phi,0);        
            levelset_projection.p*=dt;
            //Set_Levelset_Phi_From_Fine_Phi(local_mac_grid,local_phi,particle_levelset_evolution.phi,iterator.Cell_Index());
            levelset_projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
            PHYSBAM_DEBUG_WRITE_SUBSTEP("before surface solve",0,1);
            LOG::Time("Fine Projection");
            boundary->Apply_Boundary_Condition_Face(levelset_projection.p_grid,fine_face_velocities,time+dt);        
            levelset_projection.Make_Divergence_Free(fine_face_velocities,dt,time);
            //levelset_incompressible.Advance_One_Time_Step_Implicit_Part(fine_face_velocities,dt,time,false,0);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after surface solve",0,1);
            levelset_projection.p/=dt;
            levelset_projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);}
    }

    void Initialize_Bodies()
    {}

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
    {
        twist.linear=TV();
        if(time>1) twist.linear(2)=1;
        return true;
    }

    //Analytic
    
    void Initialize_Binary_Local_Solve()
    {
        local_face_velocities_array.Resize(binary_refinement_levels);
        local_face_velocities_save_array.Resize(binary_refinement_levels);
        local_projection_array.Resize(binary_refinement_levels);
        for(int i=0;i<binary_refinement_levels;i++){
            local_face_velocities_array(i)=new ARRAY<T,FACE_INDEX<TV::dimension> >(local_mac_grid);
            local_face_velocities_save_array(i)=new ARRAY<T,FACE_INDEX<TV::dimension> >(local_mac_grid);
            local_projection_array(i)=new PROJECTION_UNIFORM<GRID<TV> >(local_mac_grid,true,true);
            local_projection_array(i)->poisson->Set_Variable_beta(true);
            local_projection_array(i)->elliptic_solver->pcg.Set_Maximum_Iterations(40);
            local_projection_array(i)->elliptic_solver->pcg.Use_Modified_Incomplete_Cholesky();
        }
    }
    
    inline int Cell_Boundary_Type(const ARRAY<T,FACE_INDEX<TV::dimension> >& beta_face,const typename GRID<TV>::CELL_ITERATOR& iterator)
    {
        bool has_solid=false,has_fluid=false;
        for(int axis=0;axis<TV::dimension;axis++){
            if(beta_face.Component(axis)(iterator.First_Face_Index(axis))<1) has_solid=true;
            if(beta_face.Component(axis)(iterator.Second_Face_Index(axis))<1) has_solid=true;
            if(beta_face.Component(axis)(iterator.First_Face_Index(axis))>0) has_fluid=true;
            if(beta_face.Component(axis)(iterator.Second_Face_Index(axis))>0) has_fluid=true;}
        if(has_solid && has_fluid) return 2; // mixed with solids and fluids
        else if(has_fluid) return 1; // fluids only
        else return 0; // solids only
    }

    void Average_Velocities_From_Fine_To_Coarse_For_Analytic(ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& fine_face_velocities,const GRID<TV>& coarse_mac_grid,
        const TV_INT& fine_index_offset,const int sub_scale)
    {
        if(sub_scale==1) for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()) 
            coarse_face_velocities.Component(iterator.Axis())(iterator.Face_Index())=fine_face_velocities.Component(iterator.Axis())(iterator.Face_Index()+fine_index_offset);
        else Average_Velocities_From_Fine_To_Coarse(coarse_mac_grid,coarse_face_velocities,fine_face_velocities,sub_scale,fine_index_offset);
    }

    void Map_Coarse_To_Local_For_Cell_For_Analytic(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,
        const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities_save,typename GRID<TV>::CELL_ITERATOR& coarse_iterator,
        const ARRAY<bool,FACE_INDEX<TV::dimension> >& local_psi_N,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_beta_face,const ARRAY<T,FACE_INDEX<TV::dimension> >& local_beta_face,bool zero_out)
    {
        ARRAY<T,FACE_INDEX<TV::dimension> > sum(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()));
        if(zero_out) for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next())
            local_face_velocities.Component(local_iterator.Axis())(local_iterator.Face_Index())=0;
        FACE_INDEX<TV::dimension> fine_index,sum_index=FACE_INDEX<TV::dimension>(1,TV_INT::All_Ones_Vector());
        if(alpha!=1){
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){TV_INT offset;offset(local_iterator.Axis())++;
                sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),local_iterator.First_Boundary()?TV_INT::All_Ones_Vector():(TV_INT::All_Ones_Vector()+offset));
                if(local_psi_N(local_iterator.Full_Index())) sum(sum_index)+=local_face_velocities(local_iterator.Full_Index());}}
        FACE_INDEX<TV::dimension> coarse_index;
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
            TV_INT2 index;int offset=0;for(int i=0;i<TV::dimension;i++){if(i==local_iterator.Axis()){offset++;continue;}index(i-offset)=local_iterator.Face_Index()(i);}
            if(local_iterator.First_Boundary()) coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_iterator.First_Face_Index(local_iterator.Axis()));
            else coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_iterator.Second_Face_Index(local_iterator.Axis()));
            T area=coarse_beta_face(coarse_index),one_over_area=1./area;
            if(alpha!=1){TV_INT offset;offset(local_iterator.Axis())++;sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),local_iterator.First_Boundary()?TV_INT::All_Ones_Vector():TV_INT::All_Ones_Vector()+offset);}
            if(!local_psi_N.Component(local_iterator.Axis())(local_iterator.Face_Index()) && area!=0) 
                local_face_velocities(local_iterator.Full_Index())=(alpha!=1?((1-alpha)*one_over_area*(coarse_face_velocities(coarse_index)-sum(sum_index)/sub_scale)):0)+
                    (alpha?
                    (alpha*(local_face_velocities(local_iterator.Full_Index())+one_over_area*local_beta_face(local_iterator.Full_Index())*(coarse_face_velocities(coarse_index)-coarse_face_velocities_save(coarse_index)))):0);}
    }

    int Set_Local_Boundary_Conditions_For_Analytic(const GRID<TV>& coarse_mac_grid,const TV_INT fine_index_offset,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,const ARRAY<bool,FACE_INDEX<TV::dimension> >& fine_psi_N,
        ARRAY<T,FACE_INDEX<TV::dimension> >& beta_face,const int sub_scale)
    {
        bool contains_boundaries=false,contains_fluids=false;
        if(sub_scale==1){
            for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
                FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(iterator.Axis(),fine_index_offset+iterator.Face_Index());
                if(fine_psi_N(fine_index)){
                    beta_face.Component(iterator.Axis())(iterator.Face_Index())=0;
                    psi_N.Component(iterator.Axis())(iterator.Face_Index())=true;
                    contains_boundaries=true;}
                else{
                    beta_face.Component(iterator.Axis())(iterator.Face_Index())=1;
                    psi_N.Component(iterator.Axis())(iterator.Face_Index())=false;
                    contains_fluids=true;}}}
        else{
            for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
                beta_face.Component(iterator.Axis())(iterator.Face_Index())=0;
                psi_N.Component(iterator.Axis())(iterator.Face_Index())=false;}
            TV local_offset=TV::All_Ones_Vector()*.5*coarse_mac_grid.DX();
            GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>(-local_offset,local_offset),true);
            FACE_INDEX<TV::dimension> fine_index;
            for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
                for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
                    fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),fine_index_offset+(iterator.Cell_Index()-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
                    if(!fine_psi_N(fine_index)){
                        contains_fluids=true;
                        if(local_iterator.First_Boundary()) beta_face(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.First_Face_Index(local_iterator.Axis())))+=T(1);
                        else beta_face(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Second_Face_Index(local_iterator.Axis())))+=T(1);}
                    else contains_boundaries=true;}}
            for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
                if(beta_face.Component(iterator.Axis())(iterator.Face_Index())==0){
                    psi_N.Component(iterator.Axis())(iterator.Face_Index())=true;}
                else{
                    int factor=0;
                    if(iterator.First_Cell_Index()(iterator.Axis())>0){factor++;}
                    if(iterator.Second_Cell_Index()(iterator.Axis())<=coarse_mac_grid.Counts()(iterator.Axis())){factor++;}
                    beta_face.Component(iterator.Axis())(iterator.Face_Index())*=(factor==1)?sub_scale_face_inverse:0.5*sub_scale_face_inverse;}}}
        if(contains_boundaries && contains_fluids) return 1; // mixed with solids and fluids
        else if(!contains_boundaries && contains_fluids) return 0; // only fluids
        else return -1; // only solids
    }

    void Solve_For_Pressure_Analytically(ARRAY<T,VECTOR<int,2> >& p,ARRAY<T,VECTOR<int,2> >& div)
    {
        typedef VECTOR<int,2> TV_INT;
        p(TV_INT(2,2))=0;
        p(TV_INT(2,1))=T(-0.5)*div(TV_INT(1,1))+T(-0.25)*div(TV_INT(1,2))+T(-0.75)*div(TV_INT(2,1));
        p(TV_INT(1,2))=T(-0.5)*div(TV_INT(1,1))+T(-0.75)*div(TV_INT(1,2))+T(-0.25)*div(TV_INT(2,1));
        p(TV_INT(1,1))=T(-1)  *div(TV_INT(1,1))+T(-0.5) *div(TV_INT(1,2))+T(-0.5) *div(TV_INT(2,1));
    }

    void Solve_For_Pressure_Analytically(ARRAY<T,VECTOR<int,3> >& p,ARRAY<T,VECTOR<int,3> >& div)
    {
        typedef VECTOR<int,3> TV_INT;
        p(TV_INT(1,1,1))=T(-5./6.)*div(TV_INT(1,1,1))+T(-1./2.)*div(TV_INT(1,1,2))+T(-1./2.)*div(TV_INT(1,2,1))+T(-1./3.) *div(TV_INT(1,2,2))+T(-1./2.)*div(TV_INT(2,1,1))+T(-1./3.) *div(TV_INT(2,1,2))+T(-1./3.) *div(TV_INT(2,2,1));
        p(TV_INT(1,1,2))=T(-1./2.)*div(TV_INT(1,1,1))+T(-3./4.)*div(TV_INT(1,1,2))+T(-3./8.)*div(TV_INT(1,2,1))+T(-3./8.) *div(TV_INT(1,2,2))+T(-3./8.)*div(TV_INT(2,1,1))+T(-3./8.) *div(TV_INT(2,1,2))+T(-1./4.) *div(TV_INT(2,2,1));
        p(TV_INT(1,2,1))=T(-1./2.)*div(TV_INT(1,1,1))+T(-3./8.)*div(TV_INT(1,1,2))+T(-3./4.)*div(TV_INT(1,2,1))+T(-3./8.) *div(TV_INT(1,2,2))+T(-3./8.)*div(TV_INT(2,1,1))+T(-1./4.) *div(TV_INT(2,1,2))+T(-3./8.) *div(TV_INT(2,2,1));
        p(TV_INT(1,2,2))=T(-1./3.)*div(TV_INT(1,1,1))+T(-3./8.)*div(TV_INT(1,1,2))+T(-3./8.)*div(TV_INT(1,2,1))+T(-7./12.)*div(TV_INT(1,2,2))+T(-1./4.)*div(TV_INT(2,1,1))+T(-5./24.)*div(TV_INT(2,1,2))+T(-5./24.)*div(TV_INT(2,2,1));
        p(TV_INT(2,1,1))=T(-1./2.)*div(TV_INT(1,1,1))+T(-3./8.)*div(TV_INT(1,1,2))+T(-3./8.)*div(TV_INT(1,2,1))+T(-1./4.) *div(TV_INT(1,2,2))+T(-3./4.)*div(TV_INT(2,1,1))+T(-3./8.) *div(TV_INT(2,1,2))+T(-3./8.) *div(TV_INT(2,2,1));
        p(TV_INT(2,1,2))=T(-1./3.)*div(TV_INT(1,1,1))+T(-3./8.)*div(TV_INT(1,1,2))+T(-1./4.)*div(TV_INT(1,2,1))+T(-5./24.)*div(TV_INT(1,2,2))+T(-3./8.)*div(TV_INT(2,1,1))+T(-7./12.)*div(TV_INT(2,1,2))+T(-5./24.)*div(TV_INT(2,2,1));
        p(TV_INT(2,2,1))=T(-1./3.)*div(TV_INT(1,1,1))+T(-1./4.)*div(TV_INT(1,1,2))+T(-3./8.)*div(TV_INT(1,2,1))+T(-5./24.)*div(TV_INT(1,2,2))+T(-3./8.)*div(TV_INT(2,1,1))+T(-5./24.)*div(TV_INT(2,1,2))+T(-7./12.)*div(TV_INT(2,2,1));
        p(TV_INT(2,2,2))=T(0);
    }

    void Local_Projection_Analytic(const T dt,const T time,const int level,const GRID<TV>& coarse_mac_grid,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,
                                   const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities_save,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_beta_face,const TV_INT& parent_cell_index)
    {
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            Average_Velocities_From_Fine_To_Coarse_For_Analytic(*local_face_velocities_array(level),fine_face_velocities_save,local_mac_grid,
                ((parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index()-TV_INT::All_Ones_Vector())*(1<<level),1<<(level-1));
            Save_Velocities(local_mac_grid,*local_face_velocities_array(level),*local_face_velocities_save_array(level));
            int boundary_code=Set_Local_Boundary_Conditions_For_Analytic(local_mac_grid,((parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index()-TV_INT::All_Ones_Vector())*(1<<level),
                local_projection_array(level)->poisson->psi_N,fine_psi_N,local_projection_array(level)->poisson->beta_face,1<<(level-1));
            if(boundary_code==1){
                for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
                    local_projection_array(level)->p(local_iterator.Cell_Index())=projection.p(iterator.Cell_Index());}
                Map_Coarse_To_Local_For_Cell_For_Analytic(local_mac_grid,*local_face_velocities_array(level),coarse_face_velocities,coarse_face_velocities_save,
                    iterator,local_projection_array(level)->poisson->psi_N,coarse_beta_face,local_projection_array(level)->poisson->beta_face,false);
                local_projection_array(level)->poisson->Set_Neumann_Outer_Boundaries();
                local_projection_array(level)->poisson->beta_given_on_faces=true;
                local_projection_array(level)->Make_Divergence_Free(*local_face_velocities_array(level),dt,time);}
            else if(boundary_code==0){
                Map_Coarse_To_Local_For_Cell_For_Analytic(local_mac_grid,*local_face_velocities_array(level),coarse_face_velocities,coarse_face_velocities_save,
                    iterator,local_projection_array(level)->poisson->psi_N,coarse_beta_face,local_projection_array(level)->poisson->beta_face,false);
                for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
                    TV_INT cell_index=local_iterator.Cell_Index();
                    local_projection_array(level)->poisson->f(cell_index)=T(0);
                    for(int axis=0;axis<TV::dimension;axis++){
                        TV_INT first_face_index=local_iterator.First_Face_Index(axis),second_face_index=local_iterator.Second_Face_Index(axis);
                        local_projection_array(level)->poisson->f(cell_index)+=local_face_velocities_array(level)->Component(axis)(second_face_index)-local_face_velocities_array(level)->Component(axis)(first_face_index);}}
                Solve_For_Pressure_Analytically(local_projection_array(level)->p,local_projection_array(level)->poisson->f);
                for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
                    TV_INT first_cell_index=local_iterator.First_Cell_Index(),second_cell_index=local_iterator.Second_Cell_Index();
                    (*local_face_velocities_array(level))(local_iterator.Full_Index())-=
                        local_projection_array(level)->poisson->beta_face(local_iterator.Full_Index())*(local_projection_array(level)->p(second_cell_index)-local_projection_array(level)->p(first_cell_index));}}
            if(level>1)
                Local_Projection_Analytic(dt,time,level-1,local_mac_grid,*local_face_velocities_array(level),*local_face_velocities_save_array(level),
                    local_projection_array(level)->poisson->beta_face,(parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index());
            else for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()) 
                fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),((parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index()-TV_INT::All_Ones_Vector())*2+local_iterator.Face_Index()))=
                    (*local_face_velocities_array(level)).Component(local_iterator.Axis())(local_iterator.Face_Index());}
    }

    void Local_Projection_Analytic2(const T dt,const T time,const int level,const GRID<TV>& coarse_mac_grid,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,
        const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities_save,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_beta_face,const TV_INT& parent_cell_index)
    {
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            if(level==binary_refinement_levels && Cell_Boundary_Type(coarse_beta_face,iterator)==2){
                Average_Velocities_From_Fine_To_Coarse_For_Analytic(local_face_velocities,fine_face_velocities_save,local_mac_grid,
                    ((parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index()-TV_INT::All_Ones_Vector())*sub_scale,1);
                Set_Local_Boundary_Conditions_For_Analytic(local_mac_grid,((parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index()-TV_INT::All_Ones_Vector())*sub_scale,
                    local_projection_analytic.poisson->psi_N,fine_psi_N,local_projection_analytic.poisson->beta_face,1);
                for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
                    local_projection_analytic.p(local_iterator.Cell_Index())=projection.p(iterator.Cell_Index());}
                Map_Coarse_To_Local_For_Cell_For_Analytic(local_mac_grid,local_face_velocities,coarse_face_velocities,coarse_face_velocities_save,iterator,
                    local_projection_analytic.poisson->psi_N,coarse_beta_face,local_projection_analytic.poisson->beta_face,false);
                local_projection_analytic.poisson->Set_Neumann_Outer_Boundaries();
                local_projection_analytic.poisson->beta_given_on_faces=true;
                local_projection_analytic.Make_Divergence_Free(local_face_velocities,dt,time);
                for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next())
                    fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),((parent_cell_index-TV_INT::All_Ones_Vector())*sub_scale+iterator.Cell_Index()-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index()))=
                        (local_face_velocities).Component(local_iterator.Axis())(local_iterator.Face_Index());}
            else if(Cell_Boundary_Type(coarse_beta_face,iterator)==1){
                Average_Velocities_From_Fine_To_Coarse_For_Analytic(*local_face_velocities_array(level),fine_face_velocities_save,local_mac_grid,
                    ((parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index()-TV_INT::All_Ones_Vector())*(1<<level),1<<(level-1));
                Save_Velocities(local_mac_grid,*local_face_velocities_array(level),*local_face_velocities_save_array(level));
                Set_Local_Boundary_Conditions_For_Analytic(local_mac_grid,((parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index()-TV_INT::All_Ones_Vector())*(1<<level),
                    local_projection_array(level)->poisson->psi_N,fine_psi_N,local_projection_array(level)->poisson->beta_face,1<<(level-1));
                Map_Coarse_To_Local_For_Cell_For_Analytic(local_mac_grid,*local_face_velocities_array(level),coarse_face_velocities,coarse_face_velocities_save,iterator,
                    local_projection_array(level)->poisson->psi_N,coarse_beta_face,local_projection_array(level)->poisson->beta_face,false);
                for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
                    TV_INT cell_index=local_iterator.Cell_Index();
                    local_projection_array(level)->poisson->f(cell_index)=T(0);
                    for(int axis=0;axis<TV::dimension;axis++){
                        TV_INT first_face_index=local_iterator.First_Face_Index(axis),second_face_index=local_iterator.Second_Face_Index(axis);
                        local_projection_array(level)->poisson->f(cell_index)+=local_face_velocities_array(level)->Component(axis)(second_face_index)-local_face_velocities_array(level)->Component(axis)(first_face_index);}}
                Solve_For_Pressure_Analytically(local_projection_array(level)->p,local_projection_array(level)->poisson->f);
                for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
                    TV_INT first_cell_index=local_iterator.First_Cell_Index(),second_cell_index=local_iterator.Second_Cell_Index();
                    (*local_face_velocities_array(level)).Component(local_iterator.Axis())(local_iterator.Face_Index())-=
                        local_projection_array(level)->poisson->beta_face.Component(local_iterator.Axis())(local_iterator.Face_Index())*(local_projection_array(level)->p(second_cell_index)-local_projection_array(level)->p(first_cell_index));}
                if(level>1)
                    Local_Projection_Analytic(dt,time,level-1,local_mac_grid,*local_face_velocities_array(level),*local_face_velocities_save_array(level),
                        local_projection_array(level)->poisson->beta_face,(parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index());
                else for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next())
                    fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),((parent_cell_index-TV_INT::All_Ones_Vector())*2+iterator.Cell_Index()-TV_INT::All_Ones_Vector())*2+local_iterator.Face_Index()))=
                        (*local_face_velocities_array(level)).Component(local_iterator.Axis())(local_iterator.Face_Index());}}
    }

//#####################################################################
};
}
#endif
