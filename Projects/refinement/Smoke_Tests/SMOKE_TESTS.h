//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// 1. a plume of smoke
// 2. smoke with a static rigid sphere
// 3. smoke with a kinematic rigid sphere
//#####################################################################
#ifndef __SMOKE_TESTS__
#define __SMOKE_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_REFINEMENT_EXAMPLE.h>
#include "DIRECT_SOLVERS.h"
#include "REFINEMENT_THREADS.h"

namespace PhysBAM{

template<class TV>
class SMOKE_TESTS:public INCOMPRESSIBLE_REFINEMENT_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef VECTOR<int,TV::dimension-1> TV_INT2;
    typedef INCOMPRESSIBLE_REFINEMENT_EXAMPLE<TV> BASE;

public:
    using BASE::coarse_mac_grid;using BASE::fine_mac_grid;using BASE::coarse_face_velocities;using BASE::fine_face_velocities;using BASE::coarse_mpi_grid;using BASE::fine_mpi_grid;
    using BASE::incompressible;using BASE::projection; using BASE::output_directory;using BASE::write_substeps_level;using BASE::last_frame;using BASE::domain_boundary;using BASE::cfl;using BASE::split_dir;
    using BASE::sub_scale;using BASE::use_coarse_forces;using BASE::use_interpolated_vorticity;using BASE::kolmogorov;using BASE::domain_boundary;using BASE::restart;using BASE::rigid_geometry_collection;using BASE::thread_queue;using BASE::boundary_scalar;using BASE::density;

    int scale;
    T sub_scale_face_inverse; //tmp data to make more efficient operations    
    bool use_collisions;
    int binary_refinement_levels;
    ARRAY<T,TV_INT2> weights;
    T alpha;
    T vc;
    ARRAY<T,FACE_INDEX<TV::dimension> > fine_face_velocities_save,coarse_face_velocities_save;
    CYLINDER<T> source_box_3d;
    RANGE<TV> source_box_2d;
    ARRAY<bool,FACE_INDEX<TV::dimension> > fine_psi_N;
    bool use_direct_solver;
    int number_of_threads;
    T buoyancy_clamp;

    GRID<TV> local_mac_grid;
    ARRAY<T,FACE_INDEX<TV::dimension> > local_face_velocities;
    PROJECTION_UNIFORM<GRID<TV> > local_projection,local_projection_analytic;

    ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> >*> local_face_velocities_array,local_face_velocities_save_array;
    ARRAY<PROJECTION_UNIFORM<GRID<TV> >*> local_projection_array;

    SPARSE_MATRIX_FLAT_NXN<T> A;
    VECTOR_ND<T> b;
    ARRAY<int,TV_INT> cell_index_to_matrix_index;
    ARRAY<TV_INT,int> matrix_index_to_cell_index;
    int test_number;
    STREAM_TYPE stream_type;
    int rigid_particle_id;
    int short_kol;

    SMOKE_TESTS(const STREAM_TYPE stream_type_input,const PARSE_ARGS& parse_args)
        :INCOMPRESSIBLE_REFINEMENT_EXAMPLE<TV>(stream_type_input),use_collisions(false),binary_refinement_levels(0),
        use_direct_solver(false),number_of_threads(0),
        local_mac_grid(TV_INT::All_Ones_Vector(),RANGE<TV>::Centered_Box(),true),local_face_velocities(local_mac_grid),local_projection(local_mac_grid),local_projection_analytic(local_mac_grid),stream_type(stream_type_input)
    {
        test_number=parse_args.Get_Integer_Value("-test_number");
        last_frame=parse_args.Get_Integer_Value("-last_frame");
        restart=parse_args.Get_Integer_Value("-restart");
        alpha=parse_args.Get_Double_Value("-alpha");
        write_substeps_level=parse_args.Get_Integer_Value("-substep");
        number_of_threads=parse_args.Get_Integer_Value("-threads");
        scale=parse_args.Get_Integer_Value("-scale");
        binary_refinement_levels=parse_args.Get_Integer_Value("-binary");
        split_dir=parse_args.Get_String_Value("-split");
        buoyancy_clamp=parse_args.Get_Double_Value("-buoyancy");
        if(parse_args.Is_Value_Set("-coarse_forces")){use_coarse_forces=true;use_interpolated_vorticity=true;}
        else use_coarse_forces=false;
        if(binary_refinement_levels){
            sub_scale=1;
            for(int i=0;i<binary_refinement_levels;i++) sub_scale<<=1;}
        else sub_scale=parse_args.Get_Integer_Value("-subscale");
        if(sub_scale==1) LOG::cout<<"WARNING: Please use one_phase_incompressible instead"<<std::endl;
        vc=parse_args.Get_Double_Value("-vc");
        cfl=parse_args.Get_Double_Value("-cfl");
        kolmogorov=parse_args.Get_Double_Value("-kol");
        short_kol=parse_args.Get_Integer_Value("-short_kol");
        
        TV_INT ratio=TV_INT::All_Ones_Vector();
        if(parse_args.Get_Integer_Value("-x")!=1) ratio(0)*=parse_args.Get_Integer_Value("-x");
        if(parse_args.Get_Integer_Value("-y")!=1) ratio(1)*=parse_args.Get_Integer_Value("-y");
        if(parse_args.Get_Integer_Value("-z")!=1) ratio(2)*=parse_args.Get_Integer_Value("-z");
        TV space;space(0)=T(1);for(int i=2;i<=TV::dimension;i++) space(i)=(T)ratio(i)/(T)ratio(0);
        space*=0.5;

        T source_radius=parse_args.Get_Double_Value("-source_radius");
        if(TV::dimension==2){
            source_box_2d.min_corner=space*0.5-TV::Constant_Vector(source_radius);
            source_box_2d.max_corner=space*0.5+TV::Constant_Vector(source_radius);
            source_box_2d.min_corner.y=T(0);
            source_box_2d.max_corner.y=T(0.05);}
        else if(TV::dimension==3){
            source_box_3d.radius=source_radius;
            VECTOR<T,3> endpoint1=VECTOR<T,3>::Constant_Vector(space(0)*0.5),endpoint2=VECTOR<T,3>::Constant_Vector(space(0)*0.5);
            endpoint1(1)=T(0);endpoint2(1)=T(0.05);
            source_box_3d.Set_Endpoints(endpoint1,endpoint2);}
        else{PHYSBAM_NOT_IMPLEMENTED();}

        if(binary_refinement_levels) Initialize_Binary_Local_Solve();
        fine_mac_grid.Initialize(ratio*scale,RANGE<TV>(TV(),space),true);
        coarse_mac_grid.Initialize(ratio*scale/sub_scale,RANGE<TV>(TV(),space),true);
        local_mac_grid.Initialize(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>::Centered_Box(),true);
        local_projection.Initialize_Grid(local_mac_grid);
        local_projection_analytic.Initialize_Grid(local_mac_grid);
        local_projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
        local_projection.elliptic_solver->pcg.Use_Modified_Incomplete_Cholesky();
        if(!use_direct_solver) Initialize_Fast_Unwrapped_PCG_Solve(local_projection);
        local_face_velocities.Resize(local_mac_grid);
        weights.Resize(RANGE<TV_INT2>(TV_INT2::All_Ones_Vector(),TV_INT2::All_Ones_Vector()*sub_scale));
        Initialize_Weights(weights,sub_scale);
        if(number_of_threads) thread_queue=new THREAD_QUEUE(number_of_threads);    
        PHYSBAM_ASSERT(scale%sub_scale==0 && binary_refinement_levels>=0);
        sub_scale_face_inverse=(TV::dimension==2)?(1./sub_scale):(1./(sub_scale*sub_scale));
        incompressible.Set_Body_Force(true);
        output_directory=STRING_UTILITIES::string_sprintf("Smoke_Tests/Test_%d_%d_%d%s_%1.2f_%s_vc_%1.2f_short_%d_kol_%1.3f_%dd",test_number,scale,sub_scale,binary_refinement_levels?"_binary":"",alpha,use_coarse_forces?"coarse":"fine",vc,short_kol,kolmogorov,TV::dimension);

     }

    ~SMOKE_TESTS()
    {
        for(int i=0;i<binary_refinement_levels;i++){
            delete local_face_velocities_array(i);
            delete local_face_velocities_save_array(i);
            delete local_projection_array(i);}
    }

    void Initialize_MPI()
    {
        fine_psi_N.Resize(fine_mac_grid);
        coarse_face_velocities_save.Resize(coarse_mac_grid);
        fine_face_velocities_save.Resize(fine_mac_grid);
        incompressible.force.Resize(fine_mac_grid,1);
    }

    void Initialize_Variable_Confinement(GRID<TV>& grid)
    {
        if(vc==0) return;
        incompressible.Use_Variable_Vorticity_Confinement(grid,true);
        for(typename GRID<TV>::CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();T scale=location.y*0.5;
            if(scale) incompressible.variable_vorticity_confinement(iterator.Cell_Index())=scale*vc;}
     }

    void Initialize_Confinement()
    {
        incompressible.Set_Vorticity_Confinement(vc);
//        if(use_coarse_forces) Initialize_Variable_Confinement(coarse_mac_grid);
//        else Initialize_Variable_Confinement(fine_mac_grid);
    }

    void Initialize_Weights(ARRAY<T,VECTOR<int,0> >& weights,const int sub_scale)
    {PHYSBAM_NOT_IMPLEMENTED();}
    
    void Initialize_Weights(ARRAY<T,VECTOR<int,1> >& weights,const int sub_scale)
    {for(int i=0;i<sub_scale;i++) weights(TV_INT2(i))=1/(T)sub_scale;}
    
    void Initialize_Weights(ARRAY<T,VECTOR<int,2> >& weights,const int sub_scale)
    {for(int i=0;i<sub_scale;i++) for(int j=0;j<sub_scale;j++) weights(TV_INT2(i,j))=1/(T)(sub_scale*sub_scale);}

    void Write_Output_Files(const int frame)
    {BASE::Write_Output_Files(frame);}

    void Set_Coarse_Boundary_Conditions()
    {for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*(axis-1)+axis_side;
        if(domain_boundary(axis)(axis_side)){ //Need to check mpi as smaller solves are never using mpi (for now)
            TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);    
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(coarse_mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);local_iterator.Valid();local_iterator.Next()){TV_INT cell=local_iterator.Face_Index()+interior_cell_offset;
                TV_INT boundary_face=axis_side==1?local_iterator.Face_Index()+TV_INT::Axis_Vector(axis):local_iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if(axis!=2){ if(coarse_face_velocities.Component(axis).Valid_Index(boundary_face)){projection.elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(axis,boundary_face))=true;coarse_face_velocities(FACE_INDEX<TV::dimension>(axis,boundary_face))=0;}}
                else{projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()) projection.poisson->beta_face.Component(iterator.Axis())(iterator.Face_Index())=0;
    TV local_offset=TV::All_Ones_Vector()*.5*coarse_mac_grid.DX();
    GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>(-local_offset,local_offset),true);
    FACE_INDEX<TV::dimension> fine_index;
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
            fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(iterator.Cell_Index()-TV_INT::All_Ones_Vector())*sub_scale+local_iterator.Face_Index());
            if(!fine_psi_N(fine_index)){
                if(local_iterator.First_Boundary()) projection.poisson->beta_face(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.First_Face_Index(local_iterator.Axis())))+=T(1);
                else projection.poisson->beta_face(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Second_Face_Index(local_iterator.Axis())))+=T(1);}}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
        if(projection.poisson->beta_face.Component(iterator.Axis())(iterator.Face_Index())==0) projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
        else{
            int factor=0;
            if(iterator.First_Cell_Index()(iterator.Axis())>0){factor++;}
            if(iterator.Second_Cell_Index()(iterator.Axis())<=coarse_mac_grid.Counts()(iterator.Axis())){factor++;}
            projection.poisson->beta_face.Component(iterator.Axis())(iterator.Face_Index())*=(factor==1)?sub_scale_face_inverse:0.5*sub_scale_face_inverse;}}}

    void Set_Fine_Boundary_Conditions(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N)
    {for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++) if(domain_boundary(axis)(axis_side)){int side=2*(axis-1)+axis_side;
         for(typename GRID<TV>::FACE_ITERATOR local_iterator(fine_mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);local_iterator.Valid();local_iterator.Next()){TV_INT boundary_face=axis_side==1?local_iterator.Face_Index()+TV_INT::Axis_Vector(axis):local_iterator.Face_Index()-TV_INT::Axis_Vector(axis);
             if(axis!=2 && fine_face_velocities.Component(axis).Valid_Index(boundary_face)) fine_face_velocities(FACE_INDEX<TV::dimension>(axis,boundary_face))=0;}}
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
        //if(Source_Box_Lazy_Inside(local_iterator.Location())){
        //    psi_N.Component(local_iterator.Axis())(local_iterator.Face_Index())=true;
        //    if(local_iterator.Axis()==2)local_face_velocities.Component(local_iterator.Axis())(local_iterator.Face_Index())=.5;
        //    else local_face_velocities.Component(local_iterator.Axis())(local_iterator.Face_Index())=0;}
        if(use_collisions){
            int first_cell_in_solid=rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->Implicit_Geometry_Lazy_Inside(local_iterator.First_Cell_Center()),second_cell_in_solid=rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->Implicit_Geometry_Lazy_Inside(local_iterator.Second_Cell_Center());
            if(first_cell_in_solid+second_cell_in_solid==1 || (test_number==3 && first_cell_in_solid+second_cell_in_solid==2)){
                psi_N.Component(local_iterator.Axis())(local_iterator.Face_Index())=true;
                TV rigid_velocity=RIGID_GEOMETRY<TV>::Pointwise_Object_Velocity(rigid_geometry_collection.particles.twist(rigid_particle_id),rigid_geometry_collection.particles.frame(rigid_particle_id).t,local_iterator.Location());
                TV fluid_velocity;
                fluid_velocity(local_iterator.Axis())=local_face_velocities.Component(local_iterator.Axis())(local_iterator.Face_Index());
                for(int i=1;i<TV::dimension;i++){int axis=(local_iterator.Axis()-1+i)%TV::dimension+1;
                    fluid_velocity(axis)=(
                        local_face_velocities.Component(axis)(local_iterator.Face_Index())
                        +local_face_velocities.Component(axis)(local_iterator.Face_Index()-TV_INT::Axis_Vector(local_iterator.Axis()))
                        +local_face_velocities.Component(axis)(local_iterator.Face_Index()+TV_INT::Axis_Vector(axis))
                        +local_face_velocities.Component(axis)(local_iterator.Face_Index()+TV_INT::Axis_Vector(axis)-TV_INT::Axis_Vector(local_iterator.Axis())))*0.25;}
                TV rigid_normal=rigid_geometry_collection.Rigid_Geometry(rigid_particle_id).Implicit_Geometry_Normal(local_iterator.Location());
                T projected_velocity=Dot_Product(rigid_normal,fluid_velocity-rigid_velocity);
                if(projected_velocity<0) local_face_velocities.Component(local_iterator.Axis())(local_iterator.Face_Index())=(fluid_velocity-rigid_normal*projected_velocity)(local_iterator.Axis());}
            else if(first_cell_in_solid+second_cell_in_solid==2 && test_number==2){
                psi_N.Component(local_iterator.Axis())(local_iterator.Face_Index())=true;
                fine_face_velocities.Component(local_iterator.Axis())(local_iterator.Face_Index())=RIGID_GEOMETRY<TV>::Pointwise_Object_Velocity(rigid_geometry_collection.particles.twist(rigid_particle_id),rigid_geometry_collection.particles.frame(rigid_particle_id).t,local_iterator.Location())(local_iterator.Axis());}}}}

    void Set_Boundary_Conditions(const T time)
    {fine_psi_N.Fill(false);projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    Set_Fine_Boundary_Conditions(fine_mac_grid,fine_face_velocities,fine_psi_N);
    Set_Coarse_Boundary_Conditions();
    projection.poisson->beta_given_on_faces=true;}

    template<class TV2>
    bool Source_Box_Lazy_Inside(const TV2& location){PHYSBAM_NOT_IMPLEMENTED();return false;}

    template<class T2>
    bool Source_Box_Lazy_Inside(const VECTOR<T2,2>& location)
    {return source_box_2d.Lazy_Inside(location);}

    template<class T2>
    bool Source_Box_Lazy_Inside(const VECTOR<T2,3>& location)
    {return source_box_3d.Lazy_Inside(location);}

    void Initialize_Fields()
    {for(typename GRID<TV>::FACE_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()) fine_face_velocities(iterator.Full_Index())=0;
    for(typename GRID<TV>::CELL_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()) density(iterator.Cell_Index())=0;}

    void Get_Scalar_Field_Sources(const T time)
    {for(typename GRID<TV>::CELL_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next())
        if(Source_Box_Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=1;
    if(use_collisions){
        for(typename GRID<TV>::CELL_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next())
            if(rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->Implicit_Geometry_Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=0;}}
    
    void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::dimension> >& force,const ARRAY<T,TV_INT>& density_ghost,const T dt,const T time)
    {
        T density_buoyancy_constant=T(2)/buoyancy_clamp;
        for(typename GRID<TV>::FACE_ITERATOR iterator(fine_mac_grid,0,GRID<TV>::WHOLE_REGION,-1,1);iterator.Valid();iterator.Next()){ // y-direction forces only
            T face_density=min((density_ghost(iterator.First_Cell_Index())+density_ghost(iterator.Second_Cell_Index()))*T(.5),buoyancy_clamp);
            T density_difference=face_density;
            if(density_difference>0) force.Component(1)(iterator.Face_Index())=density_buoyancy_constant*density_difference;}
    }

    bool Map_Fine_To_Local_Boundaries_For_Cell(GRID<TV>& local_mac_grid,ARRAY<bool,FACE_INDEX<TV::dimension> >& local_psi_N,TV_INT cell_index)
    {bool has_solids=false;
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
        local_psi_N(local_iterator.Full_Index())=fine_psi_N(FACE_INDEX<TV::dimension>(local_iterator.Axis(),cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector()));
        if(local_psi_N(local_iterator.Full_Index())) has_solids=true;}
    return has_solids;}

    void Map_Coarse_To_Fine(GRID<TV>& local_mac_grid)
    {for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
        ARRAY<T,FACE_INDEX<TV::dimension> > sum(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()));sum.Fill(0);
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
            T total_area=TV::dimension==2?sub_scale:sub_scale*sub_scale;
            if(alpha!=1){TV_INT offset;offset(local_iterator.Axis())++;sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),local_iterator.First_Boundary()?TV_INT::All_Ones_Vector():TV_INT::All_Ones_Vector()+offset);}
            //if(area!=1) std::cout<<"Area "<<area<<" product "<<product<<" beta "<<beta<<" weights "<<weights(index)<<" fine save "<<fine_face_velocities_save(fine_index)<<" delta "<<coarse_face_velocities(coarse_index)-coarse_face_velocities_save(coarse_index)<<std::endl;
            if(!fine_psi_N(fine_index))
                fine_face_velocities(fine_index)=(alpha!=1?((1-alpha)*one_over_area*(coarse_face_velocities(coarse_index)-sum(sum_index)/total_area)):0)+
                    (alpha?(alpha*(fine_face_velocities_save(fine_index)+weights(index)*product*one_over_area*beta*(coarse_face_velocities(coarse_index)-coarse_face_velocities_save(coarse_index)))):0);}}}

    void Map_Fine_To_Local_Interior_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,TV_INT coarse_cell_index,bool zero_out)
    {for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
        FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector());
        if(zero_out) local_face_velocities(local_iterator.Full_Index())=0;
        else local_face_velocities(local_iterator.Full_Index())=fine_face_velocities_save(fine_index);}}

    void Map_Fine_To_Local_Boundary_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,TV_INT coarse_cell_index)
    {for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
        FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector());
        local_face_velocities(local_iterator.Full_Index())=fine_face_velocities(fine_index);}}

    void Map_Coarse_To_Local_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,
        TV_INT coarse_cell_index,bool zero_out)
    {ARRAY<T,FACE_INDEX<TV::dimension> > sum(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()));sum.Fill(0);
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
        FACE_INDEX<TV::dimension> fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector());
        if(zero_out) local_face_velocities(local_iterator.Full_Index())=0;
        else local_face_velocities(local_iterator.Full_Index())=fine_face_velocities_save(fine_index);}
    FACE_INDEX<TV::dimension> fine_index,sum_index=FACE_INDEX<TV::dimension>(1,TV_INT::All_Ones_Vector());
    if(alpha!=1){
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){TV_INT offset;offset(local_iterator.Axis())++;
            fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector());
            sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),(local_iterator.First_Boundary())?(TV_INT::All_Ones_Vector()):(TV_INT::All_Ones_Vector()+offset));
            if(fine_psi_N(fine_index)) sum(sum_index)+=fine_face_velocities_save(fine_index);}}
    //update based on (1-a)*((1/A)*(Vco-avg Vs))+a*(Vfo+(Vcn-Vco)*(1/A))
    for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
        TV_INT2 index;int offset=0;for(int i=0;i<TV::dimension;i++){if(i==local_iterator.Axis()){offset++;continue;}index(i-offset)=local_iterator.Face_Index()(i);}
        FACE_INDEX<TV::dimension> coarse_index,fine_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector());
        if(local_iterator.First_Boundary()) coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_mac_grid.First_Face_Index_In_Cell(local_iterator.Axis(),coarse_cell_index));
        else coarse_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),coarse_mac_grid.Second_Face_Index_In_Cell(local_iterator.Axis(),coarse_cell_index));
        T product=1;for(int i=1;i<TV::dimension;i++) product*=sub_scale;
        T area=projection.poisson->beta_face(coarse_index),one_over_area=1./area;
        T beta=fine_psi_N(fine_index)?0:1;
        T total_area=TV::dimension==2?sub_scale:sub_scale*sub_scale;
        if(alpha!=1){TV_INT offset;offset(local_iterator.Axis())++;sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),local_iterator.First_Boundary()?TV_INT::All_Ones_Vector():TV_INT::All_Ones_Vector()+offset);}
        if(!fine_psi_N(fine_index))
            local_face_velocities(local_iterator.Full_Index())=(alpha!=1?((1-alpha)*one_over_area*(coarse_face_velocities(coarse_index)-sum(sum_index)/total_area)):0)+
                (alpha?(alpha*(fine_face_velocities_save(fine_index)+weights(index)*product*one_over_area*beta*(coarse_face_velocities(coarse_index)-coarse_face_velocities_save(coarse_index)))):0);}}

    void Map_Local_To_Fine_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,TV_INT cell_index)
    {for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
        fine_face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),cell_index*sub_scale+local_iterator.Face_Index()-sub_scale*TV_INT::All_Ones_Vector()))=local_face_velocities(local_iterator.Full_Index());}}

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

    void Preprocess_Frame(const int frame)
    {
        if(short_kol>0){
            if(frame==short_kol) kolmogorov=0;}
    }

    void Preprocess_Projection(const T dt,const T time)
    {
        Average_Velocities_From_Fine_To_Coarse(coarse_mac_grid,coarse_face_velocities,fine_face_velocities,sub_scale);
        Save_Velocities(coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_save);
        Save_Velocities(fine_mac_grid,fine_face_velocities,fine_face_velocities_save);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("preprocess example",0,1);
    }
    
    void Postprocess_Projection(const T dt,const T time)
    {
        PHYSBAM_DEBUG_WRITE_SUBSTEP("postprocess example",0,1);
        if(binary_refinement_levels){
            if(use_direct_solver) Local_Projection_Analytic(dt,time,binary_refinement_levels,coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_save,projection.poisson->beta_face,TV_INT::All_Ones_Vector());
            else Local_Projection_Analytic2(dt,time,binary_refinement_levels,coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_save,projection.poisson->beta_face,TV_INT::All_Ones_Vector());}
        else if(number_of_threads) Threaded_Local_Projection_PCG(dt,time);
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
        /*{
            int number_of_unknowns=matrix_index_to_cell_index.m;            
            A.Negate();b*=(T)-1;
            VECTOR_ND<T> x(number_of_unknowns),q,s,r,k,z;
            for(int i=0;i<number_of_unknowns;i++) x(i)=local_projection.elliptic_solver->u(matrix_index_to_cell_index(i));
            local_projection.elliptic_solver->Find_Tolerance(b);
            local_projection.elliptic_solver->pcg.Solve(A,x,b,q,s,r,k,z,local_projection.elliptic_solver->tolerance,false);
            for(int i=0;i<number_of_unknowns;i++)local_projection.elliptic_solver->u(matrix_index_to_cell_index(i))=x(i);
        }*/
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
        Scaled_Direct_Solvers(local_projection.p,local_projection.elliptic_solver->f,scale);
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::INTERIOR_REGION);local_iterator.Valid();local_iterator.Next()){
            TV_INT first_cell_index=local_iterator.First_Cell_Index(),second_cell_index=local_iterator.Second_Cell_Index();
            face_velocities.Component(local_iterator.Axis())(local_iterator.Face_Index())-=(local_projection.p(second_cell_index)-local_projection.p(first_cell_index));}
    }

    void Threaded_Local_Projection_PCG(const T dt,const T time)
    {
        Map_Coarse_To_Fine(local_mac_grid);
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            REFINEMENT_TASK<TV> *refinement_task=new REFINEMENT_TASK<TV>();
            refinement_task->smoke_tests=this;refinement_task->cell_index=iterator.Cell_Index();refinement_task->dt=dt;refinement_task->time=time;
            thread_queue->Queue(refinement_task);}
        thread_queue->Wait();
    }

    void Local_Projection_PCG(const T dt,const T time)
    {
        Map_Coarse_To_Fine(local_mac_grid);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after mapping to fine cell",0,1);
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
            //Map_Coarse_To_Local_For_Cell(local_mac_grid,local_face_velocities,coarse_face_velocities,iterator.Cell_Index(),false);
            Map_Fine_To_Local_Boundary_For_Cell(local_mac_grid,local_face_velocities,iterator.Cell_Index());
            Map_Fine_To_Local_Interior_For_Cell(local_mac_grid,local_face_velocities,iterator.Cell_Index(),false);
            bool contains_solids=Map_Fine_To_Local_Boundaries_For_Cell(local_mac_grid,local_projection.elliptic_solver->psi_N,iterator.Cell_Index());
            for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
                local_projection.p(local_iterator.Cell_Index())=projection.p(iterator.Cell_Index());}
            local_projection.elliptic_solver->Set_Neumann_Outer_Boundaries();
            local_projection.p*=dt;        
            if(contains_solids) local_projection.Make_Divergence_Free(local_face_velocities,dt,time);
            else if(use_direct_solver) Fast_Direct_Solve(local_projection,local_face_velocities,dt,time,sub_scale);
            else Fast_Unwrapped_PCG_Solve(local_projection,local_face_velocities,dt,time);
            //else local_projection.Make_Divergence_Free(local_face_velocities,dt,time);
            local_projection.p/=dt;
            //Map_Local_To_Fine_For_Cell(local_mac_grid,local_face_velocities,iterator.Cell_Index());}
            Map_Local_To_Fine_Interior_For_Cell(local_mac_grid,local_face_velocities,iterator.Cell_Index());}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after postprocess cell",0,1);
    }

    void Initialize_Bodies()
    {
        std::string model_file_name;
        T radius;
        switch(test_number){
            case 1:
                break;
            case 2:
                use_collisions=true;
                if(TV::dimension==2) model_file_name="../../Public_Data/Rigid_Bodies_2D/circle";
                else model_file_name="../../Public_Data/Rigid_Bodies/sphere";
                radius=0.07;
                rigid_particle_id=rigid_geometry_collection.Add_Rigid_Geometry(stream_type,model_file_name,radius,true,true,true,true);
                rigid_geometry_collection.particles.frame(rigid_particle_id).t=TV::Constant_Vector(0.25);
                rigid_geometry_collection.particles.frame(rigid_particle_id).t(1)=0.5;
                rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->is_static=true;
                break;
            case 3:
                use_collisions=true;
                if(TV::dimension==2) model_file_name="../../Public_Data/Rigid_Bodies_2D/circle";
                else model_file_name="../../Public_Data/Rigid_Bodies/sphere";
                radius=0.1;
                rigid_particle_id=rigid_geometry_collection.Add_Rigid_Geometry(stream_type,model_file_name,radius,true,true,true,true);
                rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->is_static=false;
                break;}
        rigid_geometry_collection.Update_Kinematic_Particles();
    }

    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
    {
        if(test_number==3){
            frame.t=TV::Constant_Vector(0.25);frame.t(1)=0.5;
            if(time<4.0) frame.t(1)=sin(2.*pi*time/4.)*0.2+0.5;
            else frame.t(0)=sin(2*pi*time/4.)*0.1+0.25;}
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
        ARRAY<T,FACE_INDEX<TV::dimension> > sum(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT::All_Ones_Vector()));sum.Fill(0);
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
            T total_area=TV::dimension==2?sub_scale:sub_scale*sub_scale;
            if(alpha!=1){TV_INT offset;offset(local_iterator.Axis())++;sum_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),local_iterator.First_Boundary()?TV_INT::All_Ones_Vector():TV_INT::All_Ones_Vector()+offset);}
            if(!local_psi_N.Component(local_iterator.Axis())(local_iterator.Face_Index()) && area!=0) 
                local_face_velocities(local_iterator.Full_Index())=(alpha!=1?((1-alpha)*one_over_area*(coarse_face_velocities(coarse_index)-sum(sum_index)/total_area)):0)+
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
