//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_TESTS__
#define __SMOKE_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_ADAPTIVE_EXAMPLE.h>

namespace PhysBAM{

template<class TV>
class SMOKE_TESTS:public INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<TV> BASE;
    using BASE::mac_grid;using BASE::face_velocities;using BASE::write_substeps_level;using BASE::incompressible;using BASE::projection;using BASE::output_directory;
    using BASE::density;using BASE::temperature;

    bool use_collisions;
    int binary_refinement_levels;
    ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> >,TV_INT> weights;
    T alpha;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_save;
    RANGE<TV> source_box;
    SPHERE<TV> collision_sphere;
public:

    SMOKE_TESTS(const STREAM_TYPE stream_type,const PARSE_ARGS& parse_args)
        :INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<TV>(stream_type),source_box((T).25*TV::All_Ones_Vector(),(T).75*TV::All_Ones_Vector()),collision_sphere((T).5*TV::All_Ones_Vector(),.2)
    {
        use_collisions=false;
        int test_number=1;
        alpha=parse_args.Get_Double_Value("-alpha");        
        write_substeps_level=parse_args.Get_Integer_Value("-substep");
        int scale=parse_args.Get_Integer_Value("-scale");
        bool binary_refinement=parse_args.Get_Option_Value("-binary");
        int sub_scale=parse_args.Get_Integer_Value("-subscale");
        binary_refinement_levels=0;
        if(binary_refinement){int tmp=sub_scale;while(tmp>1){assert((tmp&1)==0);tmp>>=1;binary_refinement_levels++;}}
        output_directory=STRING_UTILITIES::string_sprintf("Smoke_Tests/Test_%d_%d_%d%s_%1.2f",test_number,scale,sub_scale,binary_refinement_levels?"_binary":"",alpha);        
        source_box.min_corner.y=(T)0;source_box.max_corner.y=(T)0.25;
        mac_grid.Initialize(TV_INT::All_Ones_Vector()*scale,RANGE<TV>(TV(),TV::All_Ones_Vector()),true);
        face_velocities_save.Resize(mac_grid);
        weights.Resize(mac_grid.Domain_Indices(),false);
        //resize/initialize weights
        for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            TV local_offset=TV::All_Ones_Vector()*.5*mac_grid.DX();
            mac_grid.sub_mac_grids(iterator.Cell_Index())=new GRID<TV>(TV_INT::All_Ones_Vector()*sub_scale,RANGE<TV>(iterator.Location()-local_offset,iterator.Location()+local_offset),true);
            GRID<TV>& local_mac_grid=*mac_grid.sub_mac_grids(iterator.Cell_Index());
            ARRAY<T,FACE_INDEX<TV::dimension> >& local_weights=weights(iterator.Cell_Index());local_weights.Resize(local_mac_grid);
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
                TV_INT counts=local_mac_grid.Counts();counts(local_iterator.Axis())=1;local_weights(local_iterator.Full_Index())=1/counts.Product();}}
    }

    void Write_Output_Files(const int frame)
    {BASE::Write_Output_Files(frame);}

    void Set_Boundary_Conditions(const T time,const GRID<TV>& mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >* projection)
    {if(projection){projection->elliptic_solver->psi_D.Fill(false);
        for(int axis=1;axis<=TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*(axis-1)+axis_side;
            //if(!mpi_grid || domain_boundary(axis)(axis_side)){ //Need to check mpi as smaller solves are never using mpi (for now)
            if(true){ //No mpi yet
                TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);    
                for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                    projection->elliptic_solver->psi_D(cell)=true;projection->p(cell)=0;}}}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
        if(source_box.Lazy_Inside(iterator.Location())){
            if(projection) projection->elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==2)face_velocities(iterator.Full_Index())=1;
            else face_velocities(iterator.Full_Index())=0;}
        if(use_collisions && collision_sphere.Lazy_Inside(iterator.Location())){
            if(projection) projection->elliptic_solver->psi_N(iterator.Full_Index())=true;
            face_velocities(iterator.Full_Index())=0;}}}

    void Initialize_Fields()
    {for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
        GRID<TV>& local_mac_grid=*mac_grid.sub_mac_grids(iterator.Cell_Index());ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities=*face_velocities.sub_arrays(iterator.Cell_Index());
        ARRAY<T,TV_INT>& local_density=density(iterator.Cell_Index());ARRAY<T,TV_INT>& local_temperature=temperature(iterator.Cell_Index());
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()) local_face_velocities(local_iterator.Full_Index())=0;
        for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()) local_density(local_iterator.Cell_Index())=local_temperature(local_iterator.Cell_Index())=0;}}

    void Get_Scalar_Field_Sources(const T time,const GRID<TV>& mac_grid,ARRAY<T,TV_INT>& density)
    {for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
        if(source_box.Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=1;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
        if(collision_sphere.Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=0;}

    void Map_Coarse_To_Local_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,typename GRID<TV>::CELL_ITERATOR& iterator)
    {for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
        ARRAY<T,FACE_INDEX<TV::dimension> >& local_weights=weights(iterator.Cell_Index());
        FACE_INDEX<TV::dimension> grid_index;TV_INT counts=local_mac_grid.Counts();counts(local_iterator.Axis())=1;
        if(local_iterator.First_Boundary()) grid_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.First_Face_Index(local_iterator.Axis()));
        else grid_index=FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Second_Face_Index(local_iterator.Axis()));
        local_face_velocities(local_iterator.Full_Index())=
            (1-alpha)*face_velocities(grid_index)+alpha*(local_face_velocities(local_iterator.Full_Index())+local_weights(local_iterator.Full_Index())*counts.Product()*(face_velocities(grid_index)-face_velocities_save(grid_index)));}}

    void Save_Velocities()
    {for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities_save(iterator.Full_Index())=face_velocities(iterator.Full_Index());}

    void Preprocess_Projection(const T dt,const T time)
    {
        for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
        for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            GRID<TV>& local_mac_grid=*mac_grid.sub_mac_grids(iterator.Cell_Index());ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities=*face_velocities.sub_arrays(iterator.Cell_Index());
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,0,GRID<TV>::BOUNDARY_REGION);local_iterator.Valid();local_iterator.Next()){
                if(local_iterator.First_Boundary()) face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.First_Face_Index(local_iterator.Axis())))+=local_face_velocities(local_iterator.Full_Index());
                else face_velocities(FACE_INDEX<TV::dimension>(local_iterator.Axis(),iterator.Second_Face_Index(local_iterator.Axis())))+=local_face_velocities(local_iterator.Full_Index());}}    
        for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            TV_INT counts;T factor=0;
            if(iterator.First_Cell_Index()(iterator.Axis())>0){counts=mac_grid.sub_mac_grids(iterator.First_Cell_Index())->Counts();counts(iterator.Axis())=1;factor=counts.Product();}
            if(iterator.Second_Cell_Index()(iterator.Axis())<=mac_grid.Counts()(iterator.Axis())){counts=mac_grid.sub_mac_grids(iterator.Second_Cell_Index())->Counts();counts(iterator.Axis())=1;factor+=counts.Product();}
            face_velocities(iterator.Full_Index())/=(T)factor;}
        Save_Velocities();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("preprocess example",0,1);
    }
 
    void Postprocess_Projection(const T dt,const T time)
    {
        PHYSBAM_DEBUG_WRITE_SUBSTEP("postprocess example",0,1);
        //if(binary_refinement_levels) Local_Projection_Analytic(dt,time,binary_refinement_levels,coarse_mac_grid,coarse_face_velocities,TV_INT::All_Ones_Vector());
        //else Local_Projection_PCG(dt,time);
        Local_Projection_PCG(dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after local projection example",0,1);
    }

    void Local_Projection_PCG(const T dt,const T time)
    {
        for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            GRID<TV>& local_mac_grid=*mac_grid.sub_mac_grids(iterator.Cell_Index());ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities=*face_velocities.sub_arrays(iterator.Cell_Index());
            PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > local_projection(local_mac_grid);
            INCOMPRESSIBLE_UNIFORM<GRID<TV> > local_incompressible(local_mac_grid,local_projection);
            local_incompressible.Initialize_Grids(local_mac_grid);
            local_projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
            local_projection.elliptic_solver->pcg.Use_Modified_Incomplete_Cholesky();
            VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
            local_incompressible.boundary->Set_Constant_Extrapolation(constant_extrapolation);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("before mapping cell",0,1);
            Map_Coarse_To_Local_For_Cell(local_mac_grid,local_face_velocities,face_velocities,iterator);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after mapping cell",0,1);
            for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next()){
                local_projection.p(local_iterator.Cell_Index())=projection.p(iterator.Cell_Index());}
            Set_Boundary_Conditions(time,local_mac_grid,local_face_velocities,&local_projection);
            local_projection.elliptic_solver->Set_Neumann_Outer_Boundaries();
            local_projection.p*=dt;        
            local_incompressible.Advance_One_Time_Step_Implicit_Part(local_face_velocities,dt,time,false);
            PHYSBAM_DEBUG_WRITE_SUBSTEP("after projection",0,1);
            local_projection.p/=dt;}          
    }

//#####################################################################
};
}
#endif
