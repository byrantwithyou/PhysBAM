//#####################################################################
// Copyright 2007-2008, Jon Gretarsson, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AERO_INTERFACE_1
//##################################################################### 
#ifndef __AERO_INTERFACE_1__
#define __AERO_INTERFACE_1__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_CHARBEL.h>
#include <PhysBAM_Dynamics/Particles/COMPRESSIBLE_FLUID_PARTICLES.h>
#include <PhysBAM_Tools/Utilities/scoped_ptr.h>

namespace PhysBAM{

template<class T>
class AERO_INTERFACE_1{
    typedef T RW;
    typedef GRID<TV> T_GRID;
    typedef VECTOR<T,3> TV;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<VECTOR<VECTOR<T,T_GRID::dimension+2>,4> >::TYPE T_ARRAYS_STORED_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND_LENGTH<T_GRID::dimension+2>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_DIMENSION_SCALAR::ELEMENT T_ARRAYS_ELEMENT;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR::template REBIND<TV>::TYPE T_LINEAR_INTERPOLATION_VECTOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR::template REBIND<T_ARRAYS_ELEMENT>::TYPE T_LINEAR_INTERPOLATION_DIMENSION;

private:
    STREAM_TYPE stream_type;
    T_GRID grid;
    T_GRID* global_grid;
    EULER_UNIFORM<T_GRID> euler;
    RIGID_BODY_PARTICLES<TV> rigid_body_particles;
    RIGID_BODY<TV>* rigid_body;
    scoped_ptr<FLUID_COLLISION_BODY_LIST_UNIFORM<T_GRID> > collision_body_list;
    FRAME<TV> levelset_frame;
    scoped_ptr<FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID> > inaccurate_union;
    T_ARRAYS_SCALAR phi_all_solids;
    T_ARRAYS_INT bounding_tet_id;
    T_ARRAYS_VECTOR gradient;
    T_LINEAR_INTERPOLATION_VECTOR interpolation;
    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;
    MPI_CHARBEL<T>* mpi_charbel;
    BOUNDARY_MPI<T_GRID,typename T_GRID::SCALAR,T_GRID::dimension+2>* mpi_boundary;
    BOUNDARY_UNIFORM<T_GRID,bool> done_boundary;

    T solid_extrapolation_bandwidth;
    T fluid_extrapolation_bandwidth;

    int substeps_level;
    int first_frame,output_number;
    const std::string output_directory;
public:
    AERO_INTERFACE_1(T_GRID domain,const std::string& output_directory_input="PhysBAM_output",int substeps_input=-1);
    ~AERO_INTERFACE_1();

    void Set_Substeps(const int substeps_input=-1)
    {substeps_level=substeps_input;}

    FRAME<TV>& Frame() PHYSBAM_ALWAYS_INLINE
    {return levelset_frame;}

    void Write_Triangulated_Surface(const std::string& filename, const TRIANGULATED_SURFACE<T>& triangulated_surface)
    {FILE_UTILITIES::Write_To_File<T>(filename,triangulated_surface);}

    void Write_Tetrahedralized_Volume(const std::string& filename, const TETRAHEDRALIZED_VOLUME<T>& tet_volume)
    {FILE_UTILITIES::Write_To_File<T>(filename,tet_volume);}

    void Write_Levelset(const std::string& filename){
        PHYSBAM_ASSERT(inaccurate_union.get());
        LEVELSET_IMPLICIT_OBJECT<TV> levelset_implicit_object(grid,inaccurate_union->levelset.levelset.phi);
        FILE_UTILITIES::Write_To_File<T>(filename,levelset_implicit_object);
    }

    void Write_First_Frame(const int frame) const
    {if(frame==first_frame)FILE_UTILITIES::Write_To_Text_File(output_directory+"/first_frame",first_frame,"\n");}

    void Write_Last_Frame(const int frame) const
    {FILE_UTILITIES::Write_To_Text_File(output_directory+"/last_frame",frame,"\n");}

    void Write_Solid_Output_Files(const int frame) const
    {std::string prefix=output_directory+"/";rigid_body_particles.Write(stream_type,output_directory,frame);}

public:
//#####################################################################
    const T Phi(const TV& position) const;
    const TV Gradient(const TV& position) const;
    const bool Inside_Solid(const TV& position) const;
    const bool Cell_Should_Be_Populated(const TV& position,const T& depth) const;
    void Initialize_MPI(TETRAHEDRALIZED_VOLUME<T>& tet_volume,ARRAY<int>& local_to_global_map,const int& global_particle_count,T depth,TV* dim_min=0,TV* dim_max=0);
    void Compute_Level_Set(TRIANGULATED_SURFACE<T>& triangulated_surface,COMPRESSIBLE_FLUID_PARTICLES<TV>* particles_aerof=0);
    int Closest_Boundary_Point(const TV& location,TV& closest_bondary_point,T& distance);
    void Initialize_Acceleration_Structures(TETRAHEDRALIZED_VOLUME<T>& tet_volume,T depth);
    void Compute_Ghost_Cells(TETRAHEDRALIZED_VOLUME<T>& tet_volume,COMPRESSIBLE_FLUID_PARTICLES<TV>& particles);
    TV Compute_Lift();
    void Write_Fluid_Output_Files(const int frame);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int level);
//#####################################################################
};
}
#endif
