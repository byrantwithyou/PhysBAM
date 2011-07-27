//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONVERSION_TOOLS
//#####################################################################
#ifndef __AEROF_FLUID_TO_UNIFORM_GRID__
#define __AEROF_FLUID_TO_UNIFORM_GRID__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Particles/COMPRESSIBLE_FLUID_PARTICLES.h>
#include <Deformable_Objects/DEFORMABLE_OBJECT.h>

namespace PhysBAM{
namespace CONVERSION_TOOLS{

template<class T,int d=3> void
Convert_To_Uniform_Grid(const TETRAHEDRALIZED_VOLUME<T>& ale_fluid_mesh,const COMPRESSIBLE_FLUID_PARTICLES<VECTOR<T,d> >& ale_fluid_data,const GRID<VECTOR<T,d> >& grid,
                        const typename GRID_ARRAYS_POLICY<GRID<VECTOR<T,d> > >::ARRAYS_SCALAR& phi,
                        typename GRID_ARRAYS_POLICY<GRID<VECTOR<T,d> > >::ARRAYS_SCALAR::template REBIND<VECTOR<VECTOR<T,d+2>,4> >::TYPE& U,const T conversion_thickness=(T)FLT_MAX)
{
    typedef typename GRID_ARRAYS_POLICY<GRID<VECTOR<T,d> > >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<VECTOR<VECTOR<T,d+2>,4> >::TYPE T_ARRAYS_STORED_DIMENSION;

    T_ARRAYS_BOOL done(grid); done.Fill(false);
    T_ARRAYS_INT bounding_tet_id(grid);bounding_tet_id.Fill(0);
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,d> >,T> interpolation;

    LOG::Time("Compute which ale fluid nodes have valid data"); // TODO(jontg): Is this data already available somewhere?

    ARRAY<bool> particle_not_inside_solid(tet_volume.particles.array_collection->Size());
    for(int i=1;i<=particle_not_inside_solid.m;++i)
        particle_not_inside_solid(i)=interpolation.Clamped_To_Array(grid,phi,tet_volume.particles.X(i))<0;

    LOG::Time("Computing bounding tet id");

    T phi_correction_term=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(current_grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();TV position=current_grid.X(cell_index);T phi_value=phi(cell_index);
        if(phi_value > 0 && phi_value < conversion_thickness){
            ARRAY<int> intersection_list;ale_fluid_mesh.hierarchy->Intersection_List(position,intersection_list);
            for(int i=1;i<=intersection_list.m;++i){
                TETRAHEDRON<T>& current_tet=(*ale_fluid_mesh.tetrahedron_list)(intersection_list(i));
                VECTOR<int,4>& current_tet_nodes=ale_fluid_mesh.mesh.elements(intersection_list(i));
                bool tet_has_no_ghost_nodes=true;for(int j=1;j<=4;++j) tet_has_no_ghost_nodes&=particle_not_inside_solid(current_tet_nodes(j));
                if(tet_has_no_ghost_nodes && current_tet.Inside(position)){bounding_tet_id(cell_index)=intersection_list(i);break;}
                else phi_correction_term=max(phi_correction_term,phi_value);}}}

    LOG::cout<<"\t\tphi correction term = "<<phi_with_correction_term<<std::endl;
    LOG::Time("Interpolating particle data to the Eulerian grid");

    T min_phi=fluid_extrapolation_bandwidth;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        TV position=grid.X(cell_index);T phi_value=phi(cell_index);
        if(bounding_tet_id(cell_index)){
            TETRAHEDRON<T>& current_tet=(*ale_fluid_mesh.tetrahedron_list)(bounding_tet_id(cell_index));
            VECTOR<int,4>& corners=ale_fluid_mesh.mesh.elements(bounding_tet_id(cell_index));
            TV weights=current_tet.Barycentric_Coordinates(position);

            U(cell_index)(1)=Point_From_Barycentric_Coordinates(weights,ale_fluid_data.rho(corners(1)),ale_fluid_data.rho(corners(2)),ale_fluid_data.rho(corners(3)),ale_fluid_data.rho(corners(4)));
            for(int j=1;j<=d;++j)
                U(cell_index)(j+1)=Point_From_Barycentric_Coordinates(weights,ale_fluid_data.V(corners(1))(j),ale_fluid_data.V(corners(2))(j),ale_fluid_data.V(corners(3))(j),ale_fluid_data.V(corners(4))(j));
            U(cell_index)(d+2)=Point_From_Barycentric_Coordinates(weights,ale_fluid_data.E(corners(1)),ale_fluid_data.E(corners(2)),ale_fluid_data.E(corners(3)),ale_fluid_data.E(corners(4)));
            done(cell_index)=true;}
        else if(phi_value > 0 && phi_value > conversion_thickness) done(cell_index)=true;}

    LOG::Time("Extrapolating state into solid");
    
    T_ARRAYS_SCALAR phi_with_correction_term(phi);phi_with_correction_term-=phi_correction_term-1e-16;
    EXTRAPOLATION_UNIFORM<GRID<VECTOR<T,d> >,VECTOR<T,d+2> > extrapolation(grid,phi_with_correction_term,U,3);
    
    extrapolation.Set_Custom_Seed_Done(&done);
    extrapolation.Extrapolate((T)0,false);

    LOG::Stop_Time();
}
}
#endif
