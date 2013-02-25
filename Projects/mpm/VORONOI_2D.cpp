//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include "VORONOI_2D.h"
namespace PhysBAM{
//#####################################################################
// Function Initialize_With_A_Regular_Grid
//#####################################################################
template<class T> void VORONOI_2D<T>::
Initialize_With_A_Regular_Grid(const GRID<TV>& grid)
{
    HASHTABLE<TV_INT,int> particle_index;
    int ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+grid.counts));it.Valid();it.Next()) particle_index.Get_Or_Insert(it.index)=ID++;

    HASHTABLE<TV_INT,int> segment_node_index;
    ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+grid.counts+1));it.Valid();it.Next()) segment_node_index.Get_Or_Insert(it.index)=ID++;
    int N_segmented_mesh_nodes=(grid.counts+TV_INT(1,1)).Product();

    RANGE<TV> segmented_mesh_domain(grid.domain.min_corner-(T)0.5*grid.dX,grid.domain.max_corner+(T)0.5*grid.dX);
    GRID<TV> segmented_mesh_grid(grid.counts+TV_INT(1,1),segmented_mesh_domain);

    segment_mesh_particles.Resize(N_segmented_mesh_nodes);
    segment_mesh.Set_Number_Nodes(N_segmented_mesh_nodes);
    voronoi_mesh_of_particle.Resize(grid.counts.Product());

    ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+segmented_mesh_grid.counts));it.Valid();it.Next()){
        segment_mesh_particles.X(ID++)=segmented_mesh_grid.Node(it.index);
        if(it.index(0)+1<segmented_mesh_grid.counts(0)){
            TV_INT next_index(it.index(0)+1,it.index(1));
            segment_mesh.elements.Append(TV_INT(segment_node_index.Get_Or_Insert(it.index),segment_node_index.Get_Or_Insert(next_index)));
        }
        if(it.index(1)+1<segmented_mesh_grid.counts(1)){
            TV_INT next_index(it.index(0),it.index(1)+1);
            segment_mesh.elements.Append(TV_INT(segment_node_index.Get_Or_Insert(it.index),segment_node_index.Get_Or_Insert(next_index)));
        }
    }
    
    for(int i=0;i<segment_mesh_particles.number;i++){
        segment_mesh_particles.X(i)=clamp(segment_mesh_particles.X(i),grid.domain.min_corner,grid.domain.max_corner);
    }

    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+grid.counts));it.Valid();it.Next()){
        TV_INT ll=segmented_mesh_grid.Cell(grid.Node(it.index),0);
        voronoi_mesh_of_particle(particle_index.Get_Or_Insert(it.index)).Append(TV_INT(segment_node_index.Get_Or_Insert(ll),segment_node_index.Get_Or_Insert(ll+TV_INT(1,0))));
        voronoi_mesh_of_particle(particle_index.Get_Or_Insert(it.index)).Append(TV_INT(segment_node_index.Get_Or_Insert(ll),segment_node_index.Get_Or_Insert(ll+TV_INT(0,1))));
        voronoi_mesh_of_particle(particle_index.Get_Or_Insert(it.index)).Append(TV_INT(segment_node_index.Get_Or_Insert(ll+TV_INT(1,0)),segment_node_index.Get_Or_Insert(ll+TV_INT(1,1))));
        voronoi_mesh_of_particle(particle_index.Get_Or_Insert(it.index)).Append(TV_INT(segment_node_index.Get_Or_Insert(ll+TV_INT(0,1)),segment_node_index.Get_Or_Insert(ll+TV_INT(1,1))));
    }
    

}
//#####################################################################
template class VORONOI_2D<float>;
template class VORONOI_2D<double>;
}
