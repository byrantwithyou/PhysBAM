//#####################################################################
// Copyright 2005-2008, Geoffrey Irving, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET.h>
#include <Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#ifdef USE_MPI
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#endif
namespace PhysBAM{

#ifdef USE_MPI

//#####################################################################
// Function ISend_Particles
//#####################################################################
template<class TV,class T_PARTICLES> MPI::Request
ISend_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const ARRAY<PAIR<T_PARTICLES*,int> >& particles,const int destination_rank,
    const VECTOR<int,TV::m>& destination_direction,const int tag,ARRAY<char>& buffer)
{
    int position=0;
    MPI::Comm& comm=*mpi_grid.comm;
    buffer.Resize(MPI_UTILITIES::Pack_Size(destination_direction,comm)+MPI_UTILITIES::Pack_Size(particles.m,comm)
        +(particles.m?particles.m*MPI_UTILITIES::Pack_Size(*particles(0).x,comm):0));
    MPI_UTILITIES::Pack(destination_direction,buffer,position,comm);
    MPI_UTILITIES::Pack(particles.m,buffer,position,comm);
    for(int i=0;i<particles.m;i++) MPI_UTILITIES::Pack(*particles(i).x,particles(i).y,buffer,position,comm);
    return comm.Isend(buffer.Get_Array_Pointer(),position,MPI::PACKED,destination_rank,tag);
}
//#####################################################################
// Function Recv_Particles
//#####################################################################
template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Recv_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int tag,
    const MPI::Status& probe_status,PARTICLE_LEVELSET<TV>& particle_levelset)
{
    
    typedef VECTOR<int,TV::m> TV_INT;
    MPI::Comm& comm=*mpi_grid.comm;
    T_PARTICLES* recv_particles(template_particles.Clone()); // message will unpack into this particle object
    recv_particles->Add_Element();
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int i=0;i<m;i++){
        MPI_UTILITIES::Unpack(*recv_particles,1,buffer,position,comm);
        TV& X=recv_particles->X(0);X+=wrap_offset;
        if(!domain.Lazy_Inside(X)) continue;
        TV_INT final_b=mpi_grid.local_grid.Block_Index(X,1); // TODO: check whether this is a good block
        if(!particles(final_b)) particles(final_b)=particle_levelset.Allocate_Particles(template_particles);
        particle_levelset.Copy_Particle(*recv_particles,*particles(final_b),1);}
    delete recv_particles;
}
//#####################################################################
// Function Recv_Block_Particles
//#####################################################################
template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Recv_Block_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int tag,
    const MPI::Status& probe_status,PARTICLE_LEVELSET<TV>& particle_levelset)
{
    
    typedef VECTOR<int,TV::m> TV_INT;
    MPI::Comm& comm=*mpi_grid.comm;
    T_PARTICLES* recv_particles(template_particles.Clone()); // message will unpack into this particle object
    recv_particles->Add_Element();
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int i=0;i<m;i++){
        MPI_UTILITIES::Unpack(*recv_particles,1,buffer,position,comm);
        TV& X=recv_particles->X(0);X+=wrap_offset;
        if(domain.Lazy_Inside(X)) continue;
        TV_INT final_b=mpi_grid.local_grid.Block_Index(X,1); // TODO: check whether this is a good block
        if(!particles(final_b)) particles(final_b)=particle_levelset.Allocate_Particles(template_particles);
        particle_levelset.Copy_Particle(*recv_particles,*particles(final_b),1);}
    delete recv_particles;
}
//#####################################################################
// Function Recv_Ghost_Particles
//#####################################################################
template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Recv_Ghost_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int tag,
    const MPI::Status& probe_status,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset)
{
    
    typedef VECTOR<int,TV::m> TV_INT;
    MPI::Comm& comm=*mpi_grid.comm;
    T_PARTICLES* recv_particles(template_particles.Clone()); // message will unpack into this particle object
    recv_particles->Add_Element();
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int i=0;i<m;i++){
        MPI_UTILITIES::Unpack(*recv_particles,1,buffer,position,comm);
        TV& X=recv_particles->X(0);X+=wrap_offset;
        if(domain.Lazy_Inside(X)) continue;
        TV_INT final_b=mpi_grid.local_grid.Block_Index(X,bandwidth-1); // TODO: check whether this is a good block
        if(!particles(final_b)) particles(final_b)=particle_levelset.Allocate_Particles(template_particles);
        particle_levelset.Copy_Particle(*recv_particles,*particles(final_b),1);}
    delete recv_particles;
}
//#####################################################################
// Function ISend_Particles
//#####################################################################
template<class TV,class T_PARTICLES,class T> void
ISend_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,T_PARTICLES& particles,const T ghost_distance,const int tag,ARRAY<ARRAY<char> >& buffers,ARRAY<MPI::Request>& requests)
{
    typedef VECTOR<int,TV::m> TV_INT;
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    ARRAY<RANGE<TV> > neighbor_domains(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    for(int n=0;n<neighbor_domains.m;n++){
        TV neighbor_direction=TV(GRID<TV>::One_Ring_Neighbor(TV_INT(),n));
        neighbor_domains(n)=(domain+neighbor_direction*domain.Edge_Lengths()).Thickened((T).01*mpi_grid.local_grid.dX.Min());}
    if(ghost_distance){
        domain.Change_Size(-ghost_distance);
        for(int n=0;n<neighbor_domains.m;n++)neighbor_domains(n).Change_Size(ghost_distance);}
    // send particles that have exited the domain
    buffers.Resize(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and doesn't locally delete sent particles
    for(int n=0;n<GRID<TV>::number_of_one_ring_neighbors_per_cell;n++)if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        exchange_particles(n).Preallocate(100);
        for(int i=0;i<particles.Size();i++){
            if(!domain.Lazy_Inside(particles.X(i)) && neighbor_domains(n).Lazy_Inside(particles.X(i)))
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(&particles,i));}
        requests.Append(ISend_Particles(mpi_grid,exchange_particles(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}
}
//#####################################################################
// Function Recv_Particles
//#####################################################################
template<class TV,class T_PARTICLES> void
Recv_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,T_PARTICLES& particles,const int tag)
{
    typedef VECTOR<int,TV::m> TV_INT;
    MPI::Comm& comm=*mpi_grid.comm;
    MPI::Status probe_status;
    comm.Probe(MPI::ANY_SOURCE,tag,probe_status);
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    int number=particles.Size();particles.Add_Elements(m);
    for(int i=0;i<m;i++){
        MPI_UTILITIES::Unpack(particles,++number,buffer,position,comm);
        particles.X(number)+=wrap_offset;}
}
template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Exchange_Boundary_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset)
{
    if(mpi_grid.threaded_grid){Exchange_Boundary_Particles_Threaded(*mpi_grid.threaded_grid,template_particles,particles,bandwidth,particle_levelset);return;}
    typedef VECTOR<int,TV::m> TV_INT;
    STATIC_ASSERT((is_same<T_PARTICLES,typename remove_pointer<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE>::value));
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<RANGE<TV_INT> > send_regions;
    // this way Find_Boundary_Regions will return block indices which for uniform grids are the node index not minimum corner cell
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector());
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);
    // send particles that have exited the domain
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and send a lot more corner particles than it should
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int n=0;n<send_regions.m;n++) if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        exchange_particles(n).Preallocate(100);
        for(NODE_ITERATOR<TV> iterator(mpi_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next())if(particles(iterator.Node_Index())){
            T_PARTICLES* cell_particles=particles(iterator.Node_Index());
            while(cell_particles){for(int i=0;i<cell_particles->Size();i++)if(!domain.Lazy_Inside(cell_particles->X(i))) // could be optimized further
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(cell_particles,i));cell_particles=cell_particles->next;}} // TODO: delete the particle locally?
        requests.Append(ISend_Particles(mpi_grid,exchange_particles(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}
    // probe and receive
    for(int message=0;message<requests.m;message++){
        MPI::Status probe_status;
        mpi_grid.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Particles(mpi_grid,template_particles,particles,tag,probe_status,particle_levelset);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}
template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Exchange_Overlapping_Block_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset)
{
    if(mpi_grid.threaded_grid){Exchange_Overlapping_Block_Particles_Threaded(*mpi_grid.threaded_grid,template_particles,particles,bandwidth,particle_levelset);return;}
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    STATIC_ASSERT((is_same<T_PARTICLES,typename remove_pointer<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE>::value));
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<RANGE<TV_INT> > send_regions;
    // this way Find_Boundary_Regions will return block indices which for uniform grids are the node index not minimum corner cell
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector());
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);
    // send particles that have exited the domain
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and send a lot more corner particles than it should
    RANGE<TV> block_domain=mpi_grid.local_grid.Domain();block_domain.Change_Size(-(T).5*mpi_grid.local_grid.dX);
    for(int n=0;n<send_regions.m;n++) if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        exchange_particles(n).Preallocate(100);
        for(NODE_ITERATOR<TV> iterator(mpi_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next())if(particles(iterator.Node_Index())){
            T_PARTICLES* cell_particles=particles(iterator.Node_Index());
            while(cell_particles){for(int i=0;i<cell_particles->Size();i++) if(!block_domain.Thickened(3*mpi_grid.local_grid.dX.Min()).Lazy_Inside(cell_particles->X(i))) 
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(cell_particles,i));cell_particles=cell_particles->next;}} // TODO: delete the particle locally?
        requests.Append(ISend_Particles(mpi_grid,exchange_particles(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}
    // probe and receive
    for(int message=0;message<requests.m;message++){
        MPI::Status probe_status;
        mpi_grid.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Block_Particles(mpi_grid,template_particles,particles,tag,probe_status,particle_levelset);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function Exchange_Ghost_Particles
//#####################################################################
template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Exchange_Ghost_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset)
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    STATIC_ASSERT((is_same<T_PARTICLES,typename remove_pointer<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE>::value));
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<RANGE<TV_INT> > send_regions;
    // this way Find_Boundary_Regions will return block indices which for uniform grids are the node index not minimum corner cell
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector());
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);
    // send particles that are in the ghost cells of an adjacent processor
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(GRID<TV>::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and send a lot more corner particles than it should
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int n=0;n<send_regions.m;n++) if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        exchange_particles(n).Preallocate(100);
        for(NODE_ITERATOR<TV> iterator(mpi_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next())if(particles(iterator.Node_Index())){
            T_PARTICLES* cell_particles=particles(iterator.Node_Index());
            while(cell_particles){for(int i=0;i<cell_particles->Size();i++) if(domain.Lazy_Inside(cell_particles->X(i))) // could be optimized further
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(cell_particles,i));cell_particles=cell_particles->next;}} // TODO: delete the particle locally?
        requests.Append(ISend_Particles(mpi_grid,exchange_particles(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}
    // probe and receive
    for(int message=0;message<requests.m;message++){
        MPI::Status probe_status;
        mpi_grid.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Ghost_Particles(mpi_grid,template_particles,particles,tag,probe_status,bandwidth,particle_levelset);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function Exchange_Boundary_Particles
//#####################################################################
template<class TV,class T_PARTICLES> void
Exchange_Boundary_Particles_Flat(const MPI_UNIFORM_GRID<TV>& mpi_grid,T_PARTICLES& particles,const typename TV::SCALAR ghost_distance)
{
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers;
    ISend_Particles(mpi_grid,particles,ghost_distance,tag,buffers,requests);
    for(int message=0;message<requests.m;message++)Recv_Particles(mpi_grid,particles,tag);
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################

#else

//#####################################################################
template<class TV,class T_PARTICLES> void Exchange_Boundary_Particles_Flat(const MPI_UNIFORM_GRID<TV>&,T_PARTICLES&,const typename TV::SCALAR)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES> void Exchange_Boundary_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,
    const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES> void Exchange_Overlapping_Block_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,
    const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES> void Exchange_Ghost_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,
    T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset)
{PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################

#endif

//#####################################################################
template void Exchange_Boundary_Particles<VECTOR<float,1>,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> >*,VECTOR<int,1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> >*,VECTOR<int,1> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,1> >&);
template void Exchange_Boundary_Particles<VECTOR<float,1>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> >*,VECTOR<int,1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> >*,VECTOR<int,1> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,1> >&);
template void Exchange_Boundary_Particles<VECTOR<float,2>,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >*,VECTOR<int,2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >*,VECTOR<int,2> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,2> >&);
template void Exchange_Boundary_Particles<VECTOR<float,2>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >*,VECTOR<int,2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >*,VECTOR<int,2> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,2> >&);
template void Exchange_Boundary_Particles<VECTOR<float,3>,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >*,VECTOR<int,3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >*,VECTOR<int,3> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,3> >&);
template void Exchange_Boundary_Particles<VECTOR<float,3>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> >*,VECTOR<int,3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> >*,VECTOR<int,3> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,3> >&);
template void Exchange_Boundary_Particles_Flat<VECTOR<double,3>,VORTICITY_PARTICLES<VECTOR<double,3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,VORTICITY_PARTICLES<VECTOR<double,3> >&,GRID<VECTOR<double,3> >::SCALAR);
template void Exchange_Boundary_Particles_Flat<VECTOR<float,3>,VORTICITY_PARTICLES<VECTOR<float,3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,VORTICITY_PARTICLES<VECTOR<float,3> >&,GRID<VECTOR<float,3> >::SCALAR);
template void Exchange_Ghost_Particles<VECTOR<double,1>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> >*,VECTOR<int,1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> >*,VECTOR<int,1> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,1> >&);
template void Exchange_Ghost_Particles<VECTOR<double,2>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> >*,VECTOR<int,2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> >*,VECTOR<int,2> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,2> >&);
template void Exchange_Ghost_Particles<VECTOR<double,3>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> >*,VECTOR<int,3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> >*,VECTOR<int,3> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,3> >&);
template void Exchange_Overlapping_Block_Particles<VECTOR<double,1>,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> >*,VECTOR<int,1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> >*,VECTOR<int,1> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,1> >&);
template void Exchange_Overlapping_Block_Particles<VECTOR<double,2>,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >*,VECTOR<int,2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >*,VECTOR<int,2> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,2> >&);
template void Exchange_Overlapping_Block_Particles<VECTOR<double,3>,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >*,VECTOR<int,3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >*,VECTOR<int,3> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,3> >&);
template void Exchange_Boundary_Particles<VECTOR<double,1>,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> >*,VECTOR<int,1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> >*,VECTOR<int,1> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,1> >&);
template void Exchange_Boundary_Particles<VECTOR<double,1>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> >*,VECTOR<int,1> > >(
    MPI_UNIFORM_GRID<VECTOR<double,1> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> >*,VECTOR<int,1> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,1> >&);
template void Exchange_Boundary_Particles<VECTOR<double,2>,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >*,VECTOR<int,2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >*,VECTOR<int,2> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,2> >&);
template void Exchange_Boundary_Particles<VECTOR<double,2>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> >*,VECTOR<int,2> > >(
    MPI_UNIFORM_GRID<VECTOR<double,2> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> >*,VECTOR<int,2> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,2> >&);
template void Exchange_Boundary_Particles<VECTOR<double,3>,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >*,VECTOR<int,3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >*,VECTOR<int,3> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,3> >&);
template void Exchange_Boundary_Particles<VECTOR<double,3>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> >*,VECTOR<int,3> > >(
    MPI_UNIFORM_GRID<VECTOR<double,3> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> >*,VECTOR<int,3> >&,int,
    PARTICLE_LEVELSET<VECTOR<double,3> >&);
template void Exchange_Ghost_Particles<VECTOR<float,1>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> >*,VECTOR<int,1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> >*,VECTOR<int,1> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,1> >&);
template void Exchange_Ghost_Particles<VECTOR<float,2>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >*,VECTOR<int,2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >*,VECTOR<int,2> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,2> >&);
template void Exchange_Ghost_Particles<VECTOR<float,3>,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> >,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> >*,VECTOR<int,3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> > const&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> >*,VECTOR<int,3> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,3> >&);
template void Exchange_Overlapping_Block_Particles<VECTOR<float,1>,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> >*,VECTOR<int,1> > >(
    MPI_UNIFORM_GRID<VECTOR<float,1> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> >*,VECTOR<int,1> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,1> >&);
template void Exchange_Overlapping_Block_Particles<VECTOR<float,2>,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >*,VECTOR<int,2> > >(
    MPI_UNIFORM_GRID<VECTOR<float,2> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >*,VECTOR<int,2> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,2> >&);
template void Exchange_Overlapping_Block_Particles<VECTOR<float,3>,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >*,VECTOR<int,3> > >(
    MPI_UNIFORM_GRID<VECTOR<float,3> > const&,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> > const&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >*,VECTOR<int,3> >&,int,
    PARTICLE_LEVELSET<VECTOR<float,3> >&);
}
