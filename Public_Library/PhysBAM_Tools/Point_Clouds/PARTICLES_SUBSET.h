//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLES_SUBSET
//#####################################################################
#ifndef __PARTICLES_SUBSET__
#define __PARTICLES_SUBSET__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
namespace PhysBAM{

template<class TV,class T_PARTICLES> //At the moment this only works for PARTICLE
class PARTICLES_SUBSET
{
    typedef typename TV::SCALAR T;
public:
    typedef int HAS_UNTYPED_READ_WRITE;

    T_PARTICLES& point_cloud;
    ARRAY<int> active_indices; // TODO: improve names
    ARRAY<int> subset_index_from_point_cloud_index;
    ARRAY_COLLECTION* array_collection; //hack to be compatable with point_clouds

    PARTICLES_SUBSET(T_PARTICLES& point_cloud_input)
        :point_cloud(point_cloud_input),array_collection(new ARRAY_COLLECTION)
    {array_collection->number=active_indices.m;}

    ~PARTICLES_SUBSET()
    {delete array_collection;}

    void Clean_Memory()
    {array_collection->Clean_Memory();active_indices.Clean_Memory();subset_index_from_point_cloud_index.Clean_Memory();}

    INDIRECT_ARRAY<ARRAY_VIEW<TV> > X()
    {return point_cloud.X.Subset(active_indices);}

    INDIRECT_ARRAY<ARRAY_VIEW<TV> > X() const
    {return point_cloud.X.Subset(active_indices);}

    TV& X(const int i)
    {return point_cloud.X(active_indices(i));}

    const TV& X(const int i) const
    {return point_cloud.X(active_indices(i));}

    INDIRECT_ARRAY<ARRAY_VIEW<TV> > V()
    {return point_cloud.V.Subset(active_indices);}

    TV& V(const int i)
    {return point_cloud.V(active_indices(i));}

    const TV& V(const int i) const
    {return point_cloud.V(active_indices(i));}

    INDIRECT_ARRAY<ARRAY_VIEW<T> > mass()
    {return point_cloud.mass.Subset(active_indices);}

    T& mass(const int i)
    {return point_cloud.mass(active_indices(i));}

    T mass(const int i) const
    {return point_cloud.mass(active_indices(i));}

    int Add_Element()
    {assert(subset_index_from_point_cloud_index.m==point_cloud.array_collection->Size());
    int id=active_indices.Append(point_cloud.array_collection->Add_Element());subset_index_from_point_cloud_index.Append(id);
    array_collection->number=active_indices.m;
    return id;}

    void Add_Elements(const int new_point_cloud)
    {for(int p=0;p<new_point_cloud;p++) Add_Element();}

    int Add_Existing_Element_If_Not_Already_There(const int p)
    {assert(subset_index_from_point_cloud_index.m==point_cloud.number);
    if(subset_index_from_point_cloud_index(p)<0) subset_index_from_point_cloud_index(p)=active_indices.Append(p);
    array_collection->number=active_indices.m;
    return subset_index_from_point_cloud_index(p);}

    void Copy_Element(const int from,const int to)
    {point_cloud.Copy_Element(active_indices(from),active_indices(to));}

    void Update_Subset_Index_From_Element_Index()
    {subset_index_from_point_cloud_index.Resize(point_cloud.array_collection->Size(),false,false);subset_index_from_point_cloud_index.Fill(0);
    for(int p=0;p<active_indices.m;p++)subset_index_from_point_cloud_index(active_indices(p))=p;}

    void Update_Number_Nodes()
    {subset_index_from_point_cloud_index.Resize(point_cloud.array_collection->Size());}

    void Initialize_Subset(const PARTICLES_SUBSET& subset)
    {active_indices=subset.active_indices;subset_index_from_point_cloud_index=subset.subset_index_from_point_cloud_index;array_collection->number=active_indices.m;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,active_indices);Update_Subset_Index_From_Element_Index();}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,active_indices);}
//#####################################################################
};
}
#endif
