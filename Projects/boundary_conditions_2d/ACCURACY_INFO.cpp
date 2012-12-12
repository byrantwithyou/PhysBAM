#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include "ACCURACY_INFO.h"
using namespace PhysBAM;

static void Compute_Helper(ACCURACY_INFO<1>& ai)
{
    PHYSBAM_ASSERT(ai.resolution%ai.base_resolution==0 && ai.resolution/ai.base_resolution%2==1);
    int m=ai.resolution/ai.base_resolution,n=m/2;
    VECTOR<int,1> i;
    for(i.x=ai.sample_domain.min_corner.x+m;i.x<=ai.sample_domain.max_corner.x;i.x+=m)
        ai.face_samples.Append(FACE_INDEX<1>(1,i));
    for(i.x=ai.sample_domain.min_corner.x+n;i.x<=ai.sample_domain.max_corner.x;i.x+=m)
        ai.cell_samples.Append(i);
}

static void Compute_Helper(ACCURACY_INFO<2>& ai)
{
    PHYSBAM_ASSERT(ai.resolution%ai.base_resolution==0 && ai.resolution/ai.base_resolution%2==1);
    int m=ai.resolution/ai.base_resolution,n=m/2;
    VECTOR<int,2> i;
    for(i.y=ai.sample_domain.min_corner.y+n;i.y<=ai.sample_domain.max_corner.y;i.y+=m)
        for(i.x=ai.sample_domain.min_corner.x+m;i.x<=ai.sample_domain.max_corner.x;i.x+=m)
            ai.face_samples.Append(FACE_INDEX<2>(1,i));
    for(i.x=ai.sample_domain.min_corner.x+n;i.x<=ai.sample_domain.max_corner.x;i.x+=m)
        for(i.y=ai.sample_domain.min_corner.y+m;i.y<=ai.sample_domain.max_corner.y;i.y+=m)
            ai.face_samples.Append(FACE_INDEX<2>(2,i));
    for(i.y=ai.sample_domain.min_corner.y+n;i.y<=ai.sample_domain.max_corner.y;i.y+=m)
        for(i.x=ai.sample_domain.min_corner.x+n;i.x<=ai.sample_domain.max_corner.x;i.x+=m)
            ai.cell_samples.Append(i);
}

template<int d> void ACCURACY_INFO<d>::
Compute()
{
    Compute_Helper(*this);
}

template<int d> template<class T> void ACCURACY_INFO<d>::
Print(const char* label,ARRAY<T,TV_INT>& a) const
{
    if(cell_samples.m)
        LOG::cout<<"W@Z-"<<label<<" "<<resolution<<" "<<a.Subset(cell_samples)<<std::endl;
}

template<int d> template<class T> void ACCURACY_INFO<d>::
Print(const char* label,ARRAY<T,FACE_INDEX<d> >& a) const
{
    if(face_samples.m)
        LOG::cout<<"W@Z-"<<label<<" "<<resolution<<" "<<a.Subset(face_samples)<<std::endl;
}

template<int d> template<class TV> void ACCURACY_INFO<d>::
Print_Locations(const GRID<TV>& grid) const
{
    LOG::cout<<"CELLS-I "<<cell_samples<<std::endl;
    LOG::cout<<"FACES-I "<<face_samples<<std::endl;
    LOG::cout<<"CELLS-X ";
    for(int i=0;i<cell_samples.m;i++) LOG::cout<<grid.X(cell_samples(i))<<" ";
    LOG::cout<<std::endl;
    LOG::cout<<"FACES-X ";
    for(int i=0;i<face_samples.m;i++) LOG::cout<<grid.Face(face_samples(i))<<" ";
    LOG::cout<<std::endl;
}

template struct ACCURACY_INFO<1>;
template struct ACCURACY_INFO<2>;
template void ACCURACY_INFO<1>::Print<double>(char const*,ARRAY<double,FACE_INDEX<1> >&) const;
template void ACCURACY_INFO<1>::Print<double>(char const*,ARRAY<double,VECTOR<int,1> >&) const;
template void ACCURACY_INFO<1>::Print_Locations<VECTOR<double,1> >(GRID<VECTOR<double,1> > const&) const;
template void ACCURACY_INFO<1>::Print_Locations<VECTOR<float,1> >(GRID<VECTOR<float,1> > const&) const;
template void ACCURACY_INFO<2>::Print<double>(char const*,ARRAY<double,FACE_INDEX<2> >&) const;
template void ACCURACY_INFO<2>::Print<double>(char const*,ARRAY<double,VECTOR<int,2> >&) const;
template void ACCURACY_INFO<2>::Print_Locations<VECTOR<double,2> >(GRID<VECTOR<double,2> > const&) const;
template void ACCURACY_INFO<2>::Print_Locations<VECTOR<float,2> >(GRID<VECTOR<float,2> > const&) const;
