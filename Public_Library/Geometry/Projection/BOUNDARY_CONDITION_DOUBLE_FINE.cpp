//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Math_Tools/int_div.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Intersections/GRID_SURFACE_INTERSECTION.h>
#include <Geometry/Projection/BOUNDARY_CONDITION_DOUBLE_FINE.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class TV> BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
BOUNDARY_CONDITION_DOUBLE_FINE(const GRID<TV>& mac_grid,int ghost)
    :mac_grid(mac_grid),grid(mac_grid.numbers_of_cells*2+1,mac_grid.domain,false),
    ghost(ghost),bc_type(grid.Node_Indices(2*ghost+1)),bc_type_current(grid.Node_Indices(2*ghost+1))
{
}
//#####################################################################
// Function Reset
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Reset(char type)
{
    bc_type.Fill(type);
    bc_type_current.Fill(type);
    bc_p.Remove_All();
    bc_u.Remove_All();
}
inline void set_bc(auto& bc_type,auto& bc_type_current,const auto& i,char type)
{
    enum bc_enum {bc_slip=-3,bc_noslip=-2,bc_free=-1};
    bc_type_current(i)=type;
    char& c=bc_type(i);
    if(c>=bc_free || type==bc_noslip) c=type;
}
//#####################################################################
// Function Set
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Set(const IMPLICIT_OBJECT<TV>* io,char type,std::function<T(const TV& X,int a)> f,bool thin,bool invert,T contour)
{
    PHYSBAM_ASSERT(!thin || type==bc_noslip || type==bc_slip);
    RANGE<TV_INT> domain=bc_type.domain;
    for(RANGE_ITERATOR<TV::m> it(domain);it.Valid();it.Next())
    {
        if((it.index.Sum()+TV::m)%2) continue; // don't compute if we don't need
        T phi=io->Extended_Phi(grid.Node(it.index))-contour;
        if(invert) phi=-phi;
        if(phi>0) bc_type_current(it.index)=0;
        else set_bc(bc_type,bc_type_current,it.index,type);
    }
    if(f)
        Set_Current_Values(domain,type,
            [this,io,contour,&f](const TV_INT& i0,const TV_INT& i1,int a)
            {
                TV X0=grid.Node(i0),X1=grid.Node(i1);
                T phi0=io->Extended_Phi(grid.Node(i0))-contour;
                T phi1=io->Extended_Phi(grid.Node(i1))-contour;
                T t=phi0/(phi0-phi1);
                TV Z=X0+(X1-X0)*t;
                return f(Z,a);
            },thin);
    else Set_Current_Values(domain,type,0,thin);
}
//#####################################################################
// Function Set
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Set(const T_SURFACE& surface,char type,std::function<T(const TV& X,int a)> f,bool thin)
{
    PHYSBAM_ASSERT(!thin || type==bc_noslip || type==bc_slip);

    ARRAY<bool,TV_INT> flood_fill_array(bc_type.domain.Edge_Lengths()/2+1);

    for(RANGE_ITERATOR<TV::m> it(TV_INT()+2);it.Valid();it.Next())
    {
        if((it.index.Sum()+TV::m)%2) continue; // don't compute if we don't need
        if(type!=bc_noslip && it.index!=TV_INT()+1) continue;

        TV X0=grid.Node(bc_type.domain.min_corner+it.index);
        TV X1=grid.Node(bc_type.domain.max_corner-it.index);
        TV_INT counts=bc_type.domain.Edge_Lengths()/2-it.index;
        TV_INT offset=bc_type.domain.min_corner+it.index;
        GRID<TV> vis_grid(counts,RANGE<TV>(X0,X1),false);
        // fine = 2 * vis + offset;

        HASHTABLE<EDGE_INDEX<TV::m>,GRID_SURFACE_INTERSECTION_DATA<TV> > hash;
        Grid_Surface_Intersection(hash,vis_grid,surface,!thin);

        if(!thin)
        {
            Flood_Fill(flood_fill_array,hash);
            for(RANGE_ITERATOR<TV::m> it2(counts+1);it2.Valid();it2.Next())
                if(flood_fill_array(it2.index))
                    set_bc(bc_type,bc_type_current,2*it2.index+offset,type);
        }

        for(const auto& h:hash)
        {
            TV_INT i[2]={2*h.key.index+offset};
            i[1]=i[0].Add_Axis(h.key.axis,2);
            TV X0=grid.Node(i[0]);
            TV X1=grid.Node(i[1]);
            T y=0;
            for(const auto& t:h.data.cut_elements)
            {
                TV Z=X0+(X1-X0)*t.y;
                y+=f(Z,0xDEADBEEF);
            }
            y/=h.data.cut_elements.m;
            if(type==bc_free)
            {
                for(int j=0;j<2;j++)
                    if(!h.data.in[j])
                    {
                        auto& z=bc_p.Get_Or_Insert(i[j]);
                        if(f) z.x+=y;
                        z.y++;
                    }
            }
            else
            {
                FACE_INDEX<TV::m> face;
                // TODO
                auto& z=bc_u.Get_Or_Insert(face);
                if(f) z.x+=y;
                z.y++;
                if(thin) z.z=true;
            }
        }
    }
}
//#####################################################################
// Function Set
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Set(T_OBJECT& object,char type,std::function<T(const TV& X,int a)> f,bool thin)
{
    return Set(object.Get_Boundary_Object(),type,f,thin);
}
//#####################################################################
// Function Set_Domain_Walls
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Set_Domain_Walls(int side_mask,char type,std::function<T(const TV& X,int a)> f)
{
    for(RANGE_ITERATOR<TV::m> it(bc_type.domain,grid.Node_Indices(),RI::side_mask|RI::ghost,side_mask);it.Valid();it.Next())
    {
        if((it.index.Sum()+TV::m)%2) continue; // don't compute if we don't need
        set_bc(bc_type,bc_type_current,it.index,type);
    }
    if(type==bc_free)
    {
        for(FACE_RANGE_ITERATOR<TV::m> it(mac_grid.Domain_Indices(),0,0,RF::side_mask|RF::ghost,side_mask);it.Valid();it.Next())
        {
            auto& z=bc_p.Get_Or_Insert(it.face.Cell_Index(it.side%2));
            if(f) z.x+=f(mac_grid.Face(it.face),0);
            z.y++;
        }
    }
    else
    {
        bool ns=type==bc_noslip;
        auto flags=RF::side_mask|RF::ghost;
        if(ns) flags|=RF::skip_outer;
        for(FACE_RANGE_ITERATOR<TV::m> it(mac_grid.Domain_Indices(),ns,0,flags,side_mask);it.Valid();it.Next())
        {
            auto& z=bc_u.Get_Or_Insert(it.face);
            if(f) z.x+=f(mac_grid.Face(it.face),it.face.axis);
            z.y++;
        }
    }
}
//#####################################################################
// Function Set_Current_Values
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Set_Current_Values(const RANGE<TV_INT>& domain,char type,std::function<T(const TV_INT& i0,const TV_INT& i1,int a)> f,bool thin)
{
    for(int axis=0;axis<TV::m;axis++)
    {
        int mn=bc_noslip?0:axis;
        int mx=bc_noslip?TV::m-1:axis;
        for(int a=mn;a<=mx;a++)
        {
            TV_INT off=TV_INT().Add_Axis(axis,1).Add_Axis(a,1);
            RANGE<TV_INT> range(cdiv(bc_type.domain.min_corner-1+off,2),fdiv(bc_type.domain.max_corner-off,2));
            for(FACE_RANGE_ITERATOR<TV::m> it(range,RF::none,axis);it.Valid();it.Next())
            {
                TV_INT i=(it.face.index*2+1).Add_Axis(it.face.axis,-1);
                TV_INT i0=i.Add_Axis(a,-1),i1=i.Add_Axis(a,1);
                char t0=bc_type_current(i0),t1=bc_type_current(i1);
                if(t0>=0 && t1>=0) continue;
                if(t0==type && t1==type) continue;
                if(type==bc_free)
                {
                    auto& z=bc_p.Get_Or_Insert(it.face.Cell_Index(t0>=0));
                    if(f) z.x+=f(i0,i1,it.face.axis);
                    z.y++;
                }
                else
                {
                    auto& z=bc_u.Get_Or_Insert(it.face);
                    if(f) z.x+=f(i0,i1,0);
                    z.y++;
                    z.z=thin;
                }
            }
        }
    }
}
//#####################################################################
// Function Get_Pressure_Boundary_Conditions
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Get_Pressure_Boundary_Conditions(ARRAY<bool,TV_INT>& psi_D,
    ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& p,
    ARRAY<T,FACE_INDEX<TV::m> >& u) const
{
    psi_N.Fill(false);
    for(RANGE_ITERATOR<TV::m> it(psi_D.domain);it.Valid();it.Next())
    {
        TV_INT bc_i=it.index*2+1;
        char bc_t=bc_type(bc_i);
        if(bc_t<0)
        {
            psi_D(it.index)=true;
            continue;
        }
        psi_D(it.index)=false;
        for(int a=0;a<TV::m;a++)
        {
            for(int s=0;s<2;s++)
            {
                int o=2*s-1;
                TV_INT bc_n=bc_i.Add_Axis(a,2*o);
                char bc_tn=bc_type(bc_n);
                if(bc_tn>=0) continue;
                if(bc_tn==bc_free)
                {
                    PAIR<T,int> z;
                    if(bc_p.Get(bc_n,z))
                        p(it.index.Add_Axis(a,o))=z.x/z.y;
                    else p(it.index.Add_Axis(a,o))=0;
                }
                else
                {
                    FACE_INDEX<TV::m> face(a,it.index.Add_Axis(a,s));
                    TRIPLE<T,int,bool> z;
                    if(bc_u.Get(face,z))
                        u(face)=z.x/z.y;
                    else u(face)=0;
                    psi_N(face)=true;
                }
            }
        }
    }
    for(auto z:bc_u)
        if(z.data.z)
        {
            psi_N(z.key)=true;
            u(z.key)=z.data.x/z.data.y;
        }
}
//#####################################################################
// Function Get_Viscosity_Boundary_Conditions
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Get_Viscosity_Boundary_Conditions(ARRAY<bool,TV_INT>& psi_D,
        ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& u,int axis) const
{
    TV_INT offset=(TV_INT()+1).Add_Axis(axis,-1);
    for(RANGE_ITERATOR<TV::m> it(psi_D.domain);it.Valid();it.Next())
    {
        TV_INT bc_i=it.index*2+offset;
        for(int a=0;a<TV::m;a++)
        {
            TV_INT i0=bc_i.Add_Axis(a,-1),i1=bc_i.Add_Axis(a,1);
            char bc0=bc_type(i0),bc1=bc_type(i1);
            if(bc0<0 && bc1<0)
            {
                psi_D(it.index)=true;
                continue;
            }
            int side=1;
            if(bc0<bc1){std::swap(bc0,bc1);side=0;}
            if(bc1>=0) continue;
            if(axis!=a && bc1==bc_noslip) bc1=bc_free;
            if(bc1==bc_free) psi_N(a,it.index.Add_Axis(a,side))=true;
            else
            {
                psi_D(it.index)=true;
                TRIPLE<T,int,bool> z;
                if(bc_u.Get({axis,it.index},z))
                    u(it.index)=z.x/z.y;
            }
        }
    }
    for(auto z:bc_u)
        if(z.data.z)
        {
            psi_D(z.key.index)=true;
            u(z.key.index)=z.data.x/z.data.y;
        }
}
template class BOUNDARY_CONDITION_DOUBLE_FINE<VECTOR<float,1> >;
template class BOUNDARY_CONDITION_DOUBLE_FINE<VECTOR<float,2> >;
template class BOUNDARY_CONDITION_DOUBLE_FINE<VECTOR<float,3> >;
template class BOUNDARY_CONDITION_DOUBLE_FINE<VECTOR<double,1> >;
template class BOUNDARY_CONDITION_DOUBLE_FINE<VECTOR<double,2> >;
template class BOUNDARY_CONDITION_DOUBLE_FINE<VECTOR<double,3> >;
}
