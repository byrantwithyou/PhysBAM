//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Projection/BOUNDARY_CONDITION_DOUBLE_FINE.h>
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class TV> BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
BOUNDARY_CONDITION_DOUBLE_FINE(const GRID<TV>& mac_grid,int ghost)
    :mac_grid(mac_grid),grid(mac_grid.numbers_of_cells*2+1,mac_grid.domain,false),
    ghost(ghost),bc_type(grid.Node_Indices(2*ghost+1))
{
}
//#####################################################################
// Function Reset
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Reset(char type)
{
    bc_type.Fill(type);
    bc_p.Remove_All();;
    bc_u.Remove_All();;
}
//#####################################################################
// Function Set
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Set(const LEVELSET<TV>& ls,char type,std::function<T(const TV& X)> f,bool thin,bool invert,T contour)
{
    PHYSBAM_ASSERT(!thin || type==bc_noslip || type==bc_slip);
    auto Phi=[&ls,invert,contour,this](const TV_INT& i)
        {
            T phi=ls.Extended_Phi(grid.Node(i))-contour;
            if(invert) phi=-phi;
            return phi;
        };
    if(thin)
    {
        auto Do_Face=[Phi,this,f](const FACE_INDEX<TV::m>& face,int a)
            {
                TV_INT i=(face.index*2+1).Add_Axis(a,-1);
                TV_INT i0=i.Add_Axis(a,-1),i1=i.Add_Axis(a,1);
                T phi0=Phi(i0),phi1=Phi(i1);
                if((phi0>0)!=(phi1>0))
                {
                    TV X0=grid.Node(i0),X1=grid.Node(i1);
                    T t=phi0/(phi0-phi1);
                    TV Z=X0+(X1-X0)*t;
                    auto& z=bc_u.Get_Or_Insert(face);
                    z.x+=f(Z);
                    z.y++;
                    z.z=true;
                }
            };
        for(FACE_ITERATOR<TV> it(mac_grid,ghost);it.Valid();it.Next())
        {
            if(type==bc_slip)
                Do_Face(it.Full_Index(),it.Axis());
            else
                for(int a=0;a<TV::m;a++)
                    Do_Face(it.Full_Index(),a);
        }
    }
    else
    {
        for(RANGE_ITERATOR<TV::m> it(bc_type.domain);it.Valid();it.Next())
        {
            if((it.index.Sum()+TV::m)%2) continue; // don't compute if we don't need
            if(Phi(it.index)>0) continue;
            char& c=bc_type(it.index);
            if(c>=bc_free || type==bc_noslip) c=type;
        }
    }
}
//#####################################################################
// Function Set
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Set(const T_SURFACE& surface,char type,std::function<T(const TV& X, int e)> f,bool thin)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Get_Pressure_Boundary_Conditions
//#####################################################################
template<class TV> void BOUNDARY_CONDITION_DOUBLE_FINE<TV>::
Get_Pressure_Boundary_Conditions(ARRAY<bool,TV_INT>& psi_D,
    ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& p) const
{
    for(RANGE_ITERATOR<TV::m> it(psi_D.domain);it.Valid();it.Next())
    {
        TV_INT bc_i=it.index*2+1;
        char bc_t=bc_type(bc_i);
        if(bc_t<0)
        {
            psi_D(it.index)=true;
            continue;
        }
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
                }
                else psi_N(a,it.index.Add_Axis(a,s))=true;
            }
        }
    }
    for(auto z:bc_u)
        if(z.data.z)
            psi_N(z.key)=true;
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
