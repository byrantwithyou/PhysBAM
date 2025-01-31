//#####################################################################
// Copyright 2005, Ron Fedkiw, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_MULTIPLE
//##################################################################### 
#ifndef __LEVELSET_MULTIPLE__
#define __LEVELSET_MULTIPLE__ 

#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
namespace PhysBAM{

template<class TV,class T2> class BOUNDARY;

template<class TV> class LEVELSET_CALLBACKS; // TODO: invalid dependency
template<class TV>
class LEVELSET_MULTIPLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef INTERPOLATION_UNIFORM<TV,T> T_INTERPOLATION_SCALAR;
    typedef INTERPOLATION_UNIFORM<TV,TV> T_INTERPOLATION_VECTOR;
public:
    typedef TV VECTOR_T;

    GRID<TV>& grid;
    ARRAY<ARRAY<T,TV_INT>>& phis;
    LEVELSET_CALLBACKS<TV>* levelset_callbacks;
    ARRAY<LEVELSET<TV>*> levelsets;
    bool use_external_levelsets;

    LEVELSET_MULTIPLE(GRID<TV>& grid_input,ARRAY<ARRAY<T,TV_INT>>& phis_input,const bool use_external_levelsets_input=false);
    LEVELSET_MULTIPLE(const LEVELSET_MULTIPLE&) = delete;
    void operator=(const LEVELSET_MULTIPLE&) = delete;
    ~LEVELSET_MULTIPLE();

    T Phi(const int region,const TV_INT& index) const
    {return phis(region)(index);}

    T Phi(const int region,const TV& location) const
    {return levelsets(region)->Phi(location);}

    T Region_Value(const TV_INT& index,const ARRAY<T>& value) const
    {return value(Inside_Region(index));}
    
    T Region_Value(const TV& location,const ARRAY<T>& value) const
    {return value(Inside_Region(location));}

    bool Interface(const TV_INT& index_1,const TV_INT& index_2) const
    {return Inside_Region(index_1)!=Inside_Region(index_2);}    
    
    bool Interface(const TV& location_1,const TV& location_2) const
    {return Inside_Region(location_1)!=Inside_Region(location_2);}
    
    static T Sign(const int region,const int other_region)
    {assert(region!=other_region);
    if(region < other_region) return -1;else return 1;}

    static void Project_Levelset(T& minimum_phi,T& second_minimum_phi)
    {T correction=(T).5*(minimum_phi+second_minimum_phi);
    minimum_phi-=correction;second_minimum_phi-=correction;}

    void Minimum_Regions(const TV_INT& index_1,const TV_INT& index_2,int& region_1,int& region_2,T& phi_1,T& phi_2) const
    {region_1=Inside_Region(index_1,phi_1);region_2=Inside_Region(index_2,phi_2);}

    void Minimum_Regions_Projected(const TV& location_1,const TV& location_2,int& region_1,int& region_2,T& phi_1,T& phi_2) const
    {int region_temp;T phi_temp;
    Two_Minimum_Regions(location_1,region_1,region_temp,phi_1,phi_temp);Project_Levelset(phi_1,phi_temp);
    Two_Minimum_Regions(location_2,region_2,region_temp,phi_2,phi_temp);Project_Levelset(phi_2,phi_temp);}

    T Heaviside(const TV_INT& index,const ARRAY<T>& value,const T half_width=0) const
    {int region_1,region_2;T phi_1,phi_2;Two_Minimum_Regions(index,region_1,region_2,phi_1,phi_2);
    return LEVELSET_UTILITIES<T>::Heaviside(phi_1,value(region_1),value(region_2),half_width);}

    T Heaviside(const TV& location,const ARRAY<T>& value,const T half_width=0) const
    {int region_1,region_2;T phi_1,phi_2;Two_Minimum_Regions(location,region_1,region_2,phi_1,phi_2);Project_Levelset(phi_1,phi_2);
    return LEVELSET_UTILITIES<T>::Heaviside(phi_1,value(region_1),value(region_2),half_width);}

    T Heaviside(const TV_INT& index_1,const TV_INT& index_2,const ARRAY<T>& value,const T half_width=0) const
    {int region_1,region_2;T phi_1,phi_2;Minimum_Regions(index_1,index_2,region_1,region_2,phi_1,phi_2);
    return LEVELSET_UTILITIES<T>::Heaviside((T).5*(phi_1-phi_2),value(region_1),value(region_2),half_width);}

    T Heaviside(const TV& location_1,const TV& location_2,const ARRAY<T>& value,const T half_width=0) const
    {int region_1,region_2;T phi_1,phi_2;Minimum_Regions_Projected(location_1,location_2,region_1,region_2,phi_1,phi_2);
    return LEVELSET_UTILITIES<T>::Heaviside((T).5*(phi_1-phi_2),value(region_1),value(region_2),half_width);}

    T Convex_Average(const TV_INT& index_1,const TV_INT& index_2,const ARRAY<T>& value) const
    {int region_1,region_2;T phi_1,phi_2;Minimum_Regions(index_1,index_2,region_1,region_2,phi_1,phi_2);
    return LEVELSET_UTILITIES<T>::Average(phi_1,value(region_1),-phi_2,value(region_2));}

    T Convex_Average(const TV& location_1,const TV& location_2,const ARRAY<T>& value) const
    {int region_1,region_2;T phi_1,phi_2;Minimum_Regions_Projected(location_1,location_2,region_1,region_2,phi_1,phi_2);
    return LEVELSET_UTILITIES<T>::Average(phi_1,value(region_1),-phi_2,value(region_2));}

//#####################################################################
    void Recreate_Levelsets();
    void Fill_Ghost_Cells(ARRAY<ARRAY<T,TV_INT> >& phi_ghost,const T time,const int number_of_ghost_cells);
    void Set_Custom_Boundary(BOUNDARY<TV,T>& boundary_input);
    void Set_Custom_Secondary_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input);
    void Set_Custom_Normal_Interpolation(T_INTERPOLATION_VECTOR& normal_interpolation_input);
    void Set_Custom_Curvature_Interpolation(T_INTERPOLATION_SCALAR& curvature_interpolation_input);
    void Set_Levelset_Callbacks(LEVELSET_CALLBACKS<TV>& levelset_callbacks_input);
    int Inside_Region(const TV_INT& index) const; // assumes exactly one Phi<0 on a node
    int Inside_Region(const TV_INT& index,T& phi) const; // assumes exactly one Phi<0 on a node
    int Inside_Region(const TV& location) const;
    int Inside_Region(const TV& location,T& phi) const;
    int Inside_Region_Face(const int axis,const TV_INT& face_index) const; // does not assume exactly one Phi<0
    void Two_Minimum_Regions(const TV_INT& index,int& minimum_region,int& second_minimum_region,T& minimum_phi,T& second_minimum_phi) const;
    void Two_Minimum_Regions(const TV& location,int& minimum_region,int& second_minimum_region,T& minimum_phi,T& second_minimum_phi) const;
    T CFL(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities) const;
    T CFL(const ARRAY<TV,TV_INT>& velocity) const;
    bool Is_Projected_At_Nodes();
    void Compute_Normals(const T time=0);
    void Compute_Curvature(const T time=0);
    void Fast_Marching_Method(const ARRAY<int>& local_advection_spatial_orders,int process_sign=0);
    void Project_Levelset(const int number_of_ghost_cells=0);
    void Get_Single_Levelset(const ARRAY<bool>& positive_regions,LEVELSET<TV>& levelset,const bool flood_fill_for_bubbles);
//#####################################################################
};   
}
#endif

