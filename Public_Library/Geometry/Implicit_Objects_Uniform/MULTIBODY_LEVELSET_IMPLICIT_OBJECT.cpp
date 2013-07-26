//#####################################################################
// Copyright 2002-2013, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Implicit_Objects_Uniform/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
MULTIBODY_LEVELSET_IMPLICIT_OBJECT(ARRAY<IMPLICIT_OBJECT<TV>*>* levelsets_input)
    :levelsets(levelsets_input),need_destroy_data(false)
{
    Update_Box();Update_Minimum_Cell_Size();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
~MULTIBODY_LEVELSET_IMPLICIT_OBJECT()
{
    if(need_destroy_data) Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV> MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Create()
{
    MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* object=new MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>(new ARRAY<IMPLICIT_OBJECT<TV>*>);
    object->need_destroy_data=true;
    return object;
}
//#####################################################################
// Function Delete_Pointers_And_Clean_Memory
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Delete_Pointers_And_Clean_Memory()
{
    for(int i=0;i<levelsets->m;i++) delete (*levelsets)(i);
    delete levelsets;
}
//#####################################################################
// Function Update_Box
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Update_Box()
{
    box=RANGE<TV>::Empty_Box();
    for(int i=0;i<levelsets->m;i++){
        (*levelsets)(i)->Update_Box();
        box.Enlarge_To_Include_Box((*levelsets)(i)->box);}
}
//#####################################################################
// Function Update_Minimum_Cell_Size
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Update_Minimum_Cell_Size(const int maximum_depth)
{
    minimum_cell_size=(T)FLT_MAX;
    for(int i=0;i<levelsets->m;i++){
        (*levelsets)(i)->Update_Minimum_Cell_Size(maximum_depth);
        minimum_cell_size=min((*levelsets)(i)->Minimum_Cell_Size(),minimum_cell_size);}
}
//#####################################################################
// Function Minimum_Cell_Size_Within_Box
//#####################################################################
template<class TV> typename TV::SCALAR MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const 
{
    T minimum_edge_length=(T)FLT_MAX;
    for(int i=0;i<levelsets->m;i++) minimum_edge_length=min((*levelsets)(i)->Minimum_Cell_Size_Within_Box(box),minimum_edge_length);
    return minimum_edge_length; // TODO: make this check overlap with the grids.
}
//#####################################################################
// operator()
//#####################################################################
template<class TV> typename TV::SCALAR MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
operator()(const TV& location) const
{
    T phi_value=(T)FLT_MAX;
    for(int i=0;i<levelsets->m;i++){
        if((*levelsets)(i)->box.Lazy_Outside(location)) phi_value=min((*levelsets)(i)->Extended_Phi(location),phi_value);
        else phi_value=min((*(*levelsets)(i))(location),phi_value);}
    return phi_value;
}
//#####################################################################
// Function Phi_With_Index
//#####################################################################
template<class TV> typename TV::SCALAR MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Phi_With_Index(const TV& location,int& levelset_index) const
{
    T phi_value=(T)FLT_MAX,new_phi_value;
    for(int i=0;i<levelsets->m;i++){
        if((*levelsets)(i)->box.Lazy_Outside(location)) new_phi_value=(*levelsets)(i)->Extended_Phi(location);
        else new_phi_value=(*(*levelsets)(i))(location);
        if(new_phi_value<phi_value){
            phi_value=new_phi_value;
            levelset_index=i;}}
    return phi_value;
}
//#####################################################################
// Function Extended_Phi
//#####################################################################
template<class TV> typename TV::SCALAR MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Extended_Phi(const TV& location) const
{
    T phi_value=(T)FLT_MAX;
    for(int i=0;i<levelsets->m;i++) phi_value=min((*levelsets)(i)->Extended_Phi(location),phi_value);
    return phi_value;
}
//#####################################################################
// Function Phi_Secondary
//#####################################################################
template<class TV> typename TV::SCALAR MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Phi_Secondary(const TV& location) const
{
    PHYSBAM_NOT_IMPLEMENTED(); // TODO: implement this, but need Phi_Secondary_Extended()!
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Normal(const TV& location,const int aggregate) const
{
    assert((aggregate >= 1 && aggregate <= 2*TV::m) || aggregate == -1);
    if(aggregate != -1) return box.Normal(aggregate);
    else{
        int index=-1;
        Phi_With_Index(location,index);
        return (*levelsets)(index)->Normal(location);}
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Extended_Normal(const TV& location,const int aggregate) const
{
    assert((aggregate >= 1 && aggregate <= 2*TV::m) || aggregate == -1);
    if(aggregate != -1) return box.Normal(aggregate);
    else{
        int index=-1;
        Phi_With_Index(location,index);return (*levelsets)(index)->Extended_Normal(location);}
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Compute_Normals() 
{
    for(int i=0;i<levelsets->m;i++) (*levelsets)(i)->Compute_Normals();
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    for(int i=0;i<levelsets->m;i++) (*levelsets)(i)->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
}
//#####################################################################
// Function Inflate
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Inflate(const T inflation_distance)
{
    for(int i=0;i<levelsets->m;i++) (*levelsets)(i)->Inflate(inflation_distance);
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class TV> typename TV::SCALAR MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Integration_Step(const T phi) const
{
    T step=FLT_MAX;
    for(int i=0;i<levelsets->m;i++) step=min(step,(*levelsets)(i)->Integration_Step(phi));
    return step;
}
//#####################################################################
// Function Minimum_Cell_Size
//#####################################################################
template<class TV> typename TV::SCALAR MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Minimum_Cell_Size() const 
{
    return minimum_cell_size;
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class TV> bool MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Inside(const TV& location,const T contour_value) const
{
    if(!box.Lazy_Inside(location)) return false;
    for(int i=0;i<levelsets->m;i++) if((*levelsets)(i)->Lazy_Inside(location,contour_value)) return true;
    return false;
}
//#####################################################################
// Function Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value) const
{
    if(!box.Lazy_Inside(location)) return false;
    bool levelset_value=false;
    phi_value=(T)FLT_MAX;
    T levelset_phi_value=T();
    for(int i=0;i<levelsets->m;i++){
        levelset_value=levelset_value||(*levelsets)(i)->Lazy_Inside_And_Value(location,levelset_phi_value,contour_value);
        if(levelset_phi_value<phi_value) phi_value=levelset_phi_value;}
    return levelset_value;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset
//#####################################################################
template<class TV> bool MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value) const
{
    for(int i=0;i<levelsets->m;i++) if((*levelsets)(i)->Lazy_Inside_Extended_Levelset(unclamped_X,contour_value)) return true;
    return false;
}
//#####################################################################
// Function Lazy_Inside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value) const
{
    bool levelset_value=false;
    phi_value=(T)FLT_MAX;
    T levelset_phi_value=T();
    for(int i=0;i<levelsets->m;i++){
        levelset_value=levelset_value||(*levelsets)(i)->Lazy_Inside_Extended_Levelset_And_Value(unclamped_X,levelset_phi_value,contour_value);
        if(levelset_phi_value<phi_value) phi_value=levelset_phi_value;}
    return levelset_value;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class TV> bool MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Outside(const TV& location,const T contour_value) const
{
    if(box.Lazy_Outside(location)) return true;
    bool levelset_value=true;
    for(int i=0;i<levelsets->m;i++) levelset_value=levelset_value && (*levelsets)(i)->Lazy_Outside(location,contour_value);
    return levelset_value;
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset
//#####################################################################
template<class TV> bool MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value) const
{
    bool levelset_value=true;
    for(int i=0;i<levelsets->m;i++) levelset_value=levelset_value && (*levelsets)(i)->Lazy_Outside_Extended_Levelset(unclamped_X,contour_value);
    return levelset_value;
}
//#####################################################################
// Function Lazy_Outside_Extended_Levelset_And_Value
//#####################################################################
template<class TV> bool MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value) const
{
    bool levelset_value=true;
    phi_value=(T)FLT_MAX;
    T levelset_phi_value=T();
    for(int i=0;i<levelsets->m;i++){
        levelset_value=levelset_value && (*levelsets)(i)->Lazy_Outside_Extended_Levelset_And_Value(unclamped_X,levelset_phi_value,contour_value);
        if(levelset_phi_value<phi_value) phi_value=levelset_phi_value;}
    return levelset_value;
}
//#####################################################################
// Function Min_Phi
//#####################################################################
template<class TV> typename TV::SCALAR MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Min_Phi() const
{
    T min_phi=(T)FLT_MAX;
    for(int i=0;i<levelsets->m;i++) min_phi=min(min_phi,(*levelsets)(i)->Min_Phi());
    return min_phi;
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Rescale(const T scaling_factor)
{
    for(int i=0;i<levelsets->m;i++) (*levelsets)(i)->Rescale(scaling_factor);
    Update_Box();Update_Minimum_Cell_Size();
}
//#####################################################################
// Function Translate
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Translate(const TV& translation)
{
    for(int i=0;i<levelsets->m;i++) (*levelsets)(i)->Translate(translation);
    Update_Box();Update_Minimum_Cell_Size();
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class TV> VECTOR<typename TV::SCALAR,TV::dimension-1> MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Principal_Curvatures(const TV& X) const
{
    int index=-1;
    Phi_With_Index(X,index);
    return (*levelsets)(index)->Principal_Curvatures(X);
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> std::string MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Name() const
{
    return Static_Name();
}
//#####################################################################
// Function Static_Name
//#####################################################################
template<class TV> std::string MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Static_Name()
{
    return STRING_UTILITIES::string_sprintf("MULTIBODY_LEVELSET_IMPLICIT_OBJECT<T,VECTOR<T,%d> >",TV::dimension);
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV> std::string MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Extension() const
{
    return Static_Extension();
}
//#####################################################################
// Function Static_Extension
//#####################################################################
template<class TV> std::string MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Static_Extension()
{
    return TV::dimension==2?"mphi2d":"mphi";
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Read(TYPED_ISTREAM& input) // TODO -- fix to read/write levelsets
{
    /*Read_Binary(input,levelsets);*/
    Update_Box();
    Update_Minimum_Cell_Size();
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,*levelsets);
}
namespace PhysBAM{
template class MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,1> >;
template class MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,2> >;
template class MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,3> >;
template class MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,1> >;
template class MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,2> >;
template class MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,3> >;
}
