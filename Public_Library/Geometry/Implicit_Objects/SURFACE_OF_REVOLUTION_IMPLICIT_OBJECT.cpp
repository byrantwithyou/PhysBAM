//#####################################################################
// Copyright 2007, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Images/IMAGE.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T_input> SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<T_input>::
SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT(std::string& filename,const T width_input,const T height_input,const T tolerance_input)
    :slice_levelset(slice_grid,slice_phi),width(width_input),height(height_input),tolerance(tolerance_input)
{
    ARRAY<TV,VECTOR<int,2> > slice_image;IMAGE<T>::Read(filename,slice_image);
    slice_grid=GRID<VECTOR<T,2> >(slice_image.Size(),RANGE<VECTOR<T,2> >(VECTOR<T,2>(),VECTOR<T,2>(width,height)));slice_phi.Resize(slice_grid.Domain_Indices());
    for(NODE_ITERATOR<VECTOR<T,2> > iterator(slice_grid);iterator.Valid();iterator.Next()){VECTOR<int,2> node=iterator.Node_Index();
        slice_phi(node)=2*slice_image(node).Average()-1;}
    slice_levelset.Fast_Marching_Method();slice_levelset.Compute_Normals();
    Update_Box();
}
//#####################################################################
// Function Update_Box
//#####################################################################
template<class T_input> void SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<T_input>::
Update_Box()
{
    box=RANGE<TV>(TV(-width,0,-width),TV(width,height,width));
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class T_input> T_input SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<T_input>::
Integration_Step(const T phi) const
{
    return max(phi,tolerance);
}
//#####################################################################
// Function Slice_Location
//#####################################################################
template<class T> static VECTOR<T,2>
Slice_Location(const VECTOR<T,3>& X)
{
    return VECTOR<T,2>(sqrt(sqr(X.x)+sqr(X.z)),X.y);
}
//#####################################################################
// Function operator()
//#####################################################################
template<class T_input> T_input SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<T_input>::
operator()(const TV& X) const
{
    return slice_levelset.Phi(Slice_Location(X));
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T_input> auto SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<T_input>::
Normal(const TV& X,const int aggregate) const -> TV
{
    assert((aggregate>=1 && aggregate<=6) || aggregate==-1);if(aggregate!=-1) return box.Normal(aggregate);
    VECTOR<T,2> slice_X=Slice_Location(X),horizontal_direction=VECTOR<T,2>(X.x,X.z).Normalized(); // cse should remove duplicate sqrt
    VECTOR<T,2> slice_normal=slice_levelset.Normal(slice_X);
    VECTOR<T,2> horizontal_normal=slice_normal.x*horizontal_direction;
    return TV(horizontal_normal.x,slice_normal.y,horizontal_normal.y);
}
//#####################################################################
// Function Read
//#####################################################################
template<class T_input> void SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<T_input>::
Read(TYPED_ISTREAM& input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Write
//#####################################################################
template<class T_input> void SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<T_input>::
Write(TYPED_OSTREAM& output) const
{
    PHYSBAM_FATAL_ERROR();
}
template class SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<double>;
template class SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT<float>;
}
