//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include "MAYA_SPLINE_PLUGIN_ITERATOR.h"
using namespace PhysBAM;

//#####################################################################
// Function MAYA_SPLINE_PLUGIN_ITERATOR
//#####################################################################
MAYA_SPLINE_PLUGIN_ITERATOR::
MAYA_SPLINE_PLUGIN_ITERATOR(void* geometry_input,MObjectArray& components_input) 
    :MPxGeometryIterator(geometry_input,components_input),geometry((MAYA_SPLINE_PLUGIN_GEOMETRY*)geometry_input)
{
    reset();
}
//#####################################################################
// MAYA_SPLINE_PLUGIN_ITERATOR
//#####################################################################
MAYA_SPLINE_PLUGIN_ITERATOR::
MAYA_SPLINE_PLUGIN_ITERATOR(void* geometry_input,MObject& components_input) 
    :MPxGeometryIterator(geometry_input,components_input),geometry((MAYA_SPLINE_PLUGIN_GEOMETRY*)geometry_input)
{
    reset();
}
//#####################################################################
// Function reset
//#####################################################################
void MAYA_SPLINE_PLUGIN_ITERATOR::
reset()
{
    MPxGeometryIterator::reset();
    setCurrentPoint(0);
    if(geometry) setMaxPoints(geometry->Get_Number_Vertices());
}
//#####################################################################
// Function point
//#####################################################################
MPoint MAYA_SPLINE_PLUGIN_ITERATOR::
point() const
{
    MPoint maya_point;
    if(geometry){VECTOR_3D<float> point=geometry->controls.array[index()];maya_point.x=point.x;maya_point.y=point.y;maya_point.z=point.z;}
    return maya_point;
}
//#####################################################################
// Function setPoint
//#####################################################################
void MAYA_SPLINE_PLUGIN_ITERATOR::
setPoint(const MPoint& pnt) const
{
    if(geometry) geometry->controls.array[index()]=VECTOR_3D<float>(pnt.x,pnt.y,pnt.z);
}
//#####################################################################
// Function iteratorCount
//#####################################################################
int MAYA_SPLINE_PLUGIN_ITERATOR::
iteratorCount() const
{
    return geometry->Get_Number_Vertices();
}
//#####################################################################
// Function hasPoints
//#####################################################################
bool MAYA_SPLINE_PLUGIN_ITERATOR::
hasPoints() const
{
    return true; // we have points
}
//#####################################################################
