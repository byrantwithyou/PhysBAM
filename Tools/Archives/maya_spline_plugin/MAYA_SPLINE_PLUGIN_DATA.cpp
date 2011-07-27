//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <iostream>
#include "MAYA_SPLINE_PLUGIN_DATA.h"
#include "MAYA_SPLINE_PLUGIN_ITERATOR.h"

#include <maya/MArgList.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnTripleIndexedComponent.h>

using namespace PhysBAM;

const MTypeId MAYA_SPLINE_PLUGIN_DATA::id(0x80777);
const MString MAYA_SPLINE_PLUGIN_DATA::typeName("MAYA_SPLINE_PLUGIN_DATA");
//#####################################################################
// Function MAYA_SPLINE_PLUGIN_DATA
//#####################################################################
MAYA_SPLINE_PLUGIN_DATA::
MAYA_SPLINE_PLUGIN_DATA() 
    :geometry(0)
{
    geometry=new MAYA_SPLINE_PLUGIN_GEOMETRY;
}
//#####################################################################
// Function ~MAYA_SPLINE_PLUGIN_DATA
//#####################################################################
MAYA_SPLINE_PLUGIN_DATA::
~MAYA_SPLINE_PLUGIN_DATA()
{
    delete geometry;
}
//#####################################################################
// Function copy
//#####################################################################
void MAYA_SPLINE_PLUGIN_DATA::
copy(const MPxData& raw_other)
{
    const MAYA_SPLINE_PLUGIN_DATA& other=(const MAYA_SPLINE_PLUGIN_DATA&)raw_other;
    *geometry=*other.geometry;
}
//#####################################################################
// Function typeId
//#####################################################################
MTypeId MAYA_SPLINE_PLUGIN_DATA::
typeId() const
{
    return MAYA_SPLINE_PLUGIN_DATA::id;
}
//#####################################################################
// Function name
//#####################################################################
MString MAYA_SPLINE_PLUGIN_DATA::
name() const
{
    return MAYA_SPLINE_PLUGIN_DATA::typeName;
}
//#####################################################################
// Function creator
//#####################################################################
void* MAYA_SPLINE_PLUGIN_DATA::
creator()
{
    return new MAYA_SPLINE_PLUGIN_DATA;
}
//#####################################################################
// Function iterator
//#####################################################################
MPxGeometryIterator* MAYA_SPLINE_PLUGIN_DATA::
iterator(MObjectArray& component_list,MObject& component,bool use_components)
{
    if(use_components) return new MAYA_SPLINE_PLUGIN_ITERATOR(geometry,component_list);
    else return new MAYA_SPLINE_PLUGIN_ITERATOR(geometry,component);
}
//#####################################################################
// Function iterator
//#####################################################################
MPxGeometryIterator* MAYA_SPLINE_PLUGIN_DATA::
iterator(MObjectArray& component_list,MObject& component,bool use_components,bool world) const
{
    // TODO: Figure out what world does
    if(use_components) return new MAYA_SPLINE_PLUGIN_ITERATOR(geometry,component_list);
    else return new MAYA_SPLINE_PLUGIN_ITERATOR(geometry,component);
}
//#####################################################################
// Function updateCompleteVertexGroup
//#####################################################################
bool MAYA_SPLINE_PLUGIN_DATA::updateCompleteVertexGroup(MObject& component) const
{
    MStatus stat;
    MFnSingleIndexedComponent fn_component(component,&stat);
    if(stat && geometry && fn_component.isComplete()){    
        int vertices=geometry->Get_Number_Vertices(),max_vertices;
        fn_component.getCompleteData(max_vertices);
        if (vertices != max_vertices){fn_component.setCompleteData(vertices);return true;}}
    return false;
}
//#####################################################################
