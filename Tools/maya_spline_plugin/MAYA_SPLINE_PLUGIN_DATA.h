//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef _MAYA_SPLINE_PLUGIN_DATA
#define _MAYA_SPLINE_PLUGIN_DATA

#include "MAYA_SPLINE_PLUGIN_GEOMETRY.h"
#include <maya/MPxGeometryData.h>
#include <maya/MString.h>
#include <maya/MTypeId.h>

namespace PhysBAM{

class MAYA_SPLINE_PLUGIN_DATA:public MPxGeometryData
{
public:
    static const MString typeName;
    static const MTypeId id;
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry;
    
    MAYA_SPLINE_PLUGIN_DATA();
    virtual ~MAYA_SPLINE_PLUGIN_DATA();
    virtual    void copy ( const MPxData& );
    virtual MTypeId typeId() const;
    virtual MString name() const;
    virtual MPxGeometryIterator* iterator(MObjectArray& component_list,MObject& component,bool use_components);
    virtual MPxGeometryIterator* iterator(MObjectArray& component_list,MObject& component,bool use_components,bool world) const;
    virtual bool updateCompleteVertexGroup(MObject& component) const;
    static void* creator();

//#####################################################################
};
}
#endif
