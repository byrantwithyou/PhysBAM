//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MAYA_SPLINE_PLUGIN_ITERATOR__
#define __MAYA_SPLINE_PLUGIN_ITERATOR__

#include "MAYA_SPLINE_PLUGIN_GEOMETRY.h"
#include <maya/MPoint.h>
#include <maya/MPxGeometryIterator.h>

namespace PhysBAM{

class MAYA_SPLINE_PLUGIN_ITERATOR:public MPxGeometryIterator
{
public:
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry;
    
    MAYA_SPLINE_PLUGIN_ITERATOR(void* userGeometry,MObjectArray& components);
    MAYA_SPLINE_PLUGIN_ITERATOR(void* userGeometry,MObject& components);
    virtual void reset();
    virtual MPoint point() const;
    virtual void setPoint(const MPoint&) const;
    virtual int    iteratorCount() const;
    virtual bool hasPoints() const;

//#####################################################################
};
}
#endif
