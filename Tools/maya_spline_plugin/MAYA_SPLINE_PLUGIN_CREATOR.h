//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MAYA_SPLINE_PLUGIN_CREATOR__
#define __MAYA_SPLINE_PLUGIN_CREATOR__

#include <maya/MIntArray.h>
#include <maya/MPointArray.h>
#include <maya/MPxNode.h>
#include <maya/MTypeId.h>
#include <maya/MVectorArray.h>

namespace PhysBAM{

class MAYA_SPLINE_PLUGIN_CREATOR:public MPxNode
{
public:
    static MObject input_mesh;
    static MObject output_surface;
    static MObject size_u,size_v,size_w;
    static MObject zeroed;
    static MTypeId id;

    MAYA_SPLINE_PLUGIN_CREATOR();
    virtual ~MAYA_SPLINE_PLUGIN_CREATOR();

    virtual MStatus compute(const MPlug& plug,MDataBlock& data);
    static void* creator();
    static MStatus initialize();

//#####################################################################
};
}
#endif
