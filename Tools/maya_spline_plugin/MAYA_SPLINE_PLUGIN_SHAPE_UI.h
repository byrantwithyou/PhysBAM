//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MAYA_SPLINE_PLUGIN_SHAPE_UI__
#define __MAYA_SPLINE_PLUGIN_SHAPE_UI__

#include "MAYA_SPLINE_PLUGIN_SHAPE.h"
#include <maya/MPxSurfaceShapeUI.h>

namespace PhysBAM{

class MAYA_SPLINE_PLUGIN_SHAPE_UI:public MPxSurfaceShapeUI
{
public:
    MAYA_SPLINE_PLUGIN_SHAPE_UI();
    virtual ~MAYA_SPLINE_PLUGIN_SHAPE_UI();

    virtual void getDrawRequests(const MDrawInfo & info,bool objectAndActiveOnly,MDrawRequestQueue& requests);
    virtual void draw(const MDrawRequest & request,M3dView & view) const;
    virtual bool select(MSelectInfo &selectInfo,MSelectionList &selectionList,MPointArray &worldSpaceSelectPts) const;
    // helpers
    void Draw_Vertices(const MDrawRequest & request,M3dView & view) const;
    void Draw_Wireframe(const MDrawRequest & request,M3dView & view) const;
    void Draw_Shaded(const MDrawRequest & request,M3dView & view) const;
    bool Select_Vertices(MSelectInfo &selectInfo,MSelectionList &selectionList,MPointArray &worldSpaceSelectPts) const;

    static void* creator();

private:
    enum {kDrawVertices,kDrawWireframe,kDrawWireframeOnShaded,kDrawSmoothShaded,kDrawFlatShaded,kLastToken};

//#####################################################################
};
}
#endif
