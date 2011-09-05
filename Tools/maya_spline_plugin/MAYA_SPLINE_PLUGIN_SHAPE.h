//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MAYA_SPLINE_PLUGIN_SHAPE__
#define __MAYA_SPLINE_PLUGIN_SHAPE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <maya/MPxSurfaceShape.h>

namespace PhysBAM{
template<class T> class TETRAHEDRALIZED_VOLUME;

class MAYA_SPLINE_PLUGIN_GEOMETRY;

class MAYA_SPLINE_PLUGIN_SHAPE:public MPxSurfaceShape
{
private:
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry;
    GRID_3D<float> cp_grid;
    GRID_3D<float> eval_grid; 
    ARRAYS<VECTOR<VECTOR_3D<float> ,3> > controls_ghost;
    bool debug;
public:
    ARRAYS<int> uvw_iso_mesh;
    TETRAHEDRALIZED_VOLUME<float>* tetrahedralized_volume;

public:
    static MTypeId id;
    bool has_history_on_create;

    // attributes
    static MObject input_surface;
    static MObject output_surface;
    static MObject world_surface;
    static MObject cached_surface;
    static MObject bbox_min_corner;
    static MObject bbox_max_corner;
    static MObject real_world_matrix;

    MAYA_SPLINE_PLUGIN_SHAPE();
    virtual ~MAYA_SPLINE_PLUGIN_SHAPE();
    void postConstructor();

    // compute
    virtual MStatus compute(const MPlug& plug,MDataBlock& datablock); // Overridden
    MStatus Compute_Input_Surface(const MPlug& plug,MDataBlock& datablock); 
    MStatus Compute_Output_Surface(const MPlug& plug,MDataBlock& datablock);
    MStatus Compute_World_Surface(const MPlug& plug,MDataBlock& datablock);
    MStatus Compute_Bounding_Box(MDataBlock& datablock);
    MStatus Apply_Tweaks(MDataBlock& datablock,MAYA_SPLINE_PLUGIN_GEOMETRY* geometry);
    void Update_Tesselation(MAYA_SPLINE_PLUGIN_GEOMETRY* geometry);

    // bounding boxes
    virtual bool isBounded() const;
    virtual MBoundingBox boundingBox() const;

    // getting/setting values
    virtual bool getInternalValue(const MPlug& plug,MDataHandle& handle); // Overridden
    virtual bool setInternalValue(const MPlug& plug, const MDataHandle& handle); // Overridden
    virtual MStatus connectionMade(const MPlug& plug,const MPlug& otherPlug,bool asSrc); // Overridden
    virtual MStatus connectionBroken(const MPlug& plug,const MPlug& otherPlug,bool asSrc); // Overridden
    virtual MStatus shouldSave(const MPlug& plug,bool& result); // Overridden

    virtual MObject geometryData() const; // Overridden

    // overrides to support translate/rotate/scale
    virtual void transformUsing(const MMatrix & mat,const MObjectArray & componentList);
    virtual void transformUsing(const MMatrix& mat,const MObjectArray& componentList,MVertexCachingMode cachingMode,MPointArray* pointCache);
    virtual void tweakUsing(const MMatrix & mat,const MObjectArray & componentList,MVertexCachingMode cachingMode,MPointArray* pointCache,MArrayDataHandle& handle);

    // component mapping
    virtual void componentToPlugs(MObject& component,MSelectionList& list) const; // Overridden
    virtual MPxSurfaceShape::MatchResult matchComponent(const MSelectionList& item,const MAttributeSpecArray& spec,MSelectionList& list); // Overridden
    virtual bool match(const MSelectionMask& mask,const MObjectArray& componentList) const; // Overridden
    virtual MObject createFullVertexGroup() const; // Overridden
    virtual MObject localShapeInAttr() const; // Overridden
    virtual MObject localShapeOutAttr() const; // Overridden
    virtual MObject worldShapeOutAttr() const; // Overridden
    virtual MObject cachedShapeAttr() const; // Overridden
    virtual MPxGeometryIterator* geometryIteratorSetup(MObjectArray&, MObject&,bool forReadOnly=false);
    virtual bool acceptsGeometryIterator(bool writeable=true);
    virtual bool acceptsGeometryIterator(MObject&,bool writeable=true,bool forReadOnly=false);

    // Safe helper setter getters
    bool Value(int point_index,int component_index,double& value) const;
    bool Value(int point_index,MPoint& value) const;
    bool Set_Value(int point_index,int component_index,double value);
    bool Set_Value(int point_index,const MPoint& value);
    MAYA_SPLINE_PLUGIN_GEOMETRY* Cached_Geometry();
    MAYA_SPLINE_PLUGIN_GEOMETRY* Mesh_Geometry();
    void Vertices_Changed();
    MStatus Build_Control_Points(MDataBlock& datablock,int count);

    // plugin registration
    static void* creator(); // Overridden
    static MStatus initialize(); // Overridden

    

//#####################################################################
};
}
#endif
