//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/CATMULL_ROM_SPLINE_INTERPOLATION.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include "MAYA_PLUGIN_MACROS.h"
#include "MAYA_SPLINE_PLUGIN_DATA.h"
#include "MAYA_SPLINE_PLUGIN_GEOMETRY.h"
#include "MAYA_SPLINE_PLUGIN_ITERATOR.h"
#include "MAYA_SPLINE_PLUGIN_SHAPE.h"
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>

#include <maya/MAttributeIndex.h>
#include <maya/MAttributeSpec.h>
#include <maya/MAttributeSpecArray.h>
#include <maya/MDagPath.h>
#include <maya/MFnAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnPluginData.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MMatrix.h>
#include <maya/MObjectArray.h>
#include <maya/MPointArray.h>
#include <maya/MPxNode.h>
#include <maya/MSelectionList.h>
#include <maya/MSelectionMask.h>
#include <maya/MVector.h>
 
using namespace PhysBAM; 

//#####################################################################
// static attributes and id
//#####################################################################
MObject MAYA_SPLINE_PLUGIN_SHAPE::input_surface;
MObject MAYA_SPLINE_PLUGIN_SHAPE::output_surface;
MObject MAYA_SPLINE_PLUGIN_SHAPE::cached_surface;
MObject MAYA_SPLINE_PLUGIN_SHAPE::world_surface;
MObject MAYA_SPLINE_PLUGIN_SHAPE::bbox_min_corner;
MObject MAYA_SPLINE_PLUGIN_SHAPE::bbox_max_corner;
MObject MAYA_SPLINE_PLUGIN_SHAPE::real_world_matrix;
MTypeId MAYA_SPLINE_PLUGIN_SHAPE::id(0x80099);
//#####################################################################
// Function MAYA_SPLINE_PLUGIN_SHAPE
//#####################################################################
MAYA_SPLINE_PLUGIN_SHAPE::
MAYA_SPLINE_PLUGIN_SHAPE()
    :tetrahedralized_volume(0),debug(false)
{}
//#####################################################################
// Function ~MAYA_SPLINE_PLUGIN_SHAPE
//#####################################################################
MAYA_SPLINE_PLUGIN_SHAPE::
~MAYA_SPLINE_PLUGIN_SHAPE()
{
    if(tetrahedralized_volume){tetrahedralized_volume->Destroy_Data();delete tetrahedralized_volume;}
}
//#####################################################################
// Function postConstructor
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE::
postConstructor()
{
    setRenderable(true); // allow shading group assignment
    has_history_on_create=false; // no history connecte
}
//#####################################################################
// Function compute
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
compute(const MPlug& plug,MDataBlock& datablock)
{
    if(debug) std::cerr << "MAYA_SPLINE_PLUGIN_SHAPE::compute : plug "<<plug.info()<<std::endl;
    if(plug==output_surface || plug==cached_surface) return Compute_Output_Surface(plug,datablock);
    else if(plug==world_surface) return Compute_World_Surface(plug,datablock);
    else return MS::kUnknownParameter;
}
//#####################################################################
// Function Compute_Output_Surface
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
Compute_Output_Surface(const MPlug& plug,MDataBlock& datablock)
{
    MStatus stat;
    // Copy input surface to cached surface
    if(!Compute_Input_Surface(plug,datablock)) return MS::kFailure;
    // Get cached geometry
    MDataHandle cached_handle=datablock.outputValue(cached_surface,&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: Compute_Output_Surface getting cached_surface");
    MAYA_SPLINE_PLUGIN_DATA* cached_data=(MAYA_SPLINE_PLUGIN_DATA*)cached_handle.asPluginData();
    if(cached_data==NULL) std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE: Compute_Output_Surface got NULL cached handle"<<std::endl;
    datablock.setClean(plug);
    // apply vertex offsets
    if(has_history_on_create) Apply_Tweaks(datablock,cached_data->geometry);
    else datablock.inputArrayValue(mControlPoints,&stat).setAllClean();
    // create the output surface
    MFnPluginData fn_data_creator;
    MTypeId temp_id(MAYA_SPLINE_PLUGIN_DATA::id);
    MObject new_data_object=fn_data_creator.create(temp_id,&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: compute_output_surface: error creating MAYA_SPLINE_PLUGIN_DATA");
    MAYA_SPLINE_PLUGIN_DATA* new_data=(MAYA_SPLINE_PLUGIN_DATA*)fn_data_creator.data(&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: compute_output_surface: error getting MAYA_SPLINE_PLUGIN_DATA");
    // copy the cache over
    if(cached_data) *(new_data->geometry)=*(cached_data->geometry);
    Update_Tesselation(new_data->geometry);
    MDataHandle output_handle=datablock.outputValue(output_surface);
    output_handle.set(new_data);

    stat=Compute_Bounding_Box(datablock);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_DATA: Compute_Output_Surface: Compute_Bounding_Box");
    return stat;
}
//#####################################################################
// Function Compute_World_Surface
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
Compute_World_Surface(const MPlug& plug,MDataBlock& datablock)
{
    MStatus stat;
    //std::cerr<<"In compute world surface"<<std::endl;
    Compute_Output_Surface(plug,datablock);
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=Mesh_Geometry();
    datablock.setClean(plug);
    MDataHandle matrix_handle=datablock.inputValue(real_world_matrix);
    MMatrix world_matrix=matrix_handle.asMatrix();
    MATRIX_4X4<float> matrix(world_matrix(0,0),world_matrix(0,1),world_matrix(0,2),world_matrix(0,3),
                         world_matrix(1,0),world_matrix(1,1),world_matrix(1,2),world_matrix(1,3),
                         world_matrix(2,0),world_matrix(2,1),world_matrix(2,2),world_matrix(2,3),
                         world_matrix(3,0),world_matrix(3,1),world_matrix(3,2),world_matrix(3,3));
    // create the output surface
    MFnPluginData fn_data_creator;
    MTypeId temp_id(MAYA_SPLINE_PLUGIN_DATA::id);
    MObject new_data_object=fn_data_creator.create(temp_id,&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: compute_output_surface: error creating MAYA_SPLINE_PLUGIN_DATA");
    MAYA_SPLINE_PLUGIN_DATA* new_data=(MAYA_SPLINE_PLUGIN_DATA*)fn_data_creator.data(&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: compute_output_surface: error getting MAYA_SPLINE_PLUGIN_DATA");
    std::cerr<<matrix<<std::endl;
    // copy the cache over
    if(geometry){
        *(new_data->geometry)=*geometry;
        // transform all control points into world space
        ARRAYS<VECTOR<VECTOR_3D<float> ,3> >& controls=new_data->geometry->controls;
        int vertex_count=new_data->geometry->Get_Number_Vertices();
        for(int i=0;i<vertex_count;i++) controls.array[i]=matrix*controls.array[i];}
    else std::cerr<<"got null surface on geometry"<<std::endl;
        
    MDataHandle world_handle=datablock.outputValue(world_surface);
    world_handle.set(new_data);

    return MS::kSuccess;
}
//#####################################################################
// Function Compute_World_Surface
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
Compute_Input_Surface(const MPlug& plug,MDataBlock& datablock)
{
    //std::cerr<<"In compute input surface"<<std::endl;
    MStatus stat=MS::kSuccess;
    if(has_history_on_create){ // has input surface, copy it to cache
        MDataHandle input_handle=datablock.inputValue(input_surface,&stat);
        MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: Compute_Input_Surface error getting input_surface");
        MAYA_SPLINE_PLUGIN_DATA* data=(MAYA_SPLINE_PLUGIN_DATA*)input_handle.asPluginData();
        if(!data){std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE: Compute_Input_Surface NULL input_surface data found"<<std::endl;return stat;}
        MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=((MAYA_SPLINE_PLUGIN_DATA*)data)->geometry;
        // create the cached surface
        MFnPluginData fn_data_creator;
        MTypeId temp_id(MAYA_SPLINE_PLUGIN_DATA::id);
        MObject new_data_object=fn_data_creator.create(temp_id,&stat);
        MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: compute_input_surface: error creating MAYA_SPLINE_PLUGIN_DATA");
        MAYA_SPLINE_PLUGIN_DATA* new_data=(MAYA_SPLINE_PLUGIN_DATA*)fn_data_creator.data(&stat);
        MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: compute_input_surface: error getting MAYA_SPLINE_PLUGIN_DATA");
        *(new_data->geometry)=*geometry;
        MDataHandle cached_handle=datablock.outputValue(cached_surface,&stat);
        MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: compute_input_surface: error getting cached_surface");
        cached_handle.set(new_data_object);}
    return stat;
}
//#####################################################################
// Function Compute_Bounding_Box
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
Compute_Bounding_Box(MDataBlock& datablock)
{
    //std::cerr<<"In compute bounding box"<<std::endl;
    MStatus stat=MS::kSuccess;
    // get bounding box attributes
    MDataHandle min_handle=datablock.outputValue(bbox_min_corner);
    MDataHandle max_handle=datablock.outputValue(bbox_max_corner);
    double3& min_corner=min_handle.asDouble3();
    double3& max_corner=max_handle.asDouble3();
    // get mesh pointer
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=Mesh_Geometry();
    int vertex_count=geometry->Get_Number_Vertices();
    if(vertex_count==0) return stat;
    // use physbam BOX_3D to compute bounding box
    BOX_3D<float> box;
    if(tetrahedralized_volume){
        box.Reset_Bounds(tetrahedralized_volume->particles.X(1));
        for(int i=2;i<tetrahedralized_volume->particles.number;i++) box.Enlarge_To_Include_Point(tetrahedralized_volume->particles.X(i));}
    else{
        box.Reset_Bounds(geometry->controls.array[0]);
        for(int i=1;i<vertex_count;i++) box.Enlarge_To_Include_Point(geometry->controls.array[i]);}
    
    // convert to maya attribute
    VECTOR_3D<float> min_corner_physbam=box.Minimum_Corner();
    VECTOR_3D<float> max_corner_physbam=box.Maximum_Corner();
    min_corner[0]=min_corner_physbam.x;min_corner[1]=min_corner_physbam.y;min_corner[2]=min_corner_physbam.z;
    max_corner[0]=max_corner_physbam.x;max_corner[1]=max_corner_physbam.y;max_corner[2]=max_corner_physbam.z;
    min_handle.setClean();max_handle.setClean();
    childChanged(MPxSurfaceShape::kBoundingBoxChanged);
    return stat;
}
//#####################################################################
// Function Apply_Tweaks
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
Apply_Tweaks(MDataBlock& datablock,MAYA_SPLINE_PLUGIN_GEOMETRY* geometry)
{
    MStatus stat;
    MArrayDataHandle cp_handle=datablock.inputArrayValue(mControlPoints,&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: Apply_Tweaks get cpHandle failed");
    int elements=cp_handle.elementCount();
    for(int i=0;i<elements;i++){
        int element_index=cp_handle.elementIndex();
        MDataHandle point_handle=cp_handle.outputValue();
        double3& point=point_handle.asDouble3();
        geometry->controls.array[element_index]+=VECTOR_3D<float>(point[0],point[1],point[2]);
        cp_handle.next();}
}
//#####################################################################
// Function Update_Tesselation
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE::
Update_Tesselation(MAYA_SPLINE_PLUGIN_GEOMETRY* geometry)
{
    if(tetrahedralized_volume && (cp_grid.m != geometry->grid.m || cp_grid.n != geometry->grid.n || cp_grid.mn != geometry->grid.mn)){
        tetrahedralized_volume->Destroy_Data();delete tetrahedralized_volume;tetrahedralized_volume=0;}
    if(!tetrahedralized_volume){
        VECTOR_3D<int> divisions(5,5,5);
        cp_grid=geometry->grid;eval_grid.Initialize(cp_grid.number_of_cells_x*divisions.x+1,cp_grid.number_of_cells_y*divisions.y+1,cp_grid.number_of_cells_z*divisions.z+1,cp_grid.xmin,cp_grid.xmax,cp_grid.ymin,cp_grid.ymax,cp_grid.zmin,cp_grid.zmax);
        tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<float>::Create();tetrahedralized_volume->Initialize_Cube_Mesh_And_Particles(eval_grid);
        tetrahedralized_volume->Initialize_Triangulated_Surface();tetrahedralized_volume->triangulated_surface->triangle_mesh.Initialize_Segment_Mesh();
        controls_ghost.Resize(cp_grid,2);
        ARRAYS<VECTOR<int,3> > grid_to_particle(1,eval_grid.m,1,eval_grid.n,1,eval_grid.mn);int index=1;
        for(int ij=0;ij<eval_grid.mn;ij++) for(int j=0;j<eval_grid.n;j++) for(int i=0;i<eval_grid.m;i++){
            grid_to_particle(i,j,ij)=index;index++;}
        // make an iso line mesh
        uvw_iso_mesh.Resize(2,0);// should be cp_grid.m*cp_grid.n*eval_grid.mn+cp_grid.m*eval_grid.n+eval_grid.m*cp_grid.n*cp_grid.mn);
        for(int i=0;i<cp_grid.m;i++) for(int j=0;j<cp_grid.n;j++){
            int ii=(i-1)*divisions.x+1,jj=(j-1)*divisions.y+1;for(int ij=2;ij<=eval_grid.mn;ij++) uvw_iso_mesh.Append(grid_to_particle(ii,jj,ij-1),grid_to_particle(ii,jj,ij));}
        for(int j=0;j<cp_grid.n;j++) for(int ij=0;ij<cp_grid.mn;ij++){
            int jj=(j-1)*divisions.y+1,iijj=(ij-1)*divisions.z+1;for(int i=2;i<=eval_grid.m;i++) uvw_iso_mesh.Append(grid_to_particle(i-1,jj,iijj),grid_to_particle(i,jj,iijj));}
        for(int i=0;i<cp_grid.m;i++) for(int ij=0;ij<cp_grid.mn;ij++){
            int ii=(i-1)*divisions.x+1,iijj=(ij-1)*divisions.z+1;for(int j=2;j<=eval_grid.n;j++) uvw_iso_mesh.Append(grid_to_particle(ii,j-1,iijj),grid_to_particle(ii,j,iijj));}}
    ARRAYS<VECTOR<VECTOR_3D<float> ,3> >::put(geometry->controls,controls_ghost);
    CATMULL_ROM_SPLINE_INTERPOLATION<float,VECTOR_3D<float> > interpolation;int index=1;
    for(int ij=0;ij<eval_grid.mn;ij++) for(int j=0;j<eval_grid.n;j++) for(int i=0;i<eval_grid.m;i++) {
        tetrahedralized_volume->particles.X(index)=interpolation.Clamped_To_Array(cp_grid,controls_ghost,eval_grid.X(i,j,ij));index++;}
}
//#####################################################################
// Function isBounded
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
isBounded() const
{ 
    return true;
}
//#####################################################################
// Function boundingBox
//#####################################################################
MBoundingBox MAYA_SPLINE_PLUGIN_SHAPE::
boundingBox() const
{
    MObject this_node=thisMObject();
    MPlug min_plug(this_node,bbox_min_corner);MPlug max_plug(this_node,bbox_max_corner);
    MObject min_object;MObject max_object;
    min_plug.getValue(min_object);max_plug.getValue(max_object);
    // get data from attributes
    double3 min_corner,max_corner;
    MFnNumericData fnData;
    fnData.setObject(min_object);
    fnData.getData(min_corner[0], min_corner[1], min_corner[2]);
    fnData.setObject(max_object);
    fnData.getData(max_corner[0], max_corner[1], max_corner[2]);
    // make maya points and return maya bounding box
    MPoint min_corner_point(min_corner[0], min_corner[1], min_corner[2]);
    MPoint max_corner_point(max_corner[0], max_corner[1], max_corner[2]);
    return MBoundingBox(min_corner_point,max_corner_point);
}    
//#####################################################################
// Function getInternalValue
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
getInternalValue(const MPlug& plug,MDataHandle& result)
{
    bool is_ok=true;
    if(plug==mControlPoints || plug==mControlValueX || plug==mControlValueY || plug==mControlValueZ){
        if(has_history_on_create) return MPxNode::getInternalValue(plug, result); // these are tweaks
        else{
            double val=0.0;
            if(plug==mControlPoints && !plug.isArray()){
                MPoint pnt;
                int index=plug.logicalIndex();
                Value(index,pnt);result.set(pnt[0],pnt[1],pnt[2]);}
            else if(plug==mControlValueX){
                MPlug parentPlug=plug.parent();int index=parentPlug.logicalIndex();
                Value(index,0,val);result.set(val);}
            else if(plug==mControlValueY){
                MPlug parentPlug=plug.parent();int index=parentPlug.logicalIndex();
                Value(index,1,val);result.set(val);}
            else if(plug==mControlValueZ){
                MPlug parentPlug=plug.parent();int index=parentPlug.logicalIndex();
                Value(index,2,val);result.set(val);}}}
    else if(plug==mHasHistoryOnCreate) result.set(has_history_on_create);
    else is_ok=MPxSurfaceShape::getInternalValue(plug,result);
    return is_ok;
}
//#####################################################################
// Function setInternalValue
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
setInternalValue(const MPlug& plug,const MDataHandle& handle)
{
    bool is_ok=true;

    if(plug==mControlPoints || plug==mControlValueX || plug==mControlValueY || plug==mControlValueZ){
        if(has_history_on_create){Vertices_Changed();return MPxNode::setInternalValue(plug,handle);}
        else{
            if(plug==mControlPoints && !plug.isArray()){
                int index = plug.logicalIndex();
                MPoint point;double3& ptData=handle.asDouble3();
                point.x=ptData[0];point.y=ptData[1];point.z=ptData[2];
                Set_Value(index,point);}
            else if(plug==mControlValueX){
                MPlug parentPlug=plug.parent();
                int index=parentPlug.logicalIndex();
                Set_Value(index,0,handle.asDouble());}
            else if(plug==mControlValueY){
                MPlug parentPlug=plug.parent();
                int index=parentPlug.logicalIndex();
                Set_Value(index,1,handle.asDouble());}
            else if(plug==mControlValueZ){
                MPlug parentPlug=plug.parent();
                int index=parentPlug.logicalIndex();
                Set_Value(index,2,handle.asDouble());}}}
    else if(plug==mHasHistoryOnCreate) has_history_on_create=handle.asBool();
    else is_ok=MPxSurfaceShape::setInternalValue(plug,handle);
    return is_ok;
}
//#####################################################################
// Function connectionMade
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
connectionMade(const MPlug& plug,const MPlug& otherPlug,bool asSrc)
{
    std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE: connectionMade "<<plug.info()<<" other="<<otherPlug.info()<<std::endl;
    if(plug==input_surface){
        MStatus stat;
        MObject this_object=thisMObject();
        MPlug historyPlug(this_object,mHasHistoryOnCreate);
        stat=historyPlug.setValue(true);
        if(stat!=MS::kSuccess){std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE: connectionMade: setValue(mHasHistoryOnCreate)"<<std::endl;return MS::kFailure;}}
    return MPxNode::connectionMade(plug,otherPlug,asSrc);
}
//#####################################################################
// Function connectionBroken
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
connectionBroken(const MPlug& plug,const MPlug& otherPlug,bool asSrc)
{
    if(plug==input_surface){
        MStatus stat;
        MObject this_object=thisMObject();
        MPlug historyPlug(this_object,mHasHistoryOnCreate);
        stat=historyPlug.setValue(false);
        if(stat!=MS::kSuccess){std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE: connectionBroken: setValue(mHasHistoryOnCreate)"<<std::endl;return MS::kFailure;}}
    return MPxNode::connectionBroken(plug,otherPlug,asSrc);
}
//#####################################################################
// Function shouldSave
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
shouldSave(const MPlug& plug,bool& result)
{
    MStatus status=MS::kSuccess;
    if(plug==mControlPoints || plug==mControlValueX || plug==mControlValueY || plug==mControlValueZ){
        if(has_history_on_create)  status=MPxNode::shouldSave(plug,result); // only write tweaks if they are different, based on default implementation
        else result=false;}
    else if(plug==cached_surface){
        if(has_history_on_create) result=false;
        else{
            MObject data;
            status=plug.getValue(data);
            MCHECKERROR(status,"MAYA_SPLINE_PLUGIN_SHAPE: shouldSave: getValue");
            result=!data.isNull();}}
    else status=MPxNode::shouldSave(plug,result);
    return status;
}
//#####################################################################
// Function geometryData
//#####################################################################
MObject MAYA_SPLINE_PLUGIN_SHAPE::
geometryData() const
{
    MDataBlock datablock=((MAYA_SPLINE_PLUGIN_SHAPE*)this)->forceCache();
    MDataHandle handle=datablock.inputValue(input_surface);
    return handle.data();
}
//#####################################################################
// Function transformUsing
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE::
transformUsing(const MMatrix& mat,const MObjectArray& componentList)
{
    transformUsing(mat,componentList,MPxSurfaceShape::kNoPointCaching,NULL);
}
//#####################################################################
// Function transformUsing
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE::
transformUsing(const MMatrix& mat,const MObjectArray& componentList,MVertexCachingMode cachingMode,MPointArray* pointCache)
{
    MStatus stat;
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=Mesh_Geometry();
    bool save_points=(cachingMode==MPxSurfaceShape::kSavePoints);
    unsigned int i=0,j=0;
    unsigned int length=componentList.length();

    if(cachingMode==MPxSurfaceShape::kRestorePoints){
        // restore points from the cache
        unsigned int cache_length=pointCache->length();
        if(length>0){
            for(i=0;i<length && j<cache_length;i++){
                MObject component=componentList[i];
                MFnSingleIndexedComponent fn_component(component);
                int element_count=fn_component.elementCount();
                for(int index=0;index<element_count && j<cache_length;index++,j++){
                    int element_index=fn_component.element(index);
                    MPoint maya_point=(*pointCache)[j];
                    VECTOR_3D<float> point(maya_point.x,maya_point.y,maya_point.z);
                    geometry->controls.array[element_index]=point;}}}
        else{ // copy back the entire surface
            length=geometry->Get_Number_Vertices();
            for(unsigned int index=0;index<length && j<cache_length;index++,j++){
                MPoint maya_point=(*pointCache)[j];
                VECTOR_3D<float> point(maya_point.x,maya_point.y,maya_point.z);
                geometry->controls.array[index]=point;}}}
    else{
        if(length>0){
            for(i=0;i<length;i++){
                MObject component=componentList[i];
                MFnSingleIndexedComponent fn_component(component);
                int element_count=fn_component.elementCount();
                for(int index=0;index<element_count;index++){
                    int element_index=fn_component.element(index);
                    VECTOR_3D<float> phys_point=geometry->controls.array[element_index];
                    MPoint point(phys_point.x,phys_point.y,phys_point.z);
                    if(save_points) pointCache->append(point);
                    point*=mat;
                    geometry->controls.array[element_index]=VECTOR_3D<float>(point.x,point.y,point.z);}}}
        else{
            length=geometry->Get_Number_Vertices();
            for(unsigned int index=0;index<length;index++){
                VECTOR_3D<float> phys_point=geometry->controls.array[index];
                MPoint point(phys_point.x,phys_point.y,phys_point.z);
                if(save_points) pointCache->append(point);
                point*=mat;
                geometry->controls.array[index]=VECTOR_3D<float>(point.x,point.y,point.z);}}}

    MDataBlock datablock=forceCache();
    MDataHandle cached_handle=datablock.outputValue(cached_surface,&stat);
    MCHECKERRORNORET(stat,"MAYA_SPLINE_PLUGIN_SHAPE: transformUsing failed to get cache_surface");
    MAYA_SPLINE_PLUGIN_DATA* cached=(MAYA_SPLINE_PLUGIN_DATA*)cached_handle.asPluginData();
    MDataHandle data_handle=datablock.outputValue(mControlPoints,&stat);
    MCHECKERRORNORET(stat,"MAYA_SPLINE_PLUGIN_SHAPE: transformUsing get data_handle for mControlPoints");
    if(has_history_on_create && cached != NULL){ // must create deltas and store on control points
        stat=Build_Control_Points(datablock,geometry->Get_Number_Vertices());
        MCHECKERRORNORET(stat,"MAYA_SPLINE_PLUGIN_SHAPE: transformUsing trying to create control points");
        MArrayDataHandle cp_handle(data_handle,&stat); 
        MCHECKERRORNORET(stat,"MAYA_SPLINE_PLUGIN_SHAPE: transformUsing get array data_handle for mControlPoints");
        for(int i=0;i<length;i++){
            MObject component=componentList[i];
            MFnSingleIndexedComponent fn_component(component);
            int element_count=fn_component.elementCount();
            for(int index=0;index<element_count;index++){
                int element_index=fn_component.element(index);
                cp_handle.jumpToElement(element_index);
                MDataHandle point_handle=cp_handle.outputValue();
                double3& point=point_handle.asDouble3();
                VECTOR_3D<float> phys_offset=geometry->controls.array[element_index]-cached->geometry->controls.array[element_index];
                point[0]+=phys_offset.x;point[1]+=phys_offset.y;point[2]+=phys_offset.z;}}}
    if(!cached) std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE: transformUsing NULL cachedSurface found"<<std::endl;
    else *(cached->geometry)=*geometry;

    MPlug pCPs(thisMObject(),mControlPoints);
    pCPs.setValue(data_handle);
    Compute_Bounding_Box(datablock);
    childChanged(MPxSurfaceShape::kBoundingBoxChanged);
}
//#####################################################################
// Function tweakUsing
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE::
tweakUsing(const MMatrix & mat,const MObjectArray & componentList,MVertexCachingMode cachingMode,MPointArray* pointCache,MArrayDataHandle& handle)
{
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=Mesh_Geometry();
    bool save_points=(cachingMode==MPxSurfaceShape::kSavePoints);
    bool update_points=(cachingMode==MPxSurfaceShape::kUpdatePoints);
    MArrayDataBuilder builder=handle.builder();
    MPoint delta,current_point,new_point;
    unsigned int i=0,length=componentList.length(),cache_index=0,cache_length=(pointCache?pointCache->length():0);
    if(cachingMode==MPxSurfaceShape::kRestorePoints){
        if(length>0) for(i=0;i<length;i++){ // through each component
            MObject component=componentList[i];
            MFnSingleIndexedComponent fn_component(component);
            int element_count=fn_component.elementCount();
            for(int index=0;index<element_count && cache_index<cache_length;index++,cache_index++){
                int element_index=fn_component.element(index);
                double3& pt=builder.addElement(element_index).asDouble3();
                MPoint& cache_point=(*pointCache)[cache_index];
                pt[0]+=cache_point.x;pt[1]+=cache_point.y;pt[2]+=cache_point.z;}}
        else{
            length=geometry->Get_Number_Vertices();
            for(unsigned int index=0;index<length && index<cache_length;index++,cache_index++){
                double3& pt=builder.addElement(index).asDouble3();
                MPoint& cache_point=(*pointCache)[cache_index];
                pt[0]+=cache_point.x;pt[1]+=cache_point.y;pt[2]+=cache_point.z;}}}
    else{
        if(length>0) for(i=0;i<length;i++){
            MObject component=componentList[i];
            MFnSingleIndexedComponent fn_component(component);
            int element_count=fn_component.elementCount();
            if(save_points) pointCache->setSizeIncrement(element_count);
            for(int index=0;index<element_count;index++){
                int element_index=fn_component.element(index);
                double3& pt=builder.addElement(element_index).asDouble3();
                VECTOR_3D<float> phys_point=geometry->controls.array[element_index];
                MPoint point(phys_point.x,phys_point.y,phys_point.z);
                current_point=new_point=point;
                new_point*=mat;
                delta=new_point-current_point;
                pt[0]+=delta.x;pt[1]+=delta.y;pt[2]+=delta.z;
                if(save_points) pointCache->append(delta*(-1.0));
                else if(update_points && cache_index<cache_length){
                    MPoint& cache_point=(*pointCache)[cache_index];
                    cache_point-=delta;
                    cache_index++;}}}
        else{
            length=geometry->Get_Number_Vertices();
            if(save_points) pointCache->setSizeIncrement(length);
            for(unsigned int index=0;index<length;index++){
                double3& pt=builder.addElement(index).asDouble3();
                VECTOR_3D<float> phys_point=geometry->controls.array[index];
                MPoint point(phys_point.x,phys_point.y,phys_point.z);
                current_point=new_point=point;
                new_point*=mat;
                delta=new_point-current_point;
                pt[0]+=delta.x;pt[1]+=delta.y;pt[2]+=delta.z;
                if(save_points) pointCache->append(delta*-1.0);
                else if(update_points && index<cache_length) (*pointCache)[index]-=delta;}}}
    // assign the handle to be the new array
    handle.set(builder);
    //geometry->Tesselate();
    childChanged(MPxSurfaceShape::kBoundingBoxChanged);
}
//#####################################################################
// Function componentToPlugs
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE::
componentToPlugs(MObject& component,MSelectionList& list) const
{
    if(component.hasFn(MFn::kSingleIndexedComponent)){
        MFnSingleIndexedComponent vertex_component_fn(component);
        MPlug plug(thisMObject(),mControlPoints);
        convertToTweakNodePlug(plug); // if tweak node, plug should point to tweak instead
        int elements=vertex_component_fn.elementCount();
        for(int i=0;i<elements;i++){plug.selectAncestorLogicalIndex(vertex_component_fn.element(i),plug.attribute());list.add(plug);}}
}
//#####################################################################
// Function matchComponent
//#####################################################################
MPxSurfaceShape::MatchResult MAYA_SPLINE_PLUGIN_SHAPE::
matchComponent(const MSelectionList& item,const MAttributeSpecArray& spec,MSelectionList& list)
{
    MPxSurfaceShape::MatchResult result=MPxSurfaceShape::kMatchOk;
    MAttributeSpec attribute_spec=spec[0];
    int dim=attribute_spec.dimensions();
    // Look for attributes specifications of the form :  vtx[ index ] or vtx[ lower:upper ]
    if(spec.length()==1 && dim>0 && attribute_spec.name()=="vtx"){
        int vertex_count=Mesh_Geometry()->Get_Number_Vertices();
        MAttributeIndex attribute_index=attribute_spec[0];
        
        int upper=0,lower=0;
        if(attribute_index.hasLowerBound()) attribute_index.getLower(lower);
        if(attribute_index.hasUpperBound()) attribute_index.getUpper(upper);

        if(lower>upper || upper>=vertex_count) result=MPxSurfaceShape::kMatchInvalidAttributeRange;
        else{
            MDagPath path;item.getDagPath(0, path);
            MFnSingleIndexedComponent fn_vertex_component;
            MObject vertex_component=fn_vertex_component.create(MFn::kMeshVertComponent);
            for(int i=lower;i<=upper;i++) fn_vertex_component.addElement(i);
            list.add(path,vertex_component);}}
    else return MPxSurfaceShape::matchComponent(item,spec,list);
    return result;
}
//#####################################################################
// Function match
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
match(const MSelectionMask & mask,const MObjectArray& componentList) const
{
    bool result=false;
    if(componentList.length()==0) result=mask.intersects(MSelectionMask::kSelectMeshes);
    else for (int i=0;i<(int)componentList.length();i++) if(componentList[i].apiType()==MFn::kMeshVertComponent && mask.intersects(MSelectionMask::kSelectMeshVerts)){
        result=true;break;}
    return result;
}
//#####################################################################
// Function createFullVertexGroup
//#####################################################################
MObject MAYA_SPLINE_PLUGIN_SHAPE::
createFullVertexGroup() const
{
    MFnSingleIndexedComponent fn_component;
    MObject full_component=fn_component.create(MFn::kMeshVertComponent);        
    fn_component.setCompleteData(((MAYA_SPLINE_PLUGIN_SHAPE*)this)->Mesh_Geometry()->Get_Number_Vertices());
    return full_component;
}
//#####################################################################
// Function localShapeInAttr
//#####################################################################
MObject MAYA_SPLINE_PLUGIN_SHAPE::
localShapeInAttr() const
{
    return input_surface;
}
//#####################################################################
// Function localShapeOutAttr
//#####################################################################
MObject MAYA_SPLINE_PLUGIN_SHAPE::
localShapeOutAttr() const
{
    return output_surface;
}
//#####################################################################
// Function worldShapeOutAttr
//#####################################################################
MObject MAYA_SPLINE_PLUGIN_SHAPE::
worldShapeOutAttr() const
{
    return world_surface;
}
//#####################################################################
// Function cachedShapeAttr
//#####################################################################
MObject MAYA_SPLINE_PLUGIN_SHAPE::
cachedShapeAttr() const
{
    return cached_surface;
}
//#####################################################################
// Function geometryIteratorSetup
//#####################################################################
MPxGeometryIterator* MAYA_SPLINE_PLUGIN_SHAPE::
geometryIteratorSetup(MObjectArray& componentList,MObject& components,bool forReadOnly)
{
    if(components.isNull()) return new MAYA_SPLINE_PLUGIN_ITERATOR(Mesh_Geometry(),componentList);
    else return new MAYA_SPLINE_PLUGIN_ITERATOR(Mesh_Geometry(),components);
}
//#####################################################################
// Function acceptsGeometryIterator
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
acceptsGeometryIterator(bool writeable)
{
    return true;
}
//#####################################################################
// Function acceptsGeometryIterator
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
acceptsGeometryIterator(MObject&, bool writeable,bool forReadOnly)
{
    return true;
}
//#####################################################################
// Function Value
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
Value(int point_index,int component_index,double& value) const
{
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=((MAYA_SPLINE_PLUGIN_SHAPE*)this)->Cached_Geometry();
    if(geometry){VECTOR_3D<float> vector_value=geometry->controls.array[point_index];value=vector_value[component_index+1];return true;}
    return false;
}
//#####################################################################
// Function Value
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
Value(int point_index,MPoint& value) const
{
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=((MAYA_SPLINE_PLUGIN_SHAPE*)this)->Cached_Geometry();
    if(geometry){VECTOR_3D<float> point=geometry->controls.array[point_index];value=MPoint(point.x,point.y,point.z);return true;}
    return false;
}
//#####################################################################
// Function Set_Value
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
Set_Value(int point_index,int component_index,double value)
{
    bool result=false;
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=Cached_Geometry();
    if(geometry){geometry->controls.array[point_index][component_index+1]=value;result=true;}
    Vertices_Changed();
    return result;
}
//#####################################################################
// Function Set_Value
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE::
Set_Value(int point_index,const MPoint& value)
{
    bool result=false;
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=Cached_Geometry();
    if(geometry){geometry->controls.array[point_index]=VECTOR_3D<float>(value.x,value.y,value.z);result=true;}
    Vertices_Changed();
    return result;
}
//#####################################################################
// Function Cached_Geometry
//#####################################################################
MAYA_SPLINE_PLUGIN_GEOMETRY* MAYA_SPLINE_PLUGIN_SHAPE::
Cached_Geometry()
{
    MStatus stat;
    MAYA_SPLINE_PLUGIN_GEOMETRY* result=0;
    MDataBlock datablock=forceCache(); // bypass compute to get data
    MDataHandle handle=datablock.outputValue(cached_surface);
    MFnPluginData fn_data(handle.data());
    MAYA_SPLINE_PLUGIN_DATA* data=(MAYA_SPLINE_PLUGIN_DATA*)fn_data.data(&stat);
    if(stat != MS::kSuccess){std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE  Cached_Geometry: Failed to get MAYA_SPLINE_PLUGIN_SHAPE"<<std::endl;}
    if(data != NULL) result=data->geometry;
    return result;
}
//#####################################################################
// Function Mesh_Geometry
//#####################################################################
MAYA_SPLINE_PLUGIN_GEOMETRY* MAYA_SPLINE_PLUGIN_SHAPE::
Mesh_Geometry()
{
    MStatus stat;
    MAYA_SPLINE_PLUGIN_GEOMETRY* result=0;
    MDataBlock datablock=forceCache(); // bypass compute to get data
    MDataHandle handle=datablock.inputValue(output_surface);
    MFnPluginData fn_data(handle.data());
    MAYA_SPLINE_PLUGIN_DATA* data=(MAYA_SPLINE_PLUGIN_DATA*)fn_data.data(&stat);
    if(stat != MS::kSuccess){std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE  Cached_Geometry: Failed to get MAYA_SPLINE_PLUGIN_SHAPE"<<std::endl;}
    if(data != NULL) result=data->geometry;
    return result;
}
//#####################################################################
// Function Vertices_Changed
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE::
Vertices_Changed()
{
    childChanged(MPxSurfaceShape::kBoundingBoxChanged);
    childChanged(MPxSurfaceShape::kObjectChanged);
}
//#####################################################################
// Function Build_Control_Points
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
Build_Control_Points(MDataBlock& datablock,int count)
{
    MStatus stat;
    MArrayDataHandle cp_handle=datablock.outputArrayValue(mControlPoints,&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: build_control_points error getting cp_handle");

    MArrayDataBuilder old_builder=cp_handle.builder();
    if(count != (int)old_builder.elementCount()){
        MArrayDataBuilder builder(old_builder);
        MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: create builder");

        for(int vertex=0;vertex<count;vertex++) builder.addElement(vertex).asDouble3();
        cp_handle.set(builder);}
    cp_handle.setAllClean();
    return stat;
}
//#####################################################################
// Function creator
//#####################################################################
void* MAYA_SPLINE_PLUGIN_SHAPE::
creator()
{
    return new MAYA_SPLINE_PLUGIN_SHAPE;
}
//#####################################################################
// Function initialize
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_SHAPE::
initialize()
{
    MStatus stat;
    MFnTypedAttribute typed_attr;
  
    std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE: Initialize Attributes"<<std::endl;

    // INPUTS
    input_surface=typed_attr.create("inputSurface","is",MAYA_SPLINE_PLUGIN_DATA::id,MObject::kNullObj,&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: create inputSurface attribute");
    typed_attr.setStorable(false);
    ADD_ATTRIBUTE(input_surface);
    // real world matrix... will be linked to transform node... work around for bug in Maya
    MFnMatrixAttribute matrix_fn;
    MStatus matrix_stat;
    real_world_matrix=matrix_fn.create("realWorldMatrix","rwm",MFnMatrixAttribute::kDouble,&matrix_stat);
    MCHECKERROR(matrix_stat,"addAttribute readWorldMatrix");
    matrix_fn.setHidden(false);matrix_fn.setKeyable(false);
    matrix_stat=addAttribute(real_world_matrix);

    // OUTPUTS
    MAKE_NUMERIC_ATTR(bbox_min_corner,"bbox_min_corner","bbmin",MFnNumericData::k3Double,0,false,false,false);
    MAKE_NUMERIC_ATTR(bbox_max_corner,"bbox_max_corner","bbmax",MFnNumericData::k3Double,0,false,false,false);
    // local/world output surface attributes
    output_surface=typed_attr.create("outputSurface","os",MAYA_SPLINE_PLUGIN_DATA::id,MObject::kNullObj,&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: create outputSurface attribute")
    typed_attr.setWritable(false); // this might need to go next
    ADD_ATTRIBUTE(output_surface);
    world_surface=typed_attr.create("worldSurface","ws",MAYA_SPLINE_PLUGIN_DATA::id,MObject::kNullObj,&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_SHAPE: create worldSurface attribute");
    typed_attr.setWritable(false);
    ADD_ATTRIBUTE(world_surface);
    // Cached surface used for file IO
    cached_surface=typed_attr.create("cachedSurface","cs",MAYA_SPLINE_PLUGIN_DATA::id,MObject::kNullObj,&stat);
    MCHECKERROR(stat,"create cachedSurface attribute");
    typed_attr.setReadable(true);typed_attr.setWritable(true);typed_attr.setStorable(true);
    ADD_ATTRIBUTE(cached_surface);
    
    // dependencies, surfaces affect other surfaces and bounding boxes
    ATTRIBUTE_AFFECTS(input_surface,output_surface);
    ATTRIBUTE_AFFECTS(input_surface,world_surface);
    ATTRIBUTE_AFFECTS(output_surface,world_surface);
    ATTRIBUTE_AFFECTS(input_surface,bbox_min_corner);
    ATTRIBUTE_AFFECTS(input_surface,bbox_max_corner);    
    ATTRIBUTE_AFFECTS(cached_surface,output_surface);
    ATTRIBUTE_AFFECTS(cached_surface,world_surface);
    // dependencies, control points affect the surfaces.
    ATTRIBUTE_AFFECTS(mControlPoints,output_surface);
    ATTRIBUTE_AFFECTS(mControlValueX,output_surface);
    ATTRIBUTE_AFFECTS(mControlValueY,output_surface);
    ATTRIBUTE_AFFECTS(mControlValueZ,output_surface);
    ATTRIBUTE_AFFECTS(mControlPoints,cached_surface);
    ATTRIBUTE_AFFECTS(mControlValueX,cached_surface);
    ATTRIBUTE_AFFECTS(mControlValueY,cached_surface);
    ATTRIBUTE_AFFECTS(mControlValueZ,cached_surface);
    ATTRIBUTE_AFFECTS(mControlPoints,world_surface);
    ATTRIBUTE_AFFECTS(mControlValueX,world_surface);
    ATTRIBUTE_AFFECTS(mControlValueY,world_surface);
    ATTRIBUTE_AFFECTS(mControlValueZ,world_surface);
    ATTRIBUTE_AFFECTS(real_world_matrix,world_surface);

    return MS::kSuccess;
}
//#####################################################################
