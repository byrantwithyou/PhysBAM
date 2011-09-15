//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "MAYA_PLUGIN_MACROS.h"
#include "MAYA_SPLINE_PLUGIN_CREATOR.h"
#include "MAYA_SPLINE_PLUGIN_DATA.h"
using namespace PhysBAM;

#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnPluginData.h>
#include <maya/MFnTypedAttribute.h>

MObject MAYA_SPLINE_PLUGIN_CREATOR::size_u;
MObject MAYA_SPLINE_PLUGIN_CREATOR::size_v;
MObject MAYA_SPLINE_PLUGIN_CREATOR::size_w;
MObject MAYA_SPLINE_PLUGIN_CREATOR::zeroed;
MObject MAYA_SPLINE_PLUGIN_CREATOR::input_mesh;
MObject MAYA_SPLINE_PLUGIN_CREATOR::output_surface;
MTypeId MAYA_SPLINE_PLUGIN_CREATOR::id(0x80098);


//#####################################################################
// Function MAYA_SPLINE_PLUGIN_CREATOR
//#####################################################################
MAYA_SPLINE_PLUGIN_CREATOR::
MAYA_SPLINE_PLUGIN_CREATOR()
{}
//#####################################################################
// Function ~MAYA_SPLINE_PLUGIN_CREATOR
//#####################################################################
MAYA_SPLINE_PLUGIN_CREATOR::
~MAYA_SPLINE_PLUGIN_CREATOR()
{}
//#####################################################################
// Function compute
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_CREATOR::
compute(const MPlug& plug,MDataBlock& data)
{
    MStatus stat;
    if(true) std::cerr << "MAYA_SPLINE_PLUGIN_SHAPE::compute : plug "<<plug.info()<<std::endl;
    if(plug==output_surface){
        // create new data
        MFnPluginData fn_data_creator;
        MTypeId data_id(MAYA_SPLINE_PLUGIN_DATA::id);
        MObject new_data_object=fn_data_creator.create(data_id,&stat);
        if(stat != MS::kSuccess){std::cerr<<"MAYA_SPLINE_PLUGIN_CREATOR: compute error creating MAYA_SPLINE_PLUGIN_DATA"<<std::endl;return MS::kFailure;}
        // get geometry data
        MAYA_SPLINE_PLUGIN_DATA* new_data=(MAYA_SPLINE_PLUGIN_DATA*)fn_data_creator.data(&stat);
        if(stat != MS::kSuccess){std::cerr<<"MAYA_SPLINE_PLUGIN_CREATOR: error getting at MAYA_SPLINE_PLUGIN_DATA data"<<std::endl;return MS::kFailure;}
        MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=new_data->geometry;
        // get size information
        MDataHandle size_u_handle=data.inputValue(size_u),size_v_handle=data.inputValue(size_v),size_w_handle=data.inputValue(size_w);
        VECTOR_3D<int> size(size_u_handle.asInt(),size_v_handle.asInt(),size_w_handle.asInt());
        MDataHandle zeroed_handle=data.inputValue(zeroed);bool zeroed_value=zeroed_handle.asBool();
        geometry->Initialize(GRID_3D<float>(size.x+1,size.y+1,size.z+1,0,1,0,1,0,1),zeroed_value);
        
        MDataHandle out_handle=data.outputValue(output_surface);
        out_handle.set(new_data);

        std::cerr<<"COMPUTE --- "<<geometry->grid<<std::endl;
        data.setClean(plug);
        return MS::kSuccess;}        
    else return MS::kUnknownParameter;
}
//#####################################################################
// Function creator
//#####################################################################
void* MAYA_SPLINE_PLUGIN_CREATOR::
creator()
{
    return new MAYA_SPLINE_PLUGIN_CREATOR();
}
//#####################################################################
// Function initialize
//#####################################################################
MStatus MAYA_SPLINE_PLUGIN_CREATOR::
initialize()
{
    MStatus stat;
    // make size attribute
    MFnNumericAttribute size_attribute_fn;
    MFnTypedAttribute typed_attr_fn;
    // size
    MAKE_NUMERIC_ATTR(size_u,"size_u","szu",MFnNumericData::kInt,1,false,false,true);
    MAKE_NUMERIC_ATTR(size_v,"size_v","szv",MFnNumericData::kInt,1,false,false,true);
    MAKE_NUMERIC_ATTR(size_w,"size_w","szw",MFnNumericData::kInt,1,false,false,true);
    MAKE_NUMERIC_ATTR(zeroed,"zeroed","z",MFnNumericData::kBoolean,true,false,false,true);
    // input mesh
    MAKE_TYPED_ATTR(input_mesh,"inputMesh","im",MFnData::kMesh,NULL);
    // make output surface
    output_surface=typed_attr_fn.create("outputSurface","os",MAYA_SPLINE_PLUGIN_DATA::id,MObject::kNullObj,&stat);
    MCHECKERROR(stat,"MAYA_SPLINE_PLUGIN_CREATOR: create output_surface attribute");
    typed_attr_fn.setWritable(false);
    ADD_ATTRIBUTE(output_surface);

    // set dependencies
    ATTRIBUTE_AFFECTS(size_u,output_surface);
    ATTRIBUTE_AFFECTS(size_v,output_surface);
    ATTRIBUTE_AFFECTS(size_w,output_surface);
    ATTRIBUTE_AFFECTS(zeroed,output_surface);

    return MS::kSuccess;
}
