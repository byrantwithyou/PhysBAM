#include "MAYA_SPLINE_PLUGIN_CREATE_COMMAND.h"
#include "MAYA_SPLINE_PLUGIN_CREATOR.h"
#include "MAYA_SPLINE_PLUGIN_DATA.h"
#include "MAYA_SPLINE_PLUGIN_EXPORT_COMMAND.h"
#include "MAYA_SPLINE_PLUGIN_IMPORT_COMMAND.h"
#include "MAYA_SPLINE_PLUGIN_SHAPE.h"
#include "MAYA_SPLINE_PLUGIN_SHAPE_UI.h"
#include <math.h>
#include <maya/MAttributeIndex.h>
#include <maya/MAttributeSpec.h>
#include <maya/MAttributeSpecArray.h>
#include <maya/MDagPath.h>
#include <maya/MFnAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnPlugin.h>
#include <maya/MFnPluginData.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MMatrix.h>
#include <maya/MObject.h>
#include <maya/MObjectArray.h>
//#####################################################################
// initializePlugin
//#####################################################################
MStatus initializePlugin(MObject obj)
{
    MFnPlugin plugin(obj,"Alias|Wavefront","4.5", "Any");
    // register the data node
    MStatus data_stat=plugin.registerData("MAYA_SPLINE_PLUGIN_DATA",PhysBAM::MAYA_SPLINE_PLUGIN_DATA::id,&PhysBAM::MAYA_SPLINE_PLUGIN_DATA::creator,MPxData::kGeometryData);
    if(!data_stat){std::cerr << "Failed to register geometry data : MAYA_SPLINE_PLUGIN_DATA \n"; return data_stat;}
    // register the creator node
    MStatus creator_stat=plugin.registerNode("MAYA_SPLINE_PLUGIN_CREATOR",PhysBAM::MAYA_SPLINE_PLUGIN_CREATOR::id,&PhysBAM::MAYA_SPLINE_PLUGIN_CREATOR::creator,&PhysBAM::MAYA_SPLINE_PLUGIN_CREATOR::initialize);
    if(!creator_stat){
        std::cerr<<"Failed to register creator : MAYA_SPLINE_PLUGIN_CREATOR"<<std::endl;
        plugin.deregisterData(PhysBAM::MAYA_SPLINE_PLUGIN_DATA::id);
        return creator_stat;}
    // register the shape node
    MStatus shape_stat=plugin.registerShape("MAYA_SPLINE_PLUGIN_SHAPE",PhysBAM::MAYA_SPLINE_PLUGIN_SHAPE::id,&PhysBAM::MAYA_SPLINE_PLUGIN_SHAPE::creator,&PhysBAM::MAYA_SPLINE_PLUGIN_SHAPE::initialize,&PhysBAM::MAYA_SPLINE_PLUGIN_SHAPE_UI::creator);
    if(!shape_stat){
        std::cerr<<"Failed to register shape : MAYA_SPLINE_PLUGIN_SHAPE"<<std::endl;
        plugin.deregisterData(PhysBAM::MAYA_SPLINE_PLUGIN_DATA::id);
        plugin.deregisterNode(PhysBAM::MAYA_SPLINE_PLUGIN_CREATOR::id);
        return shape_stat;}
    MStatus import_command_status=plugin.registerCommand("spline_import",PhysBAM::MAYA_SPLINE_PLUGIN_IMPORT_COMMAND::creator);
    MStatus export_command_status=plugin.registerCommand("spline_export",PhysBAM::MAYA_SPLINE_PLUGIN_EXPORT_COMMAND::creator);
    MStatus create_command_status=plugin.registerCommand("spline_create",PhysBAM::MAYA_SPLINE_PLUGIN_CREATE_COMMAND::creator);
    return MS::kSuccess;
}
//#####################################################################
// uninitializePlugin
//#####################################################################
MStatus uninitializePlugin(MObject obj)
{
    MFnPlugin plugin(obj);
    MStatus stat = plugin.deregisterData(PhysBAM::MAYA_SPLINE_PLUGIN_DATA::id);
    if(!stat){std::cerr<<"Failed to deregister data: MAYA_SPLINE_PLUGIN_DATA"<<std::endl;}
    stat=plugin.deregisterNode(PhysBAM::MAYA_SPLINE_PLUGIN_CREATOR::id);
    if(!stat){std::cerr<<"Failed to deregister creator : MAYA_SPLINE_PLUGIN_CREATOR"<<std::endl;}
    stat=plugin.deregisterNode(PhysBAM::MAYA_SPLINE_PLUGIN_SHAPE::id);
    if(!stat){std::cerr<<"Failed to deregister shape : MAYA_SPLINE_PLUGIN_SHAPE"<<std::endl;}
    plugin.deregisterCommand("spline_import");
    plugin.deregisterCommand("spline_export");
    plugin.deregisterCommand("spline_create");
    return stat;
}
