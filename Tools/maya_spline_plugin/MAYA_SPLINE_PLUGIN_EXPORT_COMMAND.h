//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MAYA_SPLINE_PLUGIN_EXPORT_COMMAND__
#define __MAYA_SPLINE_PLUGIN_EXPORT_COMMAND__

#include <maya/MArgList.h>
#include <maya/MDagModifier.h>
#include <maya/MDGContext.h>
#include <maya/MDGModifier.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnPlugin.h>
#include <maya/MItSelectionList.h>
#include <maya/MPxCommand.h>
#include <maya/MString.h>
#include <maya/MTypeId.h>

#include <Tools/Read_Write/FILE_UTILITIES.h>
#include "MAYA_PLUGIN_MACROS.h"
#include "MAYA_SPLINE_PLUGIN_CREATOR.h"
#include "MAYA_SPLINE_PLUGIN_SHAPE.h"

namespace PhysBAM{

class MAYA_SPLINE_PLUGIN_EXPORT_COMMAND:public MPxCommand
{
public:
    MAYA_SPLINE_PLUGIN_EXPORT_COMMAND()
    {}
    
    virtual ~MAYA_SPLINE_PLUGIN_EXPORT_COMMAND()
    {}

    MStatus doIt(const MArgList& arglist)
    {
        MStatus stat;
        if(arglist.length()!=1){
            std::cout<<"Usage is: spline_plugin_emport <filename>"<<std::endl;
            return MS::kFailure;}

        MSelectionList selection_list;
        stat=MGlobal::getActiveSelectionList(selection_list);
        MCHECKERROR(stat,"Import failed to get selection list");
        MTypeId temp_id(MAYA_SPLINE_PLUGIN_SHAPE::id);
        MItSelectionList selection_list_iterator(selection_list,MFn::kPluginShape);
        for (selection_list_iterator.reset();!selection_list_iterator.isDone();selection_list_iterator.next()){
            MDagPath dag;
            if(MStatus::kFailure==selection_list_iterator.getDagPath(dag)){MGlobal::displayError("MItSelectionList::getDagPath");return MStatus::kFailure;}
            MGlobal::displayInfo("Traversing "+dag.fullPathName()+" to find shapes");
            MObject node_object=dag.node(&stat);
            MCHECKERROR(stat,"  Getting MObject for Dag path");
            MFnDependencyNode shape_fn(node_object,&stat);
            MCHECKERROR(stat,"  Getting MFnDagDependencyNode");
            MPlug plug=shape_fn.findPlug("worldSurface",&stat);
            MCHECKERROR(stat,"  Failed to get worldSurface plug");
            MDGContext normal;MObject data_object;stat=plug.getValue(data_object,normal);
            MCHECKERROR(stat,"  Failed to get data_object MObject");
            MFnPluginData plugin_data(data_object);
            const MAYA_SPLINE_PLUGIN_DATA* data=(MAYA_SPLINE_PLUGIN_DATA*)plugin_data.constData(&stat);
            MCHECKERROR(stat,"  Failed to get mpxdata type turned into shape data");
            MGlobal::displayInfo(" Writing to "+dag.partialPathName());
            if(data){
                MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=data->geometry;
                if(geometry) FILE_UTILITIES::Write_To_File<float>(dag.partialPathName().asChar(),geometry->grid,geometry->controls);
                else MCHECKERROR(MS::kFailure,"Failed to get geometry pointer");}
            else MCHECKERROR(MS::kFailure,"Failed to get data pointer");}
        return MS::kSuccess;
    }
    
    static void* creator()
    {return new MAYA_SPLINE_PLUGIN_EXPORT_COMMAND;}

//#####################################################################
};
}
#endif
