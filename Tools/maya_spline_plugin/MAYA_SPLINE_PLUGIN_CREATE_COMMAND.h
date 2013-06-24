//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MAYA_SPLINE_PLUGIN_CREATE_COMMAND__
#define __MAYA_SPLINE_PLUGIN_CREATE_COMMAND__

#include <maya/MArgList.h>
#include <maya/MDagModifier.h>
#include <maya/MDGModifier.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnPlugin.h>
#include <maya/MFnTransform.h>
#include <maya/MPxCommand.h>
#include <maya/MString.h>

#include <Tools/Read_Write/FILE_UTILITIES.h>
#include "MAYA_PLUGIN_MACROS.h"
#include "MAYA_SPLINE_PLUGIN_CREATOR.h"
#include "MAYA_SPLINE_PLUGIN_SHAPE.h"

namespace PhysBAM{

class MAYA_SPLINE_PLUGIN_CREATE_COMMAND:public MPxCommand
{
public:
    MAYA_SPLINE_PLUGIN_CREATE_COMMAND()
    {}
    
    virtual ~MAYA_SPLINE_PLUGIN_CREATE_COMMAND()
    {}

    MStatus doIt(const MArgList& arglist)
    {
        std::cout<<"in doit of create"<<std::endl;
        MStatus stat;

        // build creator node
        MDagModifier dagModifier;MDGModifier dgModifier;
        MObject creator_node=dgModifier.createNode(MAYA_SPLINE_PLUGIN_CREATOR::id,&stat);
        MCHECKERROR(stat,"Create failed, createNode MAYA_SPLINE_PLUGIN_CREATOR failed");
        dgModifier.doIt();
        // build shape node and transform (traverse transform to get shape)
        MObject shape_transform_node=dagModifier.createNode(MAYA_SPLINE_PLUGIN_SHAPE::id,MObject::kNullObj,&stat);
        MCHECKERROR(stat,"Create failed, createNode MAYA_SPLINE_PLUGIN_SHAPE failed");
        dagModifier.doIt();
        MFnDagNode shape_transform_fn(shape_transform_node);
        int children=shape_transform_fn.childCount();
        MCHECKERROR(stat,"Create failed, didn't see any children on shape transform");
        if(children<1){std::cerr<<"Create failed, didn't see any children";return MS::kFailure;}
        MObject shape_node=shape_transform_fn.child(0,&stat);
        //dagModifier.renameNode(shape_node,"SHAPE_file");
        dagModifier.doIt();
        MCHECKERROR(stat,"Create failed Failed to get shape");
        // connect creator input surface to creator output surface
        MFnDependencyNode creator_fn(creator_node),shape_fn(shape_node);
        MPlug creator_output=creator_fn.findPlug("outputSurface");
        MPlug shape_input=shape_fn.findPlug("inputSurface");
        MDagModifier connectModifier;
        stat=connectModifier.connect(creator_output,shape_input);
        MCHECKERROR(stat,"Create failed, connect failed");
        stat=connectModifier.doIt();
        MCHECKERROR(stat,"Create failed, createNode failed");
        // set the grid zeroed
        MPlug zeroed_plug=creator_fn.findPlug("zeroed");
        zeroed_plug.setValue(false);
        // set the matrix
        MFnDagNode shape_dag_fn(shape_node);
        MPlug worldMatrix_output_compound=shape_transform_fn.findPlug("worldMatrix",&stat);
        MCHECKERROR(stat,"Create Failed to findPlug worldMatrix compound on transform");
        MPlug worldMatrix_output=worldMatrix_output_compound.elementByLogicalIndex(0);
        MPlug worldMatrix_input=shape_dag_fn.findPlug("realWorldMatrix",&stat);
        MCHECKERROR(stat,"Create Failed to findPlug realWorldMatrix on shape");
        MDagModifier connectTranModifier;
        stat=connectTranModifier.connect(worldMatrix_output,worldMatrix_input);
        MCHECKERROR(stat,"Create failed, connect transforms failed");
        stat=connectTranModifier.doIt();
        MCHECKERROR(stat,"Create failed, doIt on connect transforms failed");

        return MS::kSuccess;
    }
    
    static void* creator()
    {return new MAYA_SPLINE_PLUGIN_CREATE_COMMAND;}

//#####################################################################
};
}
#endif
