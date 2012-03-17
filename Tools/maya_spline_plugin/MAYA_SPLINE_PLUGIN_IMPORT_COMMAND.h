//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MAYA_SPLINE_PLUGIN_IMPORT_COMMAND__
#define __MAYA_SPLINE_PLUGIN_IMPORT_COMMAND__

#include <maya/MArgList.h>
#include <maya/MDagModifier.h>
#include <maya/MDGModifier.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnPlugin.h>
#include <maya/MPxCommand.h>
#include <maya/MString.h>

#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include "MAYA_PLUGIN_MACROS.h"
#include "MAYA_SPLINE_PLUGIN_CREATOR.h"
#include "MAYA_SPLINE_PLUGIN_SHAPE.h"

namespace PhysBAM{

class MAYA_SPLINE_PLUGIN_IMPORT_COMMAND:public MPxCommand
{
public:
    MAYA_SPLINE_PLUGIN_IMPORT_COMMAND()
    {}
    
    virtual ~MAYA_SPLINE_PLUGIN_IMPORT_COMMAND()
    {}

    MStatus doIt(const MArgList& arglist)
    {
        std::cout<<"in doit of import"<<std::endl;
        MStatus stat;
        if(arglist.length()!=1){
            std::cout<<"Usage is: spline_plugin_import <filename>"<<std::endl;
            return MS::kFailure;}

        // read the data
        std::string filename=arglist.asString(0,&stat).asChar();
        std::string dir_basename=FILE_UTILITIES::Get_Basename(filename);
        std::string::size_type lastslash=dir_basename.rfind("/");
        std::string basename=dir_basename.substr(lastslash+1);
        GRID_3D<float> grid;ARRAYS<VECTOR<VECTOR_3D<float> ,3> > controls;
        FILE_UTILITIES::Read_From_File<float>(filename,grid,controls);
        std::cout<<"read "<<filename<<" (basename="<<basename<<" and got grid "<<grid<<std::endl;
        // build creator node
        MDagModifier dagModifier;MDGModifier dgModifier;
        MObject creator_node=dgModifier.createNode(MAYA_SPLINE_PLUGIN_CREATOR::id,&stat);
        MCHECKERROR(stat,"Import failed, createNode MAYA_SPLINE_PLUGIN_CREATOR failed");
        dgModifier.renameNode(creator_node,MString(basename.c_str())+"_creator");
        dgModifier.doIt();
        // build shape node and transform (traverse transform to get shape)
        MObject shape_transform_node=dagModifier.createNode(MAYA_SPLINE_PLUGIN_SHAPE::id,MObject::kNullObj,&stat);
        MCHECKERROR(stat,"Import failed, createNode MAYA_SPLINE_PLUGIN_SHAPE failed");
        dagModifier.doIt();
        MFnDagNode shape_transform_fn(shape_transform_node);
        int children=shape_transform_fn.childCount();
        MCHECKERROR(stat,"Import failed, didn't see any children on shape transform");
        if(children<1){std::cerr<<"Import failed, didn't see any children";return MS::kFailure;}
        MObject shape_node=shape_transform_fn.child(0,&stat);
        dagModifier.renameNode(shape_node,MString(basename.c_str())+"_shape");
        dagModifier.renameNode(shape_transform_node,MString(basename.c_str())+"_transform");
        dagModifier.doIt();
        MCHECKERROR(stat,"Import failed Failed to get shape");
        // connect creator input surface to creator output surface
        MFnDependencyNode creator_fn(creator_node),shape_fn(shape_node);
        MPlug creator_output=creator_fn.findPlug("outputSurface");
        MPlug shape_input=shape_fn.findPlug("inputSurface");
        MDagModifier connectModifier;
        stat=connectModifier.connect(creator_output,shape_input);
        MCHECKERROR(stat,"Import failed, connect failed");
        stat=connectModifier.doIt();
        MCHECKERROR(stat,"Import failed, createNode failed");
        // set the grid size 
        MPlug grid_u_plug=creator_fn.findPlug("size_u"),grid_v_plug=creator_fn.findPlug("size_v"),grid_w_plug=creator_fn.findPlug("size_w");
        grid_u_plug.setValue(grid.m-1);grid_v_plug.setValue(grid.n-1);grid_w_plug.setValue(grid.mn-1);
        // set the control points
        MString shape_name=shape_fn.name(&stat);MCHECKERROR(stat,"Import failed, couldn't get name of shape node");
        for(int i=0;i<=grid.m+1;i++) for(int j=0;j<=grid.n+1;j++) for(int ij=0;ij<=grid.mn+1;ij++){
            int cp_index=controls.Standard_Index(i,j,ij);VECTOR_3D<float> point=controls(i,j,ij);
            std::cout<<"cp_index "<<cp_index<<" goes to "<<point<<std::endl;
            MGlobal::executeCommand(MString("setAttr ")+shape_name+".vtx["+cp_index+"] "+point.x+" "+point.y+" "+point.z);}
        // set the matrix
        MFnDagNode shape_dag_fn(shape_node);
        MPlug worldMatrix_output_compound=shape_transform_fn.findPlug("worldMatrix",&stat);
        MCHECKERROR(stat,"Import Failed to findPlug worldMatrix compound on transform");
        MPlug worldMatrix_output=worldMatrix_output_compound.elementByLogicalIndex(0);
        MPlug worldMatrix_input=shape_dag_fn.findPlug("realWorldMatrix",&stat);
        MCHECKERROR(stat,"Import Failed to findPlug realWorldMatrix on shape");
        MDagModifier connectTranModifier;
        stat=connectTranModifier.connect(worldMatrix_output,worldMatrix_input);
        MCHECKERROR(stat,"Import failed, connect transforms failed");
        stat=connectTranModifier.doIt();
        MCHECKERROR(stat,"Import failed, doIt on connect transforms failed");

        return MS::kSuccess;
    }
    
    static void* creator()
    {return new MAYA_SPLINE_PLUGIN_IMPORT_COMMAND;}

//#####################################################################
};
}
#endif
