//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include "MAYA_TRI_PLUGIN.h"
#include "MAYA_TRI_PLUGIN_EXTRACTOR.h"
#include <maya/MDagPath.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFnAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnPlugin.h>
#include <maya/MFnTransform.h>
#include <maya/MGlobal.h>
#include <maya/MItDag.h>
#include <maya/MItSelectionList.h>
#include <maya/MObject.h>
#include <maya/MPlug.h>
#include <maya/MPoint.h>
#include <maya/MString.h>

using namespace PhysBAM;

//#####################################################################
// Function creator 
//#####################################################################
template<class T> void* MAYA_TRI_PLUGIN<T>::
creator() 
{
    return new MAYA_TRI_PLUGIN();
}
//#####################################################################
// Function initializePlugin
//#####################################################################
MStatus
initializePlugin(MObject obj)
{
    MStatus status;
    MFnPlugin plugin(obj, "PhysBAM Triangulated Surface Importer/Exporter");//, "4.5", "Any");
    printf("register file\n");
    status=plugin.registerFileTranslator("PhysBAM Triangulated Surface",0,MAYA_TRI_PLUGIN<float>::creator,0,0,false);
    if (!status) status.perror("registerFileTranslator");
    return status;
}
//#####################################################################
// Function uninitializePlugin
//#####################################################################
MStatus
uninitializePlugin(MObject obj) 
{
    MStatus status;
    MFnPlugin plugin(obj);
    status=plugin.deregisterFileTranslator("PhysBAM Triangulated Surface");
    if (!status) status.perror("deregisterFileTranslator");
    return status;
}
//#####################################################################
// Function reader 
//#####################################################################
template<class T> MStatus MAYA_TRI_PLUGIN<T>::
reader (const MFileObject& file,const MString& optionsString,MPxFileTranslator::FileAccessMode mode) 
{
    const MString filename=file.fullName();std::string color_fileraw=filename.asChar();std::string color_filename=FILE_UTILITIES::Get_Basename(color_fileraw)+".col";
    //if(filename.substring(filename.length()-3,filename.length()-1) != "tri")
    if(!FILE_UTILITIES::File_Extension_Matches(filename.asChar(),"tri",false)){
        MGlobal::displayError(filename.substring(filename.length()-4,filename.length()-1) +": invalid extension not .tri");return MStatus::kFailure;}
    MGlobal::displayInfo("importexport::reader\n");
    MGlobal::displayInfo(color_filename.c_str());
    //std::ifstream input(filename.asChar(),std::ios::in|std::ios::binary);
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename.asChar(),true,false);
    if (!input){MGlobal::displayError(filename+": could not be opened for reading");return MS::kFailure;}
    // Read Triangulated surface
    TRIANGULATED_SURFACE<T> triangulated_surface(*(new TRIANGLE_MESH),*(new SOLIDS_PARTICLES<float,VECTOR_3D<float> >));
    triangulated_surface.template Read<T>(*input);
    delete input;
    int vertex_count=triangulated_surface.particles.number,triangle_count=triangulated_surface.triangle_mesh.triangles.m;
    // populate maya structures
    MGlobal::displayInfo(STRING_UTILITIES::string_sprintf("vertices: %d, triangles: %d",vertex_count,triangle_count).c_str());
    MStatus status=MStatus::kSuccess;
    MFloatPointArray maya_vertices(vertex_count);MIntArray maya_triangles(3*triangle_count);MIntArray poly_vertex_count(triangle_count,3);
    for(int i=0;i<vertex_count;i++){ // particles
        VECTOR_3D<T> point=triangulated_surface.particles.X(i+1);
        maya_vertices[i]=MFloatPoint(point.x,point.y,point.z);}
    for(int i=0;i<triangle_count;i++){ // connectivity
        maya_triangles[i*3+0]=triangulated_surface.triangle_mesh.triangles(1,i+1)-1;
        maya_triangles[i*3+1]=triangulated_surface.triangle_mesh.triangles(2,i+1)-1;
        maya_triangles[i*3+2]=triangulated_surface.triangle_mesh.triangles(3,i+1)-1;}
    // make maya object
    MFnMesh maya_mesh;MObject maya_transform=maya_mesh.create(vertex_count,triangle_count,maya_vertices,poly_vertex_count,maya_triangles,MObject::kNullObj,&status);
    // Read colors if they exist
    std::istream* color_input=FILE_UTILITIES::Safe_Open_Input(color_filename,true,false);
    bool have_colors=false;ARRAY<VECTOR_3D<float> > colors;
    if(color_input){
        MGlobal::displayInfo("found colors opened for reading");
        colors.template Read<T>(*color_input);
        MIntArray maya_vertices(vertex_count);MColorArray maya_colors(vertex_count,MColor(1,0,0));
        for(int i=0;i<vertex_count;i++){VECTOR_3D<T> color=colors(i+1);maya_vertices[i]=i;maya_colors[i]=MColor(color.x,color.y,color.z);}
        maya_mesh.setVertexColors(maya_colors,maya_vertices);
        delete color_input;}
    else{MGlobal::displayInfo("Didn't see color file...");}
    if(status!=MStatus::kSuccess){MGlobal::displayInfo("error creating mesh");return status;}
    MDGModifier dgModifier;
    dgModifier.renameNode(maya_transform,file.name());dgModifier.doIt();
    dgModifier.commandToExecute(MString("sets -e -fe initialShadingGroup ")+maya_mesh.name());
    MFnDagNode maya_dag_node(maya_transform,&status);
    if(status!=MStatus::kSuccess){MGlobal::displayInfo("error creating dag node");return status;}
    if(status==MStatus::kSuccess){dgModifier.commandToExecute("select "+maya_dag_node.name());dgModifier.doIt();}
    return MStatus::kSuccess;
}
//#####################################################################
// Function writer
//#####################################################################
template<class T> MStatus MAYA_TRI_PLUGIN<T>::
writer(const MFileObject& file,const MString& optionsString,MPxFileTranslator::FileAccessMode mode)
{
    MAYA_TRI_PLUGIN_EXTRACTOR<T> extractor;
    MStatus status;
    if(MPxFileTranslator::kExportAccessMode==mode) status=extractor.Extract_All();
    else if(MPxFileTranslator::kExportActiveAccessMode==mode) status=extractor.Extract_Selection();
    if(status != MStatus::kSuccess){MGlobal::displayError("Could not extract geometry");return status;}
    MString fname=file.fullName();MString fname_base=fname.substring(0,fname.length()-5);
    return extractor.Write(fname_base.asChar());
    MString suck;
    
}
//#####################################################################
// Function haveWriteMethod
//#####################################################################
template<class T> MString MAYA_TRI_PLUGIN<T>::
defaultExtension() const 
{
    return MString("tri");
}
//#####################################################################
// Function haveWriteMethod
//#####################################################################
template<class T> bool MAYA_TRI_PLUGIN<T>::
haveWriteMethod() const 
{
    return true;
}
//#####################################################################
// Function haveReadMethod
//#####################################################################
template<class T> bool MAYA_TRI_PLUGIN<T>::
haveReadMethod() const 
{
    return true;
}
//#####################################################################
// Function canBeOpened
//#####################################################################
template<class T> bool MAYA_TRI_PLUGIN<T>::
canBeOpened() const 
{
    return true;
}
//####################################################################
