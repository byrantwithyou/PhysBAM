//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//######################################################################include <maya/MString.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include "MAYA_TRI_PLUGIN_EXTRACTOR.h"
#include <maya/MDagPath.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFnAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnTransform.h>
#include <maya/MGlobal.h>
#include <maya/MItDag.h>
#include <maya/MItSelectionList.h>
#include <maya/MObject.h>
#include <maya/MPlug.h>
#include <maya/MPoint.h>

using namespace PhysBAM;

template class MAYA_TRI_PLUGIN_EXTRACTOR<float>;

template<class T> MStatus MAYA_TRI_PLUGIN_EXTRACTOR<T>::
Extract_All()
{
    MGlobal::displayError("Export All is not supported yet");
    return MStatus::kFailure;
}

template<class T> MStatus MAYA_TRI_PLUGIN_EXTRACTOR<T>::
Extract_Selection()
{
    MStatus status;
    MSelectionList selection_list;
    if (MStatus::kFailure == MGlobal::getActiveSelectionList(selection_list)){MGlobal::displayError("MGlobal::getActiveSelectionList");return MStatus::kFailure;}
    MItSelectionList selection_list_iterator(selection_list, MFn::kMesh, &status);    
    if (MStatus::kFailure == status){MGlobal::displayError("Make selection iterator");return MStatus::kFailure;}
    for (selection_list_iterator.reset();!selection_list_iterator.isDone();selection_list_iterator.next()){
        MDagPath dag;
        if(MStatus::kFailure==selection_list_iterator.getDagPath(dag)){MGlobal::displayError("MItSelectionList::getDagPath");return MStatus::kFailure;}
        if (MStatus::kFailure==Extract_Mesh(dag)) return MStatus::kFailure;}
    return MStatus::kSuccess;
}

template<class T> MStatus MAYA_TRI_PLUGIN_EXTRACTOR<T>::
Extract_Mesh(const MDagPath dag_path)
{
    MFnMesh mesh(dag_path);
    int vertex_count=mesh.numVertices(),polygon_count=mesh.numPolygons();
    MFloatPointArray vertices;mesh.getPoints(vertices, MSpace::kWorld);
    MColorArray maya_colors;mesh.getVertexColors(maya_colors);
    // colors
    colors.Resize(colors.m+vertex_count);
    for(int i=0;i<vertex_count;i++){MColor color=maya_colors[i];colors(base+i)=VECTOR_3D<T>(color.r,color.g,color.b);}
    // particles
    particles.Increase_Array_Size(vertex_count);
    for(int i=0;i<vertex_count;i++){
        MFloatPoint point=vertices[i];
        int index=particles.Add_Particle();particles.X(index)=VECTOR_3D<T>(point.x,point.y,point.z);}
    // add connectivity
    triangle_mesh.number_nodes+=vertex_count;
    MIntArray connections;
    for(int i=0;i<polygon_count;i++){
        MIntArray poly_connects;
        if(mesh.polygonVertexCount(i) != 3){MGlobal::displayError("Not a triangulated surface...");return MStatus::kFailure;}
        mesh.getPolygonVertices(i,connections);
        triangle_mesh.triangles.Append(connections[0]+base,connections[1]+base,connections[2]+base);}
    // reset base
    base+=vertex_count;
    return MStatus::kSuccess;
}

template<class T> MStatus MAYA_TRI_PLUGIN_EXTRACTOR<T>::
Write(std::string filename_base)
{
    FILE_UTILITIES::Write_To_File<T>(filename_base+".tri",triangulated_surface);
    FILE_UTILITIES::Write_To_File<T>(filename_base+".col",colors);
    return MStatus::kSuccess;
}
