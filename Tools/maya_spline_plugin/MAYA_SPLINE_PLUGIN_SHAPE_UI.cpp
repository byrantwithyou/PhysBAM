//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include "MAYA_SPLINE_PLUGIN_GEOMETRY.h"
#include "MAYA_SPLINE_PLUGIN_SHAPE_UI.h"

#include <maya/MDagPath.h>
#include <maya/MDrawData.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MIOStream.h>
#include <maya/MMaterial.h>
#include <maya/MMatrix.h>
#include <maya/MObjectArray.h>
#include <maya/MSelectionList.h>
#include <maya/MSelectionMask.h>
#include <maya/MVector.h>

using namespace PhysBAM;

// Object and component color defines
#define LEAD_COLOR             (18) // green
#define ACTIVE_COLOR           (15) // white
#define ACTIVE_AFFECTED_COLOR  (8)  // purple
#define DORMANT_COLOR          (4)  // blue
#define HILITE_COLOR           (17) // pale blue
#define DORMANT_VERTEX_COLOR   (8)  // purple
#define ACTIVE_VERTEX_COLOR    (16) // yellow

#define POINT_SIZE (6.0)

//#####################################################################
// Function MAYA_SPLINE_PLUGIN_SHAPE_UI
//#####################################################################
MAYA_SPLINE_PLUGIN_SHAPE_UI::
MAYA_SPLINE_PLUGIN_SHAPE_UI()
{
}
//#####################################################################
// Function ~MAYA_SPLINE_PLUGIN_SHAPE_UI
//#####################################################################
MAYA_SPLINE_PLUGIN_SHAPE_UI::
~MAYA_SPLINE_PLUGIN_SHAPE_UI()
{
}
//#####################################################################
// Function getDrawRequests
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE_UI::
getDrawRequests(const MDrawInfo & info,bool objectAndActiveOnly,MDrawRequestQueue& queue)
{
    MAYA_SPLINE_PLUGIN_SHAPE* shape=(MAYA_SPLINE_PLUGIN_SHAPE*)surfaceShape();
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=shape->Mesh_Geometry();
    if(!geometry) std::cerr<<"NO Stuff to Draw"<<std::endl;

    MDrawRequest request=info.getPrototype(*this);
    MDrawData data;
    getDrawData(geometry,data);
    request.setDrawData(data);

    M3dView::DisplayStyle appearance=info.displayStyle();
    M3dView::DisplayStatus display_status=info.displayStatus();

    switch(appearance){
      case M3dView::kWireFrame:{
          request.setToken(kDrawWireframe);
          M3dView::ColorTable activeColorTable=M3dView::kActiveColors;
          M3dView::ColorTable dormantColorTable=M3dView::kDormantColors;
          switch(display_status){
            case M3dView::kLead: request.setColor(LEAD_COLOR,activeColorTable);break;
            case M3dView::kActive: request.setColor(ACTIVE_COLOR,activeColorTable);break;
            case M3dView::kActiveAffected: request.setColor(ACTIVE_AFFECTED_COLOR,activeColorTable);break;
            case M3dView::kDormant: request.setColor(DORMANT_COLOR,dormantColorTable);break;
            case M3dView::kHilite: request.setColor(HILITE_COLOR,activeColorTable);break;
            default: break;}
          queue.add(request);
          break;}
      case M3dView::kGouraudShaded:{
          request.setToken(kDrawSmoothShaded);
          MDagPath path=info.multiPath();
          M3dView view=info.view();
          MMaterial material=MPxSurfaceShapeUI::material(path);
          if(!material.evaluateMaterial(view,path)) std::cerr<<"MAYA_SPLINE_PLUGIN_SHAPE_UI: Couldnt evaluate material\n";
          request.setMaterial(material);
          bool material_transparent=false;
          material.getHasTransparency(material_transparent);
          if(material_transparent) request.setIsTransparent(true);
          queue.add(request);
          
          // if highlighted or selected, draw wireframe anyways...
          if(display_status==M3dView::kActive || display_status==M3dView::kLead || display_status==M3dView::kHilite){
              MDrawRequest wireRequest=info.getPrototype(*this);
              wireRequest.setDrawData(data);wireRequest.setToken(kDrawWireframeOnShaded);wireRequest.setDisplayStyle(M3dView::kWireFrame);
              M3dView::ColorTable active_color_table=M3dView::kActiveColors;

              switch(display_status){
                case M3dView::kLead: wireRequest.setColor(LEAD_COLOR,active_color_table);break;
                case M3dView::kActive: wireRequest.setColor(ACTIVE_COLOR,active_color_table);break;
                case M3dView::kHilite: wireRequest.setColor(HILITE_COLOR,active_color_table);break;
                default: break;}
              queue.add(wireRequest);}}
      case M3dView::kFlatShaded:
        request.setToken(kDrawFlatShaded);break;
      default:break;}

    if(!objectAndActiveOnly){
        if(appearance==M3dView::kPoints || display_status==M3dView::kHilite){
            MDrawRequest vertex_request=info.getPrototype(*this);
            vertex_request.setDrawData(data);vertex_request.setToken(kDrawVertices);vertex_request.setColor(DORMANT_VERTEX_COLOR,M3dView::kActiveColors);
            queue.add(vertex_request);}
        if(surfaceShape()->hasActiveComponents()){
            MDrawRequest active_vertex_request=info.getPrototype(*this);
            active_vertex_request.setDrawData(data);active_vertex_request.setToken(kDrawVertices);active_vertex_request.setColor(ACTIVE_VERTEX_COLOR,M3dView::kActiveColors);
            MObjectArray component_list=surfaceShape()->activeComponents();
            MObject vertex_component=component_list[0];
            active_vertex_request.setComponent(vertex_component);
            queue.add(active_vertex_request);}}
}
//#####################################################################
// Function draw
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE_UI::
draw(const MDrawRequest & request,M3dView & view) const
{
    int token=request.token();
    switch(token)
    {
      case kDrawWireframe:case kDrawWireframeOnShaded: Draw_Wireframe(request,view);break;
      case kDrawSmoothShaded:case kDrawFlatShaded: Draw_Shaded(request,view);break;;
      case kDrawVertices: Draw_Vertices(request,view);break;
    }
}
//#####################################################################
// Function select
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE_UI::
select(MSelectInfo &selectInfo, MSelectionList &selectionList,MPointArray &worldSpaceSelectPts) const
{
    bool selected=false;bool componentSelected=false;bool hilited=false;
    hilited=(selectInfo.displayStatus()==M3dView::kHilite);
    if (hilited){
        componentSelected=Select_Vertices(selectInfo,selectionList,worldSpaceSelectPts);
        selected=selected||componentSelected;}

    if(!selected){
        MPoint origin;MVector direction;
        selectInfo.getLocalRay(origin,direction);if(direction.x==1 && direction.y==0 && direction.z==0) return selected;
        RAY_3D<float> ray(VECTOR_3D<float>(origin.x,origin.y,origin.z),VECTOR_3D<float>(direction.x,direction.y,direction.z));
        MAYA_SPLINE_PLUGIN_SHAPE* shape=(MAYA_SPLINE_PLUGIN_SHAPE*)surfaceShape();
        TRIANGULATED_SURFACE<float>& surface=*shape->tetrahedralized_volume->triangulated_surface;
        for(int t=0;t<surface.triangle_mesh.triangles.m;t++){
            int p1,p2,p3;surface.triangle_mesh.triangles.Get(t,p1,p2,p3);
            TRIANGLE_3D<float> triangle(surface.particles.X(p1),surface.particles.X(p2),surface.particles.X(p3));
            if(triangle.Lazy_Intersection(ray)){
                std::cout<<"intersected "<<ray<<" at X="<<ray.Point(ray.t_max)<<std::endl;
                selected=true;break;}}
        if(selected){
            MSelectionMask priority_mask(MSelectionMask::kSelectNurbsSurfaces);
            MSelectionList item;
            item.add(selectInfo.selectPath());
            MPoint transformed_point;
            if(selectInfo.singleSelection()){
                transformed_point==shape->boundingBox().center();
                transformed_point*=selectInfo.selectPath().inclusiveMatrix();}
            selectInfo.addSelection(item,transformed_point,selectionList,worldSpaceSelectPts,priority_mask,false);}}

    return selected;
}
//#####################################################################
// Function Draw_Wireframe
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE_UI::
Draw_Wireframe(const MDrawRequest & request,M3dView & view) const
{
    MDrawData data = request.drawData();
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=(MAYA_SPLINE_PLUGIN_GEOMETRY*)data.geometry();
    MAYA_SPLINE_PLUGIN_SHAPE* shape=(MAYA_SPLINE_PLUGIN_SHAPE*)surfaceShape();

    view.beginGL(); 
        bool lighting_save=glIsEnabled(GL_LIGHTING)?true:false;
        if(lighting_save) glDisable(GL_LIGHTING);
        TRIANGULATED_SURFACE<float>& surface=*(shape->tetrahedralized_volume->triangulated_surface);
        glBegin(GL_LINES);
        for(int i=0;i<shape->uvw_iso_mesh.m;i++){
            int node1,node2;shape->uvw_iso_mesh.Get(i,node1,node2);
            OpenGL_Vertex(surface.particles.X(node1));
            OpenGL_Vertex(surface.particles.X(node2));}
        glEnd();
        if(lighting_save) glEnable(GL_LIGHTING);
    view.endGL();
}
//#####################################################################
// Function Draw_Shaded
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE_UI::
Draw_Shaded(const MDrawRequest & request,M3dView & view) const
{
    MDrawData data = request.drawData();
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=(MAYA_SPLINE_PLUGIN_GEOMETRY*)data.geometry();
    MAYA_SPLINE_PLUGIN_SHAPE* shape=(MAYA_SPLINE_PLUGIN_SHAPE*)surfaceShape();
    TRIANGULATED_SURFACE<float>& surface=*(shape->tetrahedralized_volume->triangulated_surface);
    const ARRAYS<int> &triangles=surface.triangle_mesh.triangles;
    view.beginGL();;
    // setup polygon offset
#if defined(SGI) || defined(MESA)
    glEnable( GL_POLYGON_OFFSET_EXT );
#else
    glEnable( GL_POLYGON_OFFSET_FILL );
#endif
    // material
    MMaterial material=request.material();
    material.setMaterial(request.multiPath(),request.isTransparent());
    // draw triangles
    glBegin(GL_TRIANGLES);
        for(int t=0;t<triangles.m;t++){
            int i,j,k;triangles.Get(t,i,j,k);
            OpenGL_Normal(TRIANGLE_3D<float>::Clockwise_Normal(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k)));
            OpenGL_Vertex(surface.particles.X(i));OpenGL_Vertex(surface.particles.X(j));OpenGL_Vertex(surface.particles.X(k));}
    glEnd();
    view.endGL();
}
//#####################################################################
// Function Draw_Vertices
//#####################################################################
void MAYA_SPLINE_PLUGIN_SHAPE_UI::
Draw_Vertices(const MDrawRequest & request,M3dView & view) const
{
    MDrawData data = request.drawData();
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=(MAYA_SPLINE_PLUGIN_GEOMETRY*)data.geometry();

    view.beginGL(); 
        // save lighting and point size and disable lighting
        bool lighting_save=glIsEnabled(GL_LIGHTING)?true:false;
        if(lighting_save) glDisable(GL_LIGHTING);
        float save_point_size;
        glGetFloatv(GL_POINT_SIZE,&save_point_size);
        glPointSize(POINT_SIZE);
    
        MObject component=request.component();
        if(!component.isNull()){
            MFnSingleIndexedComponent fn_component(component);
            for (int i=0;i<fn_component.elementCount();i++){
                int index=fn_component.element(i);
                glBegin(GL_POINTS);
                OpenGL_Vertex(geometry->controls.array[index]);
                glEnd();
                /*char annotation[32];
                sprintf(annotation,"%d",index);
                view.drawText(annotation,);}*/}}
        else{
            for(int i=0;i<geometry->Get_Number_Vertices();i++){
                glBegin( GL_POINTS );
                OpenGL_Vertex(geometry->controls.array[i]);
                glEnd();}}
    
        // Restore the state
        if(lighting_save) glEnable(GL_LIGHTING);
        glPointSize(save_point_size);
    view.endGL(); 
}
//#####################################################################
// Function Select_Vertices
//#####################################################################
bool MAYA_SPLINE_PLUGIN_SHAPE_UI::
Select_Vertices(MSelectInfo &selectInfo,MSelectionList &selectionList,MPointArray &worldSpaceSelectPts) const
{
    bool selected=false;
    M3dView view=selectInfo.view();
    MPoint transformed_point,current_point,selection_point;
    double z,previous_z=0;
    int closest_point_vertex_index=-1;
    const MDagPath& path=selectInfo.multiPath();

    // make component for selections
    MFnSingleIndexedComponent fn_component;
    MObject surface_component=fn_component.create(MFn::kMeshVertComponent);
    int vertex_index;

    // get alignment matrix
    MMatrix alignment_matrix;MPoint single_point;
    bool single_selection=selectInfo.singleSelection();
    if(single_selection) alignment_matrix=selectInfo.getAlignmentMatrix();
    
    // get geometry
    MAYA_SPLINE_PLUGIN_SHAPE* shape=(MAYA_SPLINE_PLUGIN_SHAPE*)surfaceShape();
    MAYA_SPLINE_PLUGIN_GEOMETRY* geometry=shape->Mesh_Geometry();

    // loop ofver each point to find if selection
    int vertex_count=geometry->Get_Number_Vertices();
    for(int vertex_index=0;vertex_index<vertex_count;vertex_index++){
        VECTOR_3D<float> phys_point=geometry->controls.array[vertex_index];
        MPoint current_point(phys_point.x,phys_point.y,phys_point.z);
        view.beginSelect();
        glBegin(GL_POINTS);
        glVertex3f((float)current_point[0],(float)current_point[1],(float)current_point[2]);
        glEnd();
        if(view.endSelect()>0){
            selected=true;
            if(single_selection){ // want single point so we need to find cvlosest
                transformed_point=current_point;transformed_point.homogenize();transformed_point*=alignment_matrix;z=transformed_point.z; // get depth
                if(closest_point_vertex_index<0 || z > previous_z){previous_z=z;single_point=current_point;closest_point_vertex_index=vertex_index;}}
            else fn_component.addElement(vertex_index);}}

     // get single point and transform world
    if(selected && single_selection){
        fn_component.addElement(closest_point_vertex_index);
        selection_point=single_point;
        selection_point*=path.inclusiveMatrix();}
    // add the selection component to the selection list
    if(selected){
        MSelectionList selection_item;
        selection_item.add(path,surface_component);
        MSelectionMask mask(MSelectionMask::kSelectNurbsSurfaces);
        selectInfo.addSelection(selection_item,selection_point,selectionList,worldSpaceSelectPts,mask,true);}

    return selected;
}
//#####################################################################
// Function creator
//#####################################################################
void* MAYA_SPLINE_PLUGIN_SHAPE_UI::
creator()
{
    return new MAYA_SPLINE_PLUGIN_SHAPE_UI;
}
//#####################################################################
