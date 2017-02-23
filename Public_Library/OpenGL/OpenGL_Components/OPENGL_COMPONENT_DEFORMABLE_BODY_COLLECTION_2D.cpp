//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Fracture/TRIANGLES_OF_MATERIAL.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_B_SPLINE_2D.h>
#include <OpenGL/OpenGL/OPENGL_BEZIER_SPLINE_2D.h>
#include <OpenGL/OpenGL/OPENGL_LEVELSET_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D
//#####################################################################
template<class T> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D(STREAM_TYPE stream_type,const std::string& prefix,const int start_frame)
    :OPENGL_COMPONENT<T>(stream_type,"Deformable Object List"),prefix(prefix),
    frame_loaded(-1),valid(false),display_mode(0),draw_velocities(false),velocity_scale(0.025),
    deformable_body_collection(*new DEFORMABLE_BODY_COLLECTION<TV>(0,0)),
    velocity_field(stream_type,deformable_body_collection.particles.V,deformable_body_collection.particles.X),color_map(OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map()),
    selected_segmented_curve(-1),selected_triangulated_area(-1),selected_triangles_of_material(-1),
    selected_bezier_spline(-1),selected_embedded_curve(-1),selected_free_particles(-1),
    selected_b_spline(-1),selected_index(-1),selected_object(0)
{
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("cycle_display_mode",{[this](){Cycle_Display_Mode();},"Cycle embedded display mode"});
    // check for per frame structures
    if(FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/deformable_object_structures",prefix.c_str(),start_frame)))
        invalidate_deformable_objects_selection_each_frame=true;
    else invalidate_deformable_objects_selection_each_frame=false;

    is_animation=true;
    Initialize();
}
//#####################################################################
// Function ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D
//#####################################################################
template<class T> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D()
{
    delete color_map;
    delete &deformable_body_collection;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Initialize()
{}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Reinitialize(bool force)
{
    if(!(draw && (force || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)))) return;
    static bool first_time=true;
    std::string filename=LOG::sprintf("%s/%d/deformable_object_structures",prefix.c_str(),frame);
    bool read_static_variables=!deformable_body_collection.structures.m;
    int static_frame=-1;
    if(FILE_UTILITIES::File_Exists(filename)){static_frame=frame;read_static_variables=true;}
    if(read_static_variables && first_time) LOG::cout << "Deformable bodies static variables will be read each frame" << std::endl;
    deformable_body_collection.Read(stream_type,prefix,prefix,frame,static_frame,read_static_variables,true); // Currently this will exit if any of the files don't exist... we should
                                                                                                                     // change it to not do that

    if(read_static_variables){
        int m=deformable_body_collection.structures.m;
        active_list.Clean_Memory();active_list.Resize(m);active_list.Fill(true);
        embedded_curve_objects.Delete_Pointers_And_Clean_Memory();embedded_curve_objects.Resize(m);
        segmented_curve_objects.Delete_Pointers_And_Clean_Memory();segmented_curve_objects.Resize(m);
        bezier_spline_objects.Delete_Pointers_And_Clean_Memory();bezier_spline_objects.Resize(m);
        b_spline_objects.Delete_Pointers_And_Clean_Memory();b_spline_objects.Resize(m);
        triangulated_area_objects.Delete_Pointers_And_Clean_Memory();triangulated_area_objects.Resize(m);
        triangles_of_material_objects.Delete_Pointers_And_Clean_Memory();triangles_of_material_objects.Resize(m);
        free_particles_objects.Delete_Pointers_And_Clean_Memory();free_particles_objects.Resize(m);
        free_particles_indirect_arrays.Delete_Pointers_And_Clean_Memory();free_particles_indirect_arrays.Resize(m);
        phi_list.Delete_Pointers_And_Clean_Memory();phi_list.Resize(m);
        int color_map_index=14;
        for(int i=0;i<m;i++){
            STRUCTURE<TV>* structure=deformable_body_collection.structures(i);
            if(EMBEDDED_MATERIAL_SURFACE<TV,2>* embedding=dynamic_cast<EMBEDDED_MATERIAL_SURFACE<TV,2>*>(structure)){
                has_embedded_objects=true;
                if(first_time) LOG::cout<<"object "<<i<<": embedded triangulated area\n";
                triangles_of_material_objects(i)=new OPENGL_TRIANGULATED_AREA<T>(stream_type,embedding->material_surface,false,OPENGL_COLOR::Red(),OPENGL_COLOR::Blue());
                embedding->embedded_object.simplicial_object.mesh.Initialize_Segment_Mesh();
                triangulated_area_objects(i)=new OPENGL_TRIANGULATED_AREA<T>(stream_type,embedding->embedded_object.simplicial_object,true);
                if(first_time) LOG::cout<<"object "<<i<<": embedded segmented curve\n";
                embedded_curve_objects(i)=new OPENGL_SEGMENTED_CURVE_2D<T>(stream_type,embedding->embedded_object.embedded_object);
                embedded_curve_objects(i)->draw_velocities=draw_velocities;
                embedded_curve_objects(i)->velocity_scale=velocity_scale;} // apply current parameters
            else if(SEGMENTED_CURVE_2D<T>* segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": segmented curve\n";
                segmented_curve_objects(i)=new OPENGL_SEGMENTED_CURVE_2D<T>(stream_type,*segmented_curve,color_map->Lookup(color_map_index--));
                segmented_curve_objects(i)->draw_velocities=draw_velocities;
                segmented_curve_objects(i)->velocity_scale=velocity_scale;} // apply current parameters
            else if(BEZIER_SPLINE<TV,3>* bezier_spline=dynamic_cast<BEZIER_SPLINE<TV,3>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": bezier spline\n";
                bezier_spline_objects(i)=new OPENGL_BEZIER_SPLINE_2D<T,3>(stream_type,*bezier_spline,color_map->Lookup(color_map_index--));
                bezier_spline_objects(i)->draw_velocities=draw_velocities;
                bezier_spline_objects(i)->velocity_scale=velocity_scale;} // apply current parameters
            else if(B_SPLINE<TV,3>* b_spline=dynamic_cast<B_SPLINE<TV,3>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": b-spline\n";
                b_spline_objects(i)=new OPENGL_B_SPLINE_2D<T,3>(stream_type,*b_spline,color_map->Lookup(color_map_index--));
                b_spline_objects(i)->draw_velocities=draw_velocities;
                b_spline_objects(i)->velocity_scale=velocity_scale;} // apply current parameters
            else if(TRIANGULATED_AREA<T>* triangulated_area=dynamic_cast<TRIANGULATED_AREA<T>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": triangulated area\n";
                triangulated_area->mesh.Initialize_Segment_Mesh(); // to enable segment selection
                triangulated_area_objects(i)=new OPENGL_TRIANGULATED_AREA<T>(stream_type,*triangulated_area,true);}
            else PHYSBAM_FATAL_ERROR(LOG::sprintf("Weird object %d",i));}}
    for(int i=0;i<deformable_body_collection.structures.m;i++){
        std::string suffix=LOG::sprintf("_%d",i);
        std::string frame_prefix=LOG::sprintf("%s/%d",prefix.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(frame_prefix+"stress_map_of_triangulated_area"+suffix)){
            if(first_time) LOG::cout<<"adding stress map to triangulated area"<<std::endl;
            ARRAY<OPENGL_COLOR > *color_map=new ARRAY<OPENGL_COLOR >;
            FILE_UTILITIES::Read_From_File(stream_type,frame_prefix+"stress_map_of_triangulated_area"+suffix,*color_map);
            triangulated_area_objects(i)->Set_Color_Map(color_map);}}
    frame_loaded=frame;
    valid=true;
    first_time=false;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(LOG::sprintf("%s%d/deformable_object_particles",prefix.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    if(draw_input) active_list.Fill(true);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Display() const
{
    if(!draw||!valid) return;
    bool draw_triangulated_areas=false,draw_triangles_of_material=false,draw_embedded_curves=false;
    switch(display_mode){
        case 0:draw_triangulated_areas=true;draw_embedded_curves=true;break;
        case 1:draw_triangles_of_material=true;break;
        case 2:draw_triangles_of_material=true;break;
        case 3:draw_triangulated_areas=true;break;
        case 4:draw_triangulated_areas=true;draw_triangles_of_material=true;break;}

    glPushName(0);
    glPushName(0);
    for(int i=0;i<segmented_curve_objects.m;i++){
        if(!active_list(i)) continue;
        if(segmented_curve_objects(i)){
            glLoadName(i);
            segmented_curve_objects(i)->Display();}}
    glPopName();
    glPopName();

    glPushName(1);
    glPushName(0);
    if(draw_triangulated_areas){
        for(int i=0;i<triangulated_area_objects.m;i++){
            if(!active_list(i)) continue;
            if(triangulated_area_objects(i)){
                glLoadName(i);
                triangulated_area_objects(i)->Display();}}}
    glPopName();
    glPopName();

    glPushName(2);
    glPushName(0);
    if(draw_triangles_of_material){
        for(int i=0;i<triangles_of_material_objects.m;i++){
            if(!active_list(i)) continue;
            if(triangles_of_material_objects(i)){
                glLoadName(i);
                triangles_of_material_objects(i)->Display();}}}
    glPopName();
    glPopName();
    
    glPushName(3);
    glPushName(0);
    for(int i=0;i<bezier_spline_objects.m;i++){
        if(!active_list(i)) continue;
        if(bezier_spline_objects(i)){
            glLoadName(i);
            bezier_spline_objects(i)->Display();}}
    glPopName();
    glPopName();

    glPushName(6);
    glPushName(0);
    for(int i=0;i<b_spline_objects.m;i++){
        if(!active_list(i)) continue;
        if(b_spline_objects(i)){
            glLoadName(i);
            b_spline_objects(i)->Display();}}
    glPopName();
    glPopName();

    glPushName(4);
    glPushName(0);
    if(draw_embedded_curves){
        for(int i=0;i<segmented_curve_objects.m;i++){
            if(!active_list(i)) continue;
            if(embedded_curve_objects(i)){
                glLoadName(i);
                embedded_curve_objects(i)->Display();}}}
    glPopName();
    glPopName();

    glPushName(5);
    glPushName(0);
    for(int i=0;i<segmented_curve_objects.m;i++){
        if(!active_list(i)) continue;
        if(free_particles_objects(i) && display_mode!=1){
            glLoadName(i);
            free_particles_objects(i)->Display();}}
    glPopName();
    glPopName();

    if(draw_velocities && deformable_body_collection.particles.store_velocity) velocity_field.Display();
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Increase_Vector_Size()
{
    T magnitude_adjustment=(T)1.1;
    velocity_scale*=magnitude_adjustment;
    for(int i=0;i<segmented_curve_objects.m;i++)
        if(active_list(i) && segmented_curve_objects(i))
            segmented_curve_objects(i)->velocity_scale=velocity_scale; // apply to objects which already exist
    for(int i=0;i<bezier_spline_objects.m;i++)
        if(active_list(i) && bezier_spline_objects(i))
            bezier_spline_objects(i)->velocity_scale=velocity_scale;
    for(int i=0;i<b_spline_objects.m;i++)
        if(active_list(i) && b_spline_objects(i))
            b_spline_objects(i)->velocity_scale=velocity_scale;
    for(int i=0;i<triangulated_area_objects.m;i++)
        if(active_list(i) && triangulated_area_objects(i))
            triangulated_area_objects(i)->velocity_scale=velocity_scale;
    velocity_field.Scale_Vector_Size(magnitude_adjustment);
    for(int i=0;i<embedded_curve_objects.m;i++)
        if(active_list(i) && embedded_curve_objects(i))
            embedded_curve_objects(i)->velocity_scale=velocity_scale; // apply to objects which already exist
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Decrease_Vector_Size()
{
    T magnitude_adjustment=(T)1/(T)1.1;
    velocity_scale*=magnitude_adjustment;
    for(int i=0;i<segmented_curve_objects.m;i++)
        if(active_list(i) && segmented_curve_objects(i))
            segmented_curve_objects(i)->velocity_scale=velocity_scale; // apply to objects which already exist
    for(int i=0;i<bezier_spline_objects.m;i++)
        if(active_list(i) && bezier_spline_objects(i))
            bezier_spline_objects(i)->velocity_scale=velocity_scale;
    for(int i=0;i<b_spline_objects.m;i++)
        if(active_list(i) && b_spline_objects(i))
            b_spline_objects(i)->velocity_scale=velocity_scale;
    for(int i=0;i<triangulated_area_objects.m;i++)
        if(active_list(i) && triangulated_area_objects(i))
            triangulated_area_objects(i)->velocity_scale=velocity_scale;
    velocity_field.Scale_Vector_Size(magnitude_adjustment);
    for(int i=0;i<embedded_curve_objects.m;i++)
        if(active_list(i) && embedded_curve_objects(i))
            embedded_curve_objects(i)->velocity_scale=velocity_scale; // apply to objects which already exist
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
    for(int i=0;i<segmented_curve_objects.m;i++)
        if(active_list(i) && segmented_curve_objects(i))
            segmented_curve_objects(i)->draw_velocities=draw_velocities;
    for(int i=0;i<bezier_spline_objects.m;i++)
        if(active_list(i) && bezier_spline_objects(i))
            bezier_spline_objects(i)->draw_velocities=draw_velocities;
    for(int i=0;i<b_spline_objects.m;i++)
        if(active_list(i) && b_spline_objects(i))
            b_spline_objects(i)->draw_velocities=draw_velocities;
    for(int i=0;i<triangulated_area_objects.m;i++)
        if(active_list(i) && triangulated_area_objects(i))
            triangulated_area_objects(i)->draw_velocities=draw_velocities;
    for(int i=0;i<embedded_curve_objects.m;i++)
        if(active_list(i) && embedded_curve_objects(i))
            embedded_curve_objects(i)->draw_velocities=draw_velocities;
}
//#####################################################################
// Function Cycle_Display_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Cycle_Display_Mode()
{
    display_mode=(display_mode+1)%5;
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Use_Bounding_Box() const
{
    return draw && valid && deformable_body_collection.structures.m>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Bounding_Box() const
{
    RANGE<VECTOR<T,3> > box;
    if(draw && valid && deformable_body_collection.structures.m>0){
        for(int i=0;i<segmented_curve_objects.m;i++)
            if(active_list(i) && segmented_curve_objects(i))
                box.Enlarge_To_Include_Box(segmented_curve_objects(i)->Bounding_Box());
        for(int i=0;i<bezier_spline_objects.m;i++)
            if(active_list(i) && bezier_spline_objects(i))
                box.Enlarge_To_Include_Box(bezier_spline_objects(i)->Bounding_Box());
        for(int i=0;i<b_spline_objects.m;i++)
            if(active_list(i) && b_spline_objects(i))
                box.Enlarge_To_Include_Box(b_spline_objects(i)->Bounding_Box());
        for(int i=0;i<triangulated_area_objects.m;i++)
            if(active_list(i) && triangulated_area_objects(i))
                box.Enlarge_To_Include_Box(triangulated_area_objects(i)->Bounding_Box());
        for(int i=0;i<triangles_of_material_objects.m;i++)
            if(active_list(i) && triangles_of_material_objects(i))
                box.Enlarge_To_Include_Box(triangles_of_material_objects(i)->Bounding_Box());}
    return box;
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    switch(indices(0)){
        case 0:return segmented_curve_objects(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
        case 1:return triangulated_area_objects(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
        case 2:return triangles_of_material_objects(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
        case 3:return bezier_spline_objects(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
        case 4:return embedded_curve_objects(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
        case 5:return free_particles_objects(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
        case 6:return b_spline_objects(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
        default:return 10;}
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(modifiers&GLUT_ACTIVE_CTRL){
        active_list(indices(1))=false;
        return false;}

    selected_index=indices(1);
    switch(indices(0)){
        case 0:
            selected_segmented_curve=indices(1);
            selected_object=segmented_curve_objects(indices(1));
            break;
        case 1:
            selected_triangulated_area=indices(1);
            selected_object=triangulated_area_objects(indices(1));
            break;
        case 2:
            selected_triangles_of_material=indices(1);
            selected_object=triangles_of_material_objects(indices(1));
            break;
        case 3:
            selected_bezier_spline=indices(1);
            selected_object=bezier_spline_objects(indices(1));
            break;
        case 4:
            selected_embedded_curve=indices(1);
            selected_object=embedded_curve_objects(indices(1));
            break;
        case 5:
            selected_free_particles=indices(1);
            selected_object=free_particles_objects(indices(1));
            break;
        case 6:
            selected_b_spline=indices(1);
            selected_object=b_spline_objects(indices(1));
            break;
        default:PHYSBAM_FATAL_ERROR();}
    return selected_object->Set_Selection(indices.Array_View(2,indices.m-2),modifiers);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Clear_Selection()
{
    for(int i=0;i<triangulated_area_objects.m;i++)if(triangulated_area_objects(i))triangulated_area_objects(i)->Clear_Selection();
    for(int i=0;i<triangles_of_material_objects.m;i++)if(triangles_of_material_objects(i))triangles_of_material_objects(i)->Clear_Selection();
    for(int i=0;i<free_particles_objects.m;i++)if(free_particles_objects(i))free_particles_objects(i)->Clear_Selection();
    for(int i=0;i<embedded_curve_objects.m;i++)if(embedded_curve_objects(i))embedded_curve_objects(i)->Clear_Selection();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Print_Selection_Info(std::ostream& output_stream) const
{
    output_stream<<"Deformable object "<<selected_index<<std::endl;
    selected_object->Print_Selection_Info(output_stream);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Destroy_Selection_After_Frame_Change()
{
    if(invalidate_deformable_objects_selection_each_frame){
        Clear_Selection();
        return true;}
    return false;
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>::
Selection_Bounding_Box() const
{
    return World_Space_Box(selected_object->Bounding_Box());
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<float>;
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<double>;
}
