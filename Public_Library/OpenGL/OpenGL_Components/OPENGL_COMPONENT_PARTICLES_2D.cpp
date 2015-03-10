//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
#include <sstream>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_PARTICLES_2D<T>::
OPENGL_COMPONENT_PARTICLES_2D(STREAM_TYPE stream_type,const std::string &filename_input, const std::string &filename_set_input, bool use_ids_input, bool particles_stored_per_cell_uniform_input)
    :OPENGL_COMPONENT<T>(stream_type,"Particles 2D"), particles(new GEOMETRY_PARTICLES<TV>),opengl_points(new OPENGL_POINTS_2D<T>(stream_type,*(new ARRAY<TV>))),
    opengl_vector_field(stream_type,*(new ARRAY<TV>), opengl_points->points, OPENGL_COLOR::Cyan()), 
      filename(filename_input), filename_set(filename_set_input), frame_loaded(-1), set(0), set_loaded(-1),number_of_sets(0),use_sets(false),valid(false),
      draw_velocities(false), have_velocities(false), use_ids(use_ids_input), 
      particles_stored_per_cell_uniform(particles_stored_per_cell_uniform_input),
      draw_multiple_particle_sets(false)
{
    viewer_callbacks.Set("toggle_draw_point_numbers",{[this](){Toggle_Draw_Point_Numbers();},"Toggle draw point numbers"});
    viewer_callbacks.Set("toggle_draw_radii",{[this](){Toggle_Draw_Radii();},"Toggle draw radii"});
    viewer_callbacks.Set("toggle_draw_velocities",{[this](){Toggle_Draw_Velocities();},"Toggle draw velocities"});
    viewer_callbacks.Set("command_prompt",{[this](){Command_Prompt();},"Command prompt"});
    viewer_callbacks.Set("next_set",{[this](){Next_Set();},"Switch to next set"});
    viewer_callbacks.Set("previous_set",{[this](){Previous_Set();},"Switch to previous set"});
    viewer_callbacks.Set("toggle_draw_multiple_particle_sets",{[this](){Toggle_Draw_Multiple_Particle_Sets();},"Toggle drawing multiple particle sets"});

    number_of_sets=0;
    while(filename_set!=""){
        std::string filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,number_of_sets);
        LOG::cout<<"Checking "<<filename<<std::endl;
        if(FILE_UTILITIES::File_Exists(filename)) number_of_sets++;
        else break;}
    if(number_of_sets>0){use_sets=true;draw_multiple_particle_sets=true;}
    else number_of_sets=1;

    particles_multiple.Resize(number_of_sets);opengl_points_multiple.Resize(number_of_sets);selected_ids.Resize(number_of_sets);
    particles_multiple(0)=particles;opengl_points_multiple(0)=opengl_points;
    OPENGL_INDEXED_COLOR_MAP* color_map=OPENGL_INDEXED_COLOR_MAP::Particle_Multiple_Color_Map();
    opengl_points->color=color_map->Lookup(0);
    for(int i=1;i<number_of_sets;i++){
        particles_multiple(i)=new GEOMETRY_PARTICLES<TV>;
        opengl_points_multiple(i)=new OPENGL_POINTS_2D<T>(stream_type,*(new ARRAY<TV>),color_map->Lookup(i));}
    delete color_map;
        
    is_animation=FILE_UTILITIES::Is_Animated(filename);
    // Don't need to call Reinitialize here because it will be called in first call to Set_Frame
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_PARTICLES_2D<T>::
~OPENGL_COMPONENT_PARTICLES_2D()
{
    for(int i=0;i<number_of_sets;i++){delete particles_multiple(i);delete &opengl_points_multiple(i)->points;delete opengl_points_multiple(i);}
    delete &opengl_vector_field.vector_field;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_PARTICLES_2D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Display() const
{
    if(valid && draw){
        GLint mode=0;
        glGetIntegerv(GL_RENDER_MODE, &mode);
        if(mode == GL_SELECT){
            if(draw_multiple_particle_sets){
                glPushName(0);
                for(int i=0;i<number_of_sets;i++){
                    glLoadName(i);
                    opengl_points_multiple(i)->Display();}
                glPopName();}
            else opengl_points->Display();
            if(draw_velocities && have_velocities) opengl_vector_field.Display();}
        else
        {
            if(draw_multiple_particle_sets) for(int i=0;i<number_of_sets;i++) opengl_points_multiple(i)->Display();
            else opengl_points->Display();
            if(draw_velocities && have_velocities) opengl_vector_field.Display();}}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_PARTICLES_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_points->Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_COMPONENT_PARTICLES_2D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    int particle_set;int point_index;
    if(buffer_size==1){particle_set=set;point_index=buffer[0];}
    else if(buffer_size==2){particle_set=buffer[0];point_index=buffer[1];}
    else return 0;

    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *selection=new OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>(this);

    // We have the OPENGL_POINTS_2D index but need to find the particle index
    int particle_index=0;
    int active_count=0;
    for(particle_index=0;particle_index<particles_multiple(particle_set)->Size();particle_index++){
        if(active_count++==point_index) break;}
    selection->index=particle_index;
    selection->particle_set=particle_set;
    selection->location=particles_multiple(particle_set)->X(particle_index);
    ARRAY_VIEW<int>* ids=0;if(use_ids) ids=Get_Particles_Id_Array(particle_set);
    if(ids){
        selection->has_id=true;
        selection->id=(*ids)(particle_index);}
    else selection->has_id=false;
    return selection;
}
//#####################################################################
// Function Get_Selection_By_Id
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_COMPONENT_PARTICLES_2D<T>::
Get_Selection_By_Id(int id,int particle_set)
{
    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *selection=0;
    ARRAY_VIEW<int>* ids=0;if(use_ids) ids=Get_Particles_Id_Array(particle_set);
    if(ids){
        int index;
        if(ids->Find(id,index)){
            selection=new OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>(this);
            selection->index=index;
            selection->has_id=true;
            selection->id=id;
            selection->particle_set=particle_set;
            selection->location=particles_multiple(particle_set)->X(index);}}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    if(selection->type != OPENGL_SELECTION<T>::COMPONENT_PARTICLES_2D) return;
    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;
    int particle_index=real_selection->index;
    ARRAY_VIEW<int>* ids=0;
    if(use_ids && (ids=Get_Particles_Id_Array(real_selection->particle_set))){
        Select_Particle_By_Id((*ids)(particle_index),real_selection->particle_set);}
    else{
        int point_index=0;
        for(int i=0;i<particle_index;i++)
            point_index++;
        opengl_points_multiple(real_selection->particle_set)->Select_Point(point_index);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Clear_Highlight()
{
    if(use_ids) Clear_Id_Selection();
    else for(int i=0;i<opengl_points_multiple.m;i++) opengl_points_multiple(i)->Clear_Selection();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const
{
    if(selection && selection->type == OPENGL_SELECTION<T>::COMPONENT_PARTICLES_2D && selection->object == this){
        OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;
        
        if(!draw_multiple_particle_sets && real_selection->particle_set!=set) return;
    
        output_stream<<"Selected particle in ["<<component_name<<"("<<real_selection->particle_set<<")] (total number = "<<particles->Size()<<")"<<std::endl;
    
        int current_index=-1;
        if(real_selection->has_id){
            output_stream<<"  Selected by id "<<real_selection->id<<std::endl;
            ARRAY_VIEW<int>* ids=Get_Particles_Id_Array(real_selection->particle_set);
            if(!ids || !ids->Find(real_selection->id,current_index))
                output_stream<<"  Doesn't exist"<<std::endl;}
        else{
            if(real_selection->index<particles_multiple(real_selection->particle_set)->Size()){
                output_stream<<"  Selected by index "<<real_selection->index<<std::endl;
                current_index=real_selection->index;}}

        if(current_index>=0){
            // real_selection->index is index into particles array at time of selection.  Not very useful. current_index is more useful
            output_stream<<"current index = "<<current_index<<std::endl;
            particles_multiple(real_selection->particle_set)->Print(output_stream,current_index);}}
}
//#####################################################################
// Function Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_COMPONENT_PARTICLES_2D<T>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION<T>* old_selection,bool& delete_selection)
{
    // TODO: reimplement transfering particles between objects.
    if(old_selection && old_selection->object == this){
        OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)old_selection;
        if(real_selection->has_id){
            OPENGL_SELECTION<T>* new_selection=Get_Selection_By_Id(real_selection->id,real_selection->particle_set);
            return new_selection;}
        else delete_selection=true;
    }
    return 0;
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_PARTICLES_2D<T>::
Selection_Bounding_Box(OPENGL_SELECTION<T>* selection) const
{
    int current_index=Get_Current_Index_Of_Selection(selection);
    if(current_index != -1) return World_Space_Box(RANGE<TV>(particles->X(current_index)));
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Get_Current_Index_Of_Selection
//#####################################################################
template<class T> int OPENGL_COMPONENT_PARTICLES_2D<T>::
Get_Current_Index_Of_Selection(OPENGL_SELECTION<T>* selection) const
{
    PHYSBAM_ASSERT(selection->type == OPENGL_SELECTION<T>::COMPONENT_PARTICLES_2D && selection->object == this);
    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;

    int current_index=-1;

    if(real_selection->has_id){
        ARRAY_VIEW<int>* ids=Get_Particles_Id_Array(real_selection->particle_set);
        if(ids) ids->Find(real_selection->id,current_index);}
    else if(real_selection->index<particles_multiple(real_selection->particle_set)->Size()) current_index=real_selection->index;

    return current_index;
}
//#####################################################################
// Function Get_Particle_Set_Of_Selection
//#####################################################################
template<class T> GEOMETRY_PARTICLES<VECTOR<T,2> >* OPENGL_COMPONENT_PARTICLES_2D<T>::
Get_Particle_Set_Of_Selection(OPENGL_SELECTION<T>* selection) const
{
    PHYSBAM_ASSERT(selection->type == OPENGL_SELECTION<T>::COMPONENT_PARTICLES_2D && selection->object == this);
    OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;
    PHYSBAM_ASSERT(real_selection->particle_set<particles_multiple.m);
    return particles_multiple(real_selection->particle_set);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Reinitialize(bool force)
{
    if(!draw || !(force || !valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))) return;
    valid=true;have_velocities=false;

    for(int i=0;i<number_of_sets;i++){
        std::string frame_filename;
        if(use_sets) frame_filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,i);
        else frame_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
        
        try{
            std::istream* input_file=FILE_UTILITIES::Safe_Open_Input(frame_filename);
            TYPED_ISTREAM typed_input(*input_file,stream_type);
            if(particles_stored_per_cell_uniform){
                ARRAY<GEOMETRY_PARTICLES<TV>*,VECTOR<int,2> > particles_per_cell;
                Read_Binary(typed_input,particles_per_cell);
                ARRAY<PARTICLES<TV>*> initialization_array(particles_per_cell.array.Size());
                for(int j=0;j<particles_per_cell.array.Size();j++){
                    if(particles_per_cell.array(j)) initialization_array(j)=particles_per_cell.array(j);
                    else initialization_array(j)=0;}
                ARRAY_VIEW<const PARTICLES<TV>* const> initialization_array_view(initialization_array.Size(),initialization_array.Get_Array_Pointer());
                particles_multiple(i)->Initialize(initialization_array_view);
                particles_per_cell.Delete_Pointers_And_Clean_Memory();}
            else{
                Read_Binary(typed_input,*particles_multiple(i));}
            delete input_file;
            opengl_points_multiple(i)->Set_Points_From_Particles(*particles_multiple(i),true,Get_Particles_Id_Array(i)!=0);}
        catch(FILESYSTEM_ERROR&){valid=false;}}
    frame_loaded=frame;
#if 0
    if(draw_velocities && particles->store_velocity){
        have_velocities=true;
        opengl_vector_field.vector_field.Resize(particles->Size());
        for(int i=0;i<particles->Size();i++)opengl_vector_field.vector_field(i)=particles->V(i);}
#endif
    if(use_ids) Apply_Id_Selection();
}
//#####################################################################
// Function Toggle_Draw_Point_Numbers
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Toggle_Draw_Point_Numbers()
{
    for(int i=0;i<opengl_points_multiple.m;i++)
        opengl_points_multiple(i)->draw_point_numbers=!opengl_points_multiple(i)->draw_point_numbers;
}
//#####################################################################
// Function Toggle_Draw_Radii
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Toggle_Draw_Radii()
{
    for(int i=0;i<opengl_points_multiple.m;i++)
        opengl_points_multiple(i)->draw_radii=!opengl_points_multiple(i)->draw_radii;
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Toggle_Draw_Velocities()
{
    draw_velocities=!draw_velocities;
    if(draw_velocities) Reinitialize(true);
}
//#####################################################################
// Function Command_Prompt_Response
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Command_Prompt_Response()
{
    if(!OPENGL_WORLD<T>::Singleton()->prompt_response.empty()){
        std::string command;
        std::istringstream sstream(OPENGL_WORLD<T>::Singleton()->prompt_response);
        sstream >> command;
        if(command == "s")
        {
            ARRAY<int> indices;
            int index;
            while(sstream >> index) indices.Append(index);
            if(use_ids && Get_Particles_Id_Array(set))
            {
                Clear_Id_Selection();
                Select_Particles_By_Ids(indices);
            }
            else 
            {
                opengl_points->Clear_Selection();
                opengl_points->Select_Points(indices);
            }
        }
    }
}
//#####################################################################
// Function Command_Prompt
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Command_Prompt()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Command: ",{[this](){Command_Prompt_Response();},""});
}
//#####################################################################
// Function Next_Set
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Next_Set()
{
    set=min(set+1,number_of_sets-1);
    LOG::cout<<"viewing particle set "<<set<<std::endl;
    particles=particles_multiple(set);opengl_points=opengl_points_multiple(set);
    Reinitialize();
}
//#####################################################################
// Function Previous_Set
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Previous_Set()
{
    set=max(set-1,0);
    LOG::cout<<"viewing particle set "<<set<<std::endl;
    particles=particles_multiple(set);opengl_points=opengl_points_multiple(set);
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Multiple_Particle_Sets
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Toggle_Draw_Multiple_Particle_Sets()
{
    draw_multiple_particle_sets=!draw_multiple_particle_sets;
}
//#####################################################################
// Function Select_Particle_By_Id
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Select_Particle_By_Id(int id, int particle_set)
{
    selected_ids(particle_set).Append(id);
    Apply_Id_Selection();
}
//#####################################################################
// Function Select_Particles_By_Ids
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Select_Particles_By_Ids(const ARRAY<int> &ids)
{
    // TODO: implement for multiple particle sets
    for(int i=0;i<ids.m;i++) selected_ids(0).Append(ids(i));
    Apply_Id_Selection();
}
//#####################################################################
// Function Clear_Id_Selection
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Clear_Id_Selection()
{
    for(int i=0;i<selected_ids.m;i++) selected_ids(i).Remove_All();
    Apply_Id_Selection();
}
//#####################################################################
// Function Apply_Id_Selection
//#####################################################################
template<class T> void OPENGL_COMPONENT_PARTICLES_2D<T>::
Apply_Id_Selection()
{
    for(int current_set=0;current_set<number_of_sets;current_set++){
        opengl_points_multiple(current_set)->Clear_Selection();
        if(use_ids && valid && opengl_points_multiple(current_set)->points.m == particles_multiple(current_set)->Size())
        {
            ARRAY_VIEW<int>* ids=Get_Particles_Id_Array(current_set);
            if(!ids) continue;
            int idx=0;
            for(int i=0;i<particles_multiple(current_set)->Size();i++){
                int dummy;
                if(selected_ids(current_set).Find((*ids)(i),dummy)) opengl_points_multiple(current_set)->Select_Point(idx);
                idx++;}
        }
    }
}
//#####################################################################
// Function Get_Particles_Id_Array
//#####################################################################
template<class T> ARRAY_VIEW<int>* OPENGL_COMPONENT_PARTICLES_2D<T>::
Get_Particles_Id_Array(int set_number) const
{
    if(set_number<0) set_number=set;
    ARRAY_VIEW<int>* ids=particles_multiple(set_number)->template Get_Array<int>(ATTRIBUTE_ID_ID);
    if(ids && ids->Size() && (*ids)(0)) return ids; // A hack to ignore ids if the first one equals zero
    return 0;
}
//#####################################################################
// Selection object functions
//#####################################################################
template<class T>
RANGE<VECTOR<T,3> > OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>::
Bounding_Box() const
{
    return object->Selection_Bounding_Box((OPENGL_SELECTION<T>*)this);
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_PARTICLES_2D<float>;
template class OPENGL_COMPONENT_PARTICLES_2D<double>;
}
