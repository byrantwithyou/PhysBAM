//#####################################################################
// Copyright 2005, Andrew Selle, William Fong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VIEWER_WINDOW__
#define __VIEWER_WINDOW__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include "IMAGE_GL_WINDOW.h"
#include <Fl/Fl_Button.h>
#include <FL/Fl_Counter.h>
#include <Fl/Fl_Double_Window.h>
#include <FL/Fl_File_Chooser.H>
#include <Fl/Fl_Hold_Browser.h>
#include <Fl/Fl_Output.h>
#include <FL/Fl_Value_Input.h>
#include <Fl/Fl_Value_Output.h>
#include <FL/Fl_Window.h>

namespace PhysBAM{

template<class T>
class VIEWER_WINDOW:public Fl_Window
{
public:
    Fl_Counter counter;
    Fl_Value_Output current_index;
    Fl_Output current_label;
    Fl_Value_Input box_xmin_input,box_xmax_input;
    Fl_Value_Input box_ymin_input,box_ymax_input;
    Fl_Value_Input box_zmin_input,box_zmax_input;
    Fl_Hold_Browser tissue_browser;
    Fl_Button reset_boundaries_button;
    Fl_Value_Input granularity_x,granularity_y,granularity_z;
    Fl_Value_Input percent_granularity_input;
    Fl_Button add_button,delete_button;
    Fl_Button load_button,save_button;
    Fl_Button create_levelset_button;
    Fl_Button load_all,load_visible;
    VH_LEVELSET_BUILDER<T> visible_human;
    ARRAY<int> tissues;
    HASHTABLE<int,int> tissues_hash;
    VECTOR<int,2> current_point;
    RANGE<VECTOR<int,3> > box_size;
    std::string lastfile;

    IMAGE_GL_WINDOW<T> gl;
    VIEWER_WINDOW()
        :Fl_Window(0,0,1280,1024,"Visible Human Segmenter"),counter(5,10,200,20,"Slice"),
         current_index(40,70,200,20,"Index"),current_label(40,45,200,20,"Label"),
         box_xmin_input(75,220,40,20,"xmin:"),box_xmax_input(180,220,40,20,"xmax:"),
         box_ymin_input(75,250,40,20,"ymin:"),box_ymax_input(180,250,40,20,"ymax:"),
         box_zmin_input(75,280,40,20,"zmin:"),box_zmax_input(180,280,40,20,"zmax:"),
         tissue_browser(10,370,230,400,"Tissue"),reset_boundaries_button(40,195,195,20,"Reset Boundaries"),
         granularity_x(75,310,40,20,"granularity:"),granularity_y(130,310,40,20,NULL),granularity_z(180,310,40,20,NULL),
         percent_granularity_input(75,340,40,20,"percent:"),
         add_button(40,95,90,20,"Add"),delete_button(140,95,90,20,"Delete"),load_button(40,120,90,20,"Load"),
         save_button(140,120,90,20,"Save"),create_levelset_button(40,170,195,20,"Create Level Set"),
         load_all(40,145,90,20,"Load All"),load_visible(140,145,90,20,"Load Visible"),
         visible_human("slices"),
         gl(240,0,1040,1024,*this)
    {
        counter.bounds(1,visible_human.image_resolution.y);counter.callback(Slice_Change,this);counter.value(170);counter.step(1);counter.lstep(10);
        add_button.callback(Add_Tissue,this);delete_button.callback(Delete_Tissue,this);load_button.callback(Load_Tissue,this);
        save_button.callback(Save_Tissue,this);create_levelset_button.callback(Create_Levelset,this);reset_boundaries_button.callback(Reset_Box_Granularity,this);
        load_all.callback(Add_All_Tissues,this);load_visible.callback(Add_All_Visible_Tissues,this);percent_granularity_input.callback(Change_Gran_Percentage,this);
        box_xmin_input.callback(Update_Image_Granularity,this);box_xmax_input.callback(Update_Image_Granularity,this);box_ymin_input.callback(Update_Image_Granularity,this);
        box_ymax_input.callback(Update_Image_Granularity,this);box_zmin_input.callback(Update_Image_Granularity,this);box_zmax_input.callback(Update_Image_Granularity,this);
        Update_Image_Granularity(NULL,this);
        visible_human.Read_Slice((int)counter.value());
        percent_granularity_input.value(100);
        lastfile="levelset.phi";
        this->show();    
        gl.show();
    }

    void Update_Hover_Point(VECTOR<int,2> value)
    {current_point=value;current_index.value(visible_human.image(value));current_label.value(visible_human.labels(visible_human.image(value)).c_str());}

    void Update_Tissues()
    {tissues.Remove_All();tissues_hash.Delete_All_Entries();
    for(int i=1;i<=tissue_browser.size();i++){tissues.Append((int)(unsigned long long int)tissue_browser.data(i));tissues_hash.Insert((int)(unsigned long int)tissue_browser.data(i),0);}
    gl.Update_Image();}

private:
    static void Add_All_Tissues(Fl_Widget* widget,void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    window.visible_human.All_Tissues(&window.tissue_browser, window.tissues_hash);
    window.Update_Tissues();
    Reset_Box_Granularity(NULL,data);}

    static void Add_All_Visible_Tissues(Fl_Widget* widget,void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    HASHTABLE<int,int> *temp = new HASHTABLE<int,int>;
    for (int x= 1; x<= window.visible_human.image_resolution.x;x++){
        for (int y = 1; y<=window.visible_human.image_resolution.z;y++){
            int tissue_id=window.visible_human.image(VECTOR<int,2>(x,y));int dummy;
            if(!temp->Get(tissue_id,dummy) &&tissue_id != 0 && !window.tissues_hash.Get(tissue_id,dummy)
                         && window.visible_human.bounding_boxes(tissue_id).min_corner.x != -1){
                std::string s=str(boost::format("(%d) %s")%tissue_id%window.visible_human.labels(tissue_id).c_str());
                window.tissue_browser.add(s.c_str(),(void*)tissue_id);
                temp->Insert(tissue_id, 0);}}}
        delete temp;window.Update_Tissues();Reset_Box_Granularity(NULL,data);}
    
    static void Change_Gran_Percentage(Fl_Widget* widget,void* data)
    {Update_Image_Granularity(NULL,data);}

    static void Reset_Box_Granularity(Fl_Widget* widget,void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    Update_Bounding_Box(NULL,data);
    Update_Image_Granularity(NULL,data);
    window.percent_granularity_input.value(100);}

    static void Update_Bounding_Box(Fl_Widget* widget,void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    window.visible_human.original_box_size=window.visible_human.Get_Bounding_Box(window.tissues);
    window.box_xmin_input.value(window.visible_human.original_box_size.min_corner.x);window.box_ymin_input.value(window.visible_human.original_box_size.min_corner.y);window.box_zmin_input.value(window.visible_human.original_box_size.min_corner.z);
    window.box_xmax_input.value(window.visible_human.original_box_size.max_corner.x);window.box_ymax_input.value(window.visible_human.original_box_size.max_corner.y);window.box_zmax_input.value(window.visible_human.original_box_size.max_corner.z);}

    static void Update_Image_Granularity(Fl_Widget* widget,void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    RANGE<VECTOR<int,3> > image_box((int)window.box_xmin_input.value(),(int)window.box_xmax_input.value(),(int)window.box_ymin_input.value(),(int)window.box_ymax_input.value(),(int)window.box_zmin_input.value(),(int)window.box_zmax_input.value());
    VECTOR<int,3> grid_size=window.visible_human.Get_Grid_Size(image_box);
    double scale_factor=window.percent_granularity_input.value()/100;
    window.granularity_x.value((int)(grid_size.x*scale_factor));
    window.granularity_y.value((int)(grid_size.y*scale_factor));
    window.granularity_z.value((int)(grid_size.z*scale_factor));}

    static void Slice_Change(Fl_Widget* widget,void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    window.visible_human.Read_Slice((int)window.counter.value());
    window.gl.update_image=true;window.gl.damage(1);}

    static void Add_Tissue(Fl_Widget* widget,void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    int tissue_id=window.visible_human.image(window.current_point);
    int dummy;if(window.tissues_hash.Get(tissue_id,dummy)) return;
    std::string s=str(boost::format("(%d) %s")%tissue_id%window.visible_human.labels(tissue_id).c_str());
    window.tissue_browser.add(s.c_str(),(void*)tissue_id);
    window.Update_Tissues();
    Reset_Box_Granularity(NULL,data);}

    static void Delete_Tissue(Fl_Widget* widget,void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    if(window.tissue_browser.value() != -1) window.tissue_browser.remove(window.tissue_browser.value());
    window.Update_Tissues();
    Reset_Box_Granularity(NULL,data);}

    static void Load_Tissue(Fl_Widget* widget, void* data)
    {ARRAY<int> tissues_id_list;
    VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    char *filename=fl_file_chooser("Open Tissue List", "*.dat.gz","tissuelist.dat");
    if(!filename) return;
    FILE_UTILITIES::Read_From_File<T>(filename,tissues_id_list);
    int dummy;
    for(int i=0;i<tissues_id_list.m;i++){
        if(window.tissues_hash.Get(tissues_id_list(i),dummy)) continue;
        std::string s=str(boost::format("(%d) %s")%tissues_id_list(i)%window.visible_human.labels(tissues_id_list(i)).c_str());
        window.tissue_browser.add(s.c_str(),(void*)tissues_id_list(i));
        window.tissues_hash.Insert(tissues_id_list(i),0);}
    window.Update_Tissues();
    Reset_Box_Granularity(NULL,data);}

    static void Save_Tissue(Fl_Widget* widget, void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    char *filename=fl_file_chooser("Filename For Tissue List", "*.dat.gz","tissuelist.dat");
    if(filename) FILE_UTILITIES::Write_To_File<int>(filename,window.tissues);}

    static void Create_Levelset(Fl_Widget* widget,void* data)
    {VIEWER_WINDOW& window=*(VIEWER_WINDOW*)data;
    char *filename=fl_file_chooser("Filename For Levelset", "*.phi",window.lastfile.c_str());
    if (!filename) return;
    window.lastfile=filename;
    RANGE<VECTOR<int,3> > box((int)window.box_xmin_input.value(),(int)window.box_xmax_input.value(),(int)window.box_ymin_input.value(),(int)window.box_ymax_input.value(),(int)window.box_zmin_input.value(),(int)window.box_zmax_input.value());
    VECTOR<int,3> granularity((int)window.granularity_x.value(),(int)window.granularity_y.value(),(int)window.granularity_z.value());
    if(filename) window.visible_human.Create_Levelset(window.tissues,box,granularity,filename);}

//#####################################################################
};
}
#endif
