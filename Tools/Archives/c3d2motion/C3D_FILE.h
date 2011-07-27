//##########################################################################################
// Class C3D_FILE
// The parser and holder for the data
//##########################################################################################
#ifndef _C3D_FILE_h
#define _C3D_FILE_h

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <cassert>
#include <map>
#include <string>
#include "C3D_PARAMETER_GROUP.h"
#include "C3D_PARSE_ERROR.h"
#include "C3D_POINT.h"

class C3D_FILE{
    PhysBAM::ARRAY<C3D_PARAMETER_GROUP> groups;//There can only be 128 groups,and they can have any number of parameters in them
    std::map<std::string,int> group_name_to_id;
    static const int MAX_GROUPS=128;
public:
    std::string error_string;//stores the parse error that was received
    PhysBAM::ARRAYS<VECTOR<C3D_POINT,2> > frames;//frames(i in [1,number_frames],j in [1,number_channels])
    C3D_BYTE original_machine_code;
    
    C3D_FILE()
    {
        groups.Preallocate(MAX_GROUPS);
    }
    
    void Initialize(const char *filename) //makes the c3d structure and parses the c3d file
        try
    {Read(filename);}
    catch (C3D_PARSE_ERROR error){error_string=error.error_string;}
    
    C3D_PARAMETER_GROUP* Get_Group(const char *name)
    {int i=group_name_to_id[name];return (1<=i&&i<=groups.m)?&groups(i):0;}
    
    // Get basic primitive data
    int Number_Frames(){return Get_Group("POINT")->Get("FRAMES")->Get_Int();}
    int Number_Channels(){return Get_Group("POINT")->Get("USED")->Get_Int();}
    float Frame_Rate(){return Get_Group("POINT")->Get("RATE")->Get_Float();}
    
    // For getting the channel labels and description
    void Get_Labels(PhysBAM::ARRAY<std::string>& label_array){return Get_Group("POINT")->Get("LABELS")->Fill_String_Array(label_array);}
    
    void Get_Descriptions(PhysBAM::ARRAY<std::string>& description_array)
    {C3D_PARAMETER_VALUE* value=Get_Group("POINT")->Get("DESCRIPTIONS");
    if(value) value->Fill_String_Array(description_array);else description_array.Resize(0);}
    
    void Display_Data()
    {printf("Parameter Data: %d groups\n",groups.m);printf("--------------\n");
    for(int i=1;i<=groups.m;++i){printf("Group 0x%02x \"%s\"\n",i,groups(i).Get_Name());groups(i).Display_Data();printf("\n");}}
    
    bool Valid_Frame_Index(const int i) const
    {return frames.m_start<=i && i<=frames.m_end;}
    
    bool Valid_Channel_Index(const int j) const
    {return frames.n_start<=j && j<=frames.n_end;}
    
    void Get_Frame_Data(const int frame,PhysBAM::ARRAY<C3D_POINT>& data_array)
    {int channels=Number_Channels();data_array.Resize(channels);for(int j=1;j<=channels;++j) data_array(j)=frames(frame,j);}
    
    int Number_Visible_Markers(const int frame) const
    {if(!Valid_Frame_Index(frame)) return 0;
    int n=0;for(int j=frames.n_start;j<=frames.n_end;++j) if(frames(frame,j).camera_visible) ++n;return n;}
    
    float Max_Residual(const int frame) const
    {if(!Valid_Frame_Index(frame)) return 0;
    float r=0;for(int j=frames.n_start;j<=frames.n_end;++j) r=PhysBAM::max(r,frames(frame,j).residual);return r;}
    
//##########################################################################################
private:
    bool Read(const char *filename);// does the real parsing job
    void Read_Main_Header(FILE *fp);
    void Read_Parameter_Section(FILE *fp);
    void Read_Point_Section(FILE *fp);
//##########################################################################################
};


#endif
