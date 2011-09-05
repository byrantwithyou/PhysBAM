//##########################################################################################
// File: C3D_FILE.cpp
//##########################################################################################
#include "C3D_FILE.h"
#include "C3D_PARSE_ERROR.h"
#include <math.h>

typedef unsigned short C3D_WORD;

// Reads a string used by description and name in the parameter section
static char * Read_String(int string_length,FILE *fp)
{
    char *name_buffer=new char[string_length+1];
    fread(name_buffer,1,string_length,fp);name_buffer[string_length] = '\0';return name_buffer;
}

//##########################################################################################
// Function  Read
// actually does the parsing
//##########################################################################################
bool C3D_FILE::
Read(const char *filename)
{
    FILE *fp=fopen(filename,"rb");if(!fp) return false;
    try{
        Read_Main_Header(fp);
        Read_Parameter_Section(fp);
        Read_Point_Section(fp);
        // finish with success
        fclose(fp);return true;}
    catch(C3D_PARSE_ERROR error){
        // die with failure
        fclose(fp);
        throw C3D_PARSE_ERROR(error.error_string);
        return false;}
}
//##########################################################################################
// Function  Read_Main_Header
// reads the main data header and verifies that it is a c3d
//##########################################################################################
void C3D_FILE::
Read_Main_Header(FILE *fp)
{  
    C3D_BYTE c3d_header[2];fread(c3d_header,sizeof(c3d_header[0]),2,fp);// Read first two bytes
    if (c3d_header[1]!=0x50) throw C3D_PARSE_ERROR("Invalid header flag (expected byte 2 to be 0x50)");// Sanity check
    // Find and skip to parameter section
    unsigned int parameter_section_offset=512*(c3d_header[0]-1);fseek(fp,parameter_section_offset,SEEK_SET);
}
//##########################################################################################
// Function  Read_Parameter_Section
// reads the parameter section
//##########################################################################################
void C3D_FILE::
Read_Parameter_Section(FILE *fp)
{  
    // Read the parameter header
    C3D_BYTE parameter_header[4];
    fread(parameter_header,sizeof(parameter_header[0]),4,fp);
    // first two bytes are crap, third is number of 512byte blocks, fourth is type
    //printf("Size of Parameter section %d\n",parameter_header[2]);
    //printf("Machine type 0x%0.2x (%d) \n",parameter_header[3],parameter_header[3]);
    original_machine_code = parameter_header[3];
    
    int count=0;
    while(1){ 
        ++count;//printf("%d\n",count);
        C3D_SIGNED_BYTE record_head[2];fread(record_head,sizeof(record_head[0]),2,fp);// read the group/param header to figure out what we got
        int name_length=abs(record_head[0]);if(name_length<=0) break;// check if this is the end record
        int group_id=1+abs(record_head[1]);// Figure out what group # we are
        bool is_group_definition=record_head[1]<0;// Figure out whether we are a param entry or group entry
        char *name_buffer=Read_String(name_length,fp);// Get string & terminate it
        C3D_WORD offset;fread(&offset,sizeof(offset),1,fp);// read offset to next
        
        groups.Resize(PhysBAM::max(groups.m,group_id));C3D_PARAMETER_ENTRY* read_entry;
        if (is_group_definition){// Make the group
            read_entry=&groups(group_id);group_name_to_id[name_buffer]=group_id;}
        else{//make a member for the group
            C3D_PARAMETER_VALUE* value=new C3D_PARAMETER_VALUE;value->Read_Data_Part(fp,original_machine_code);
            groups(group_id).Add_Member(name_buffer,value);
            read_entry=value;}
        // Read the description
        C3D_BYTE description_length;
        fread(&description_length,sizeof(description_length),1,fp);
        char *description_buffer=Read_String(description_length,fp);
        // set the basic info
        read_entry->Set_Info(name_buffer,description_buffer);
        delete[] name_buffer;delete[] description_buffer;
    }
}

//#######################
// point reading infrastructure
//#######################
//This is for the int structure
typedef struct _point_int_format{
    C3D_INT x,y,z;
    C3D_BYTE camera_mask,residual;
}point_int_format;

//This is for the float structure
typedef struct _point_float_format{
    C3D_FLOAT x,y,z;
    C3D_BYTE camera_mask,residual;
}point_float_format;

//##########################################################################################
// Function  Read_Point_Section
//##########################################################################################
void C3D_FILE::
Read_Point_Section(FILE *fp)
{  
    // first find the offset to the point section
    C3D_PARAMETER_GROUP& points_group=groups(group_name_to_id["POINT"]);
    int block_offset=points_group.Get("DATA_START")->Get_Int();
    int true_offset=(block_offset-1)*512;
    int number_frames=points_group.Get("FRAMES")->Get_Int();
    int number_channels=points_group.Get("USED")->Get_Int();
    float scale=points_group.Get("SCALE")->Get_Float();
    
    frames.Resize(1,1,number_frames,1,number_channels,true,false);//Allocate the data
    fseek(fp,true_offset,SEEK_SET);// Go to the point data section
    
    // Put the data in our structure (convert everything to floats,while we're at it)
    if(scale<0){//It is in float format
        point_float_format *temp_pts=new point_float_format[number_channels];
        for(int i=1;i<=number_frames;++i){//Read a frame
            fread(temp_pts,sizeof(temp_pts[0]),number_channels,fp);
            for(int j=1;j<=number_channels;++j){//copy a channel
                const point_float_format& temp=temp_pts[j-1];C3D_POINT& frame=frames(i,j);
                frame.x=temp.x;
                frame.y=temp.y;
                frame.z=temp.z;
                frame.camera_visible=temp.camera_mask;
                frame.residual=(float)temp.residual*scale;}}
        delete[] temp_pts;}
    else{//It is in int format
        point_int_format *temp_pts=new point_int_format[number_channels];
        for(int i=1;i<=number_frames;++i){//Read a frame
            fread(temp_pts,sizeof(temp_pts[0]),number_channels,fp);
            for(int j=1;j<=number_channels;++j){//copy a channel
                const point_int_format& temp=temp_pts[j-1];C3D_POINT& frame=frames(i,j);
                frame.x=temp.x*scale;
                frame.y=temp.y*scale;
                frame.z=temp.z*scale;
                frame.camera_visible=temp.camera_mask;
                frame.residual=(float)temp.residual;}}
        delete[] temp_pts;}
    
    printf("Read %d channels for %d frames\n",number_channels,number_frames);
}




