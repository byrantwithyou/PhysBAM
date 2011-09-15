//##########################################################################################
// File: C3D_PARAMETER_VALUE.cpp
//##########################################################################################
#include "C3D_FILE.h"
#include "C3D_PARAMETER_VALUE.h"
//##########################################################################################
// Function  Read_Data_Part
// Reads a data value using the magic...
//##########################################################################################
void C3D_PARAMETER_VALUE::
Read_Data_Part(FILE* fp,const C3D_BYTE original_machine_code)
{  
    assert(this);
    C3D_SIGNED_BYTE element_data_info[2];fread(element_data_info,sizeof(element_data_info[0]),2,fp);
    
    type=element_data_info[0];
    dimension=element_data_info[1];
    if(dimension>7) throw C3D_PARSE_ERROR("Cannot have dimension greater than 7");
    else if(dimension>0) fread(lengths,sizeof(lengths[0]),dimension,fp);
    else{dimension=1;lengths[0]=1;}
    if(dimension>2) printf("Wow dimension > 2\n");
    
    // Read in data
    const int number_entries=Get_Number_Entries(),val_size=Value_Size();
    data.dat_byte=new C3D_BYTE[number_entries*val_size]; //allocate space
    fread(data.dat_byte,sizeof(C3D_BYTE),val_size*number_entries,fp); //now do the deed
    
    // convert DEC "lame ass" format to "every sane person uses" intel format
    if (type==C3D_PARAM_TYPE_FLOAT && original_machine_code==C3D_MACHINE_DEC){
    for(int i=0;i<number_entries;i++){
        char c[4]; // temp storage
        c[0] = data.dat_char[i*4+2];
        c[1] = data.dat_char[i*4+3];
        c[2] = data.dat_char[i*4+0];
        c[3] = data.dat_char[i*4+1];
        if(c[0] || c[1] || c[2] || c[3]) --c[3]; //adjust exponent
        data.dat_float[i] = *(float*)c;}}
    else if (type==C3D_PARAM_TYPE_FLOAT && original_machine_code==C3D_MACHINE_MIPS)
        throw C3D_PARSE_ERROR("We don't know how to deal with MIPS format!!");
}
//##########################################################################################
// Function  Display_Data
//##########################################################################################
void C3D_PARAMETER_VALUE::
Display_Data()
{
    printf("  %-30s ",Get_Name());
    if(type==C3D_PARAM_TYPE_BYTE) printf("type=%s ","byte  ");
    else if(type==C3D_PARAM_TYPE_INT) printf("type=%s ","int   ");
    else if(type==C3D_PARAM_TYPE_FLOAT) printf("type=%s ","float ");
    else if(type==C3D_PARAM_TYPE_CHAR) printf("type=%s ","char  ");
    for(int i=0;i<dimension;i++){printf("%d",lengths[i]);if (i!=dimension-1) printf("x");}
    if(dimension==1 && lengths[0]==1){ 
        switch(type){
        case C3D_PARAM_TYPE_FLOAT: printf(" %f",data.dat_float[0]); break;
        case C3D_PARAM_TYPE_INT: printf(" %d",data.dat_int[0]); break;
        case C3D_PARAM_TYPE_CHAR: printf(" %c",data.dat_char[0]); break;
        case C3D_PARAM_TYPE_BYTE: printf(" %u",data.dat_byte[0]); break;}}
    printf("\n");
}








