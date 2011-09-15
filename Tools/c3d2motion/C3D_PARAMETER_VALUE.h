//##########################################################################################
// Class  C3D_PARAMETER_VALUE
// A value in the parameter table
//##########################################################################################
#ifndef _C3D_PARAMETER_VALUE_h
#define _C3D_PARAMETER_VALUE_h
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include "C3D_PARAMETER_ENTRY.h"
#include "C3D_PARSE_ERROR.h"
#include <assert.h>
#include <stdio.h>

// magic c3d types
typedef unsigned char C3D_BYTE;
typedef char C3D_SIGNED_BYTE;
typedef char C3D_CHAR;
typedef float C3D_FLOAT;
typedef short C3D_INT;

class C3D_PARAMETER_VALUE:public C3D_PARAMETER_ENTRY{
    int dimension;
    C3D_SIGNED_BYTE type;
    C3D_BYTE lengths[7];
    // We need this so we can grab it whatever way we want it.
    union{
        C3D_FLOAT* dat_float;
        C3D_INT* dat_int;
        C3D_CHAR* dat_char;
        C3D_BYTE* dat_byte;
    }data;
    // Type codes
    static const C3D_SIGNED_BYTE    C3D_PARAM_TYPE_CHAR = -1;
    static const C3D_SIGNED_BYTE    C3D_PARAM_TYPE_BYTE = 1;
    static const C3D_SIGNED_BYTE    C3D_PARAM_TYPE_INT = 2;
    static const C3D_SIGNED_BYTE    C3D_PARAM_TYPE_FLOAT = 4;
    // Machine codes
    static const C3D_BYTE    C3D_MACHINE_INTEL = 84;
    static const C3D_BYTE    C3D_MACHINE_DEC = 85;
    static const C3D_BYTE    C3D_MACHINE_MIPS = 86;
public:
    
    C3D_PARAMETER_VALUE()
        :C3D_PARAMETER_ENTRY()
    { 
        dimension=0;memset(lengths,0,7);data.dat_byte=0;
    }
    
    virtual ~C3D_PARAMETER_VALUE()
    {if(data.dat_byte) delete[] data.dat_byte;}
    
    //C3D_PARAMETER_VALUE& operator=(const C3D_PARAMETER_VALUE& input)
    //{if(this==&input) return;
    //dimension=input.dimension;type=input.type;memcpy(lengths,input.lengths,7*sizeof(C3D_BYTE));
    //delete dat_byte;dat_byte=new C3D_BYTE[number_entries*val_size];memcpy(dat_byte,input.data.dat_byte,number_entries*val_size);}
    
    int Get_Number_Entries() const //Gets the number of entries total i.e. prod(lengths[i],i=0,...,dim-1)
    {int entries=1;for(int i=0;i<dimension;i++){entries*=lengths[i];}return entries;}
    
    int Value_Size() const
    {switch(type){
    case C3D_PARAM_TYPE_BYTE: return sizeof(C3D_BYTE);
    case C3D_PARAM_TYPE_CHAR: return sizeof(C3D_CHAR);
    case C3D_PARAM_TYPE_INT: return sizeof(C3D_INT);
    case C3D_PARAM_TYPE_FLOAT: return sizeof(C3D_FLOAT);
    default: throw C3D_PARSE_ERROR("We do not understand this type");}}
    
    int Get_Number_Rows() //Calculate how many entries in the major dimension
    {return dimension>0? lengths[dimension-1] : 0;}
    
    int Get_Row_Size() //calculate how many entries per major row 
    {if(dimension>0){int rowsize=1;for (int i=0;i<dimension-1;i++) rowsize*=lengths[i];return rowsize;}
    else return 0;}
        
    int Get_Int()
    {assert(dimension==1 && lengths[0]==1 && type==C3D_PARAM_TYPE_INT);
    return data.dat_int[0];}
    
    float Get_Float()
    {assert(dimension==1 && lengths[0]==1 && type==C3D_PARAM_TYPE_FLOAT);
    return data.dat_float[0];}
    
    //lengths[1] - number of strings
    //lengths[0] - length of each string
    void Fill_String_Array(PhysBAM::ARRAY<std::string>& array)
    {if(this==0) return;
        assert(dimension==2 && type==C3D_PARAM_TYPE_CHAR);
    const int number_strings=lengths[1],string_length=lengths[0];array.Resize(number_strings);
    for(int i=1;i<=number_strings;++i) array(i)=std::string(data.dat_char+(i-1)*string_length,data.dat_char+i*string_length);}
    
//##########################################################################################
    void Read_Data_Part(FILE *fp,const C3D_BYTE original_machine_code);
    void Display_Data();
//##########################################################################################
};
#endif
