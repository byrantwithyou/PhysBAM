//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Struct RGB_HEADER
//#####################################################################
// A group of fields with the layout of a Windows rgb file header. Used by the class RGB_FILE.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __RGB_HEADER__
#define __RGB_HEADER__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <cstring>
namespace PhysBAM{

struct RGB_HEADER
{
    typedef int HAS_UNTYPED_READ_WRITE;

    unsigned short magic_number; // must be 474
    unsigned char compression; // 0 or 1
    unsigned char bytes_per_channel; // 1
    unsigned short dimensions; // 2 for a single 2D image
    unsigned short width,height;
    unsigned short channels; // 3 is RGB, 4 is RGBA
    unsigned int pixel_minimum_value; // typically 0
    unsigned int pixel_maximum_value; // typically 255
    unsigned int dummy; // ignored
    char name[80];
    unsigned int colormap;
    char dummy2[404];
    
    void Initialize(const int width_input,const int height_input)
    {magic_number=474;compression=0;bytes_per_channel=1;dimensions=2;width=width_input;height=height_input;channels=3;pixel_maximum_value=0;pixel_maximum_value=255;strcpy(name,"PhysBAM");colormap=0;}

    void Swap_Endian() const
    {Swap_Endianity(magic_number);Swap_Endianity(compression);Swap_Endianity(bytes_per_channel);Swap_Endianity(dimensions);
    Swap_Endianity(width);Swap_Endianity(height);Swap_Endianity(channels);Swap_Endianity(pixel_minimum_value);Swap_Endianity(pixel_maximum_value);Swap_Endianity(colormap);}        

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,magic_number,compression,bytes_per_channel,dimensions,width,height);
    Read_Binary<RW>(input,channels,pixel_minimum_value,pixel_maximum_value,dummy);
    Read_Binary_Array<RW>(input,name,80);Read_Binary<RW>(input,colormap);Read_Binary_Array<RW>(input,dummy2,404);
    Swap_Endian();
    
    //check validity
    if(magic_number!=474) throw READ_ERROR(STRING_UTILITIES::string_sprintf("illegal file type %d, should be 474",magic_number));
    if(channels != 3) throw READ_ERROR(STRING_UTILITIES::string_sprintf("illegal number of channels: %d (must be 3)",channels));
    if(bytes_per_channel != 1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("only 1 bytes per channel supported (you gave %d)",bytes_per_channel));
    if(dimensions != 2 && dimensions != 3) throw READ_ERROR(STRING_UTILITIES::string_sprintf("dimensions field must be 2 or 3 (you gave %d)",dimensions));}

    template<class RW> void Write(std::ostream& output) const
    {RGB_HEADER swapped=*this;swapped.Swap_Endian();
    Write_Binary<RW>(output,swapped.magic_number,swapped.compression,swapped.bytes_per_channel,swapped.dimensions,swapped.width,swapped.height,swapped.channels);
    Write_Binary<RW>(output,swapped.pixel_minimum_value,swapped.pixel_maximum_value,swapped.dummy);
    Write_Binary_Array<RW>(output,name,80);Write_Binary<RW>(output,colormap);Write_Binary_Array<RW>(output,dummy2,404);}
};
}
#endif
#endif
