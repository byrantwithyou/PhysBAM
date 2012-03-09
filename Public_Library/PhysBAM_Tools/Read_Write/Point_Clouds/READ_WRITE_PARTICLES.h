//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_PARTICLES
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_PARTICLES__
#define __READ_WRITE_PARTICLES__

#include <PhysBAM_Tools/Point_Clouds/PARTICLES.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<PARTICLES<TV>,RW>
{
public:
    static void Read(std::istream& input,PARTICLES<TV>& object)
    {int version;
    Read_Binary<RW>(input,version);
    if(version!=1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized particle version %d",(int)version));
    Read_Binary<RW>(input,*object.array_collection);}

    static void Write(std::ostream& output,const PARTICLES<TV>& object)
    {Write_Binary<RW>(output,1,*object.array_collection);}

    static void Print(std::ostream& output,const PARTICLES<TV>& object,const int p)
    {Read_Write<ARRAY_COLLECTION,RW>::Print(output,*object.array_collection,p);}
};

template<class RW,class T_PARTICLES>
class Read_Write<T_PARTICLES,RW,typename ENABLE_IF<AND<IS_BASE_OF<PARTICLES<typename T_PARTICLES::VECTOR_T> ,T_PARTICLES>::value,NOT<IS_SAME<T_PARTICLES,PARTICLES<typename T_PARTICLES::VECTOR_T> >::value>::value>::value>::TYPE>:public Read_Write<PARTICLES<typename T_PARTICLES::VECTOR_T>,RW>
{
};
}
#endif
#endif
