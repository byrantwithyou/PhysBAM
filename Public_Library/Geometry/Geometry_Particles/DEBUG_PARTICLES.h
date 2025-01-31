//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEBUG_PARTICLES
//#####################################################################
#ifndef __DEBUG_PARTICLES__
#define __DEBUG_PARTICLES__

#include <Core/Arrays/ATTRIBUTE_ID.h>
#include <Core/Read_Write/TYPED_STREAM.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV> class GEOMETRY_PARTICLES;
template<class TV> class GRID;
template<class TV> class IMPLICIT_OBJECT;
class VIEWER_DIR;

template<class TV>
struct DEBUG_OBJECT
{
    typedef typename TV::SCALAR T;
    enum TYPE {segment=2,triangle=3} type;
    VECTOR<TV,3> X;
    VECTOR<T,3> color,bgcolor;
    T separation;
    bool draw_vertices;

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,X,type,color,bgcolor,separation,draw_vertices);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,X,type,color,bgcolor,separation,draw_vertices);}
};

template<class TV>
struct DEBUG_TEXT
{
    typedef typename TV::SCALAR T;
    TV X;
    std::string text;
    VECTOR<T,3> color;

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,X,text,color);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,X,text,color);}
};

template<class TV>
class DEBUG_PARTICLES
{
public:
    typedef typename TV::SCALAR T;
    DEBUG_PARTICLES();
    ~DEBUG_PARTICLES();

    GEOMETRY_PARTICLES<TV>& debug_particles;
    mutable ARRAY<DEBUG_OBJECT<TV> > debug_objects;
    mutable ARRAY<DEBUG_TEXT<TV> > debug_text;
    T edge_separation;

    static DEBUG_PARTICLES<TV>* Store_Debug_Particles(DEBUG_PARTICLES<TV>* particle=0);
    void Write_Debug_Particles(STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const;
    void Clear_Debug_Particles() const;
};
template<class TV> DEBUG_PARTICLES<TV>& Get_Debug_Particles();
template<class TV,class ATTR> void Debug_Particle_Set_Attribute(const std::string& name,const ATTR& attr);
template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);
template<class TV,int d> inline void Add_Debug_Object(const VECTOR<TV,d>& object,const VECTOR<typename TV::SCALAR,3>& color){Add_Debug_Object(object,color,color);}
template<class TV,int d> void Add_Debug_Object(const VECTOR<TV,d>& object,const VECTOR<typename TV::SCALAR,3>& color,const VECTOR<typename TV::SCALAR,3>& bgcolor);
template<class T_SURFACE,class T> void Dump_Surface(const T_SURFACE& surface,const VECTOR<T,3>& color){Dump_Surface(surface,color,color);}
template<class TV,class TV_INT,class T> void Dump_Levelset(const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi,const VECTOR<T,3>& color){Dump_Levelset(grid,phi,color,color);}
template<class TV,class TV_INT,class T> void Dump_Levelset(const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi,const VECTOR<T,3>& color,T contour_value){Dump_Levelset(grid,phi,color,color,contour_value);}
template<class T_SURFACE,class T> void Dump_Surface(const T_SURFACE& surface,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor);
template<class TV,class TV_INT,class T> void Dump_Levelset(const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor);
template<class TV,class TV_INT,class T> void Dump_Levelset(const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor,T contour_value);
template<class TV,class T> void Dump_Levelset(const GRID<TV>& grid,const IMPLICIT_OBJECT<TV>& phi,const VECTOR<T,3>& color){Dump_Levelset(grid,phi,color,color);}
template<class TV,class T> void Dump_Levelset(const GRID<TV>& grid,const IMPLICIT_OBJECT<TV>& phi,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor);
template<class TV> void Add_Debug_Text(const TV& X,const std::string& text,const VECTOR<typename TV::SCALAR,3>& color);
}
#endif
