//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
using namespace PhysBAM;
OPENGL_INDEXED_COLOR_MAP::
OPENGL_INDEXED_COLOR_MAP()
    :index_mode(PERIODIC)
{
}
OPENGL_INDEXED_COLOR_MAP::
~OPENGL_INDEXED_COLOR_MAP()
{
}
OPENGL_COLOR OPENGL_INDEXED_COLOR_MAP::
Lookup(int index) const
{
    if(index_mode==PERIODIC){
        int i=index%color_map.m;
        if(i<0) i+=color_map.m;
        return color_map(i);}
    if(index>=color_map.m) return color_map.Last();
    if(index<=0) return color_map(0);
    return color_map(index);
}
OPENGL_INDEXED_COLOR_MAP* OPENGL_INDEXED_COLOR_MAP::
Basic_16_Color_Map()
{
    OPENGL_INDEXED_COLOR_MAP* map=new OPENGL_INDEXED_COLOR_MAP;
    map->Set_Color(0,OPENGL_COLOR::Black());
    map->Set_Color(1,OPENGL_COLOR::Red(0.5));
    map->Set_Color(2,OPENGL_COLOR::Green(0.5));
    map->Set_Color(3,OPENGL_COLOR::Yellow(0.5));
    map->Set_Color(4,OPENGL_COLOR::Blue(0.5));
    map->Set_Color(5,OPENGL_COLOR::Magenta(0.5));
    map->Set_Color(6,OPENGL_COLOR::Cyan(0.5));
    map->Set_Color(7,OPENGL_COLOR::Gray(0.25));
    map->Set_Color(8,OPENGL_COLOR::Gray(0.75));
    map->Set_Color(9,OPENGL_COLOR::Red());
    map->Set_Color(10,OPENGL_COLOR::Green());
    map->Set_Color(11,OPENGL_COLOR::Yellow());
    map->Set_Color(12,OPENGL_COLOR::Blue());
    map->Set_Color(13,OPENGL_COLOR::Magenta());
    map->Set_Color(14,OPENGL_COLOR::Cyan());
    map->Set_Color(15,OPENGL_COLOR::White());
    return map;
}
OPENGL_INDEXED_COLOR_MAP* OPENGL_INDEXED_COLOR_MAP::
Levelset_Multiple_Color_Map()
{
    OPENGL_INDEXED_COLOR_MAP* map=new OPENGL_INDEXED_COLOR_MAP;
    map->Set_Color(0,OPENGL_COLOR::Red(0.6));
    map->Set_Color(1,OPENGL_COLOR::Green(0.6));
    map->Set_Color(2,OPENGL_COLOR::Blue(0.6));
    map->Set_Color(3,OPENGL_COLOR::Magenta(0.6));
    map->Set_Color(4,OPENGL_COLOR::Cyan(0.6));
    map->Set_Color(5,OPENGL_COLOR::Gray(0.25));
    map->Set_Color(6,OPENGL_COLOR::Gray(0.75));
    map->Set_Color(7,OPENGL_COLOR::Yellow(0.6));
    return map;
}
OPENGL_INDEXED_COLOR_MAP* OPENGL_INDEXED_COLOR_MAP::
Particle_Multiple_Color_Map()
{
    OPENGL_INDEXED_COLOR_MAP* map=Levelset_Multiple_Color_Map();
    for(int i=0;i<map->color_map.m;i++){for(int j=0;j<3;j++) map->color_map(i).rgba[j]+=.35;map->color_map(i).rgba[3]=1;}
    return map;
}
OPENGL_INDEXED_COLOR_MAP* OPENGL_INDEXED_COLOR_MAP::
Rigid_Body_Color_Map()
{
    OPENGL_INDEXED_COLOR_MAP* map=new OPENGL_INDEXED_COLOR_MAP;
    map->Set_Color(0,OPENGL_COLOR(1,72.f/255,72.f/255));
    map->Set_Color(1,OPENGL_COLOR(1,128.f/255,72.f/255));
    map->Set_Color(2,OPENGL_COLOR(1,191.f/255,72.f/255));
    map->Set_Color(3,OPENGL_COLOR(1,239.f/255,72.f/255));
    map->Set_Color(4,OPENGL_COLOR(172.f/255,239.f/255,72.f/255));
    map->Set_Color(5,OPENGL_COLOR(80.f/255,239.f/255,72.f/255));
    map->Set_Color(6,OPENGL_COLOR(80.f/255,239.f/255,1));
    map->Set_Color(7,OPENGL_COLOR(80.f/255,172.f/255,1));
    map->Set_Color(8,OPENGL_COLOR(80.f/255,102.f/255,1));
    map->Set_Color(9,OPENGL_COLOR(80.f/255,21.f/255,1));
    map->Set_Color(10,OPENGL_COLOR(138.f/255,21.f/255,1));
    map->Set_Color(11,OPENGL_COLOR(200.f/255,21.f/255,1));
    return map;
}
OPENGL_INDEXED_COLOR_MAP* OPENGL_INDEXED_COLOR_MAP::
Rigid_Body_Back_Color_Map()
{
    OPENGL_INDEXED_COLOR_MAP* map=Rigid_Body_Color_Map();
    for(int i=0;i<map->color_map.m;i++){for(int j=0;j<3;j++) map->color_map(i).rgba[j]+=.35;map->color_map(i).rgba[3]=1;}
    return map;
}
