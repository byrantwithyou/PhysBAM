#ifndef __REQUEST_LOGGING__
#define __REQUEST_LOGGING__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include "libmain.h"

extern "C" void enable_logging(const char* file);
extern "C" void finish_logging();

struct request_log_simulation_info
{
    int id;
    int next_object;
    std::map<const physbam_base*,int> object_map;

    request_log_simulation_info();
    ~request_log_simulation_info();
};

struct request_log_info
{
    std::string request_log_file;
    std::ofstream out;
    int next_simulation_id;
    std::map<physbam_simulation*,request_log_simulation_info> simulation_map;

    request_log_info(const char* file);
    ~request_log_info();
};

extern request_log_info * log_info;

template<class T,int d> inline void
Read(std::istream& input,std::vector<data_exchange::fixed_vector<T,d> >& v)
{
    int size = -1;
    PhysBAM::Read_Binary<float>(input,size);
    v.resize(size);
    PhysBAM::Read_Binary_Array<float>(input,(T*)&v[0],v.size()*d);
}

template<class T,int d> inline void
Write(std::ostream& output,const std::vector<data_exchange::fixed_vector<T,d> >& v)
{
    PhysBAM::Write_Binary<float>(output,(int)v.size());
    PhysBAM::Write_Binary_Array<float>(output,(T*)&v[0],v.size()*d);
}

inline void
Read(std::istream& input,std::vector<int>& v)
{
    int size = -1;
    PhysBAM::Read_Binary<float>(input,size);
    v.resize(size);
    PhysBAM::Read_Binary_Array<float>(input,&v[0],v.size());
}

inline void
Write(std::ostream& output,const std::vector<int>& v)
{
    PhysBAM::Write_Binary<float>(output,(int)v.size());
    PhysBAM::Write_Binary_Array<float>(output,&v[0],v.size());
}

template<class T,int d> inline void
Read(std::istream& input,data_exchange::fixed_vector<T,d>& v)
{
    PhysBAM::Read_Binary_Array<float>(input,v.data,d);
}

template<class T,int d> inline void
Write(std::ostream& output,const data_exchange::fixed_vector<T,d>& v)
{
    PhysBAM::Write_Binary_Array<float>(output,v.data,d);
}

inline void
Read(std::istream& input,data_exchange::polygon_mesh& pm)
{
    Read(input,pm.polygons);
    Read(input,pm.polygon_counts);
}

inline void
Write(std::ostream& output,const data_exchange::polygon_mesh& pm)
{
    Write(output,pm.polygons);
    Write(output,pm.polygon_counts);
}

void Read(std::istream& input,data_exchange::simulation_object *& so)
{
    PhysBAM::Read_Binary<float>(input,so->id);
    if(so->id==data_exchange::deformable_body::fixed_id())
    {
        data_exchange::deformable_body * body = new data_exchange::deformable_body;
        so = body;
        Read(input,body->mesh);
        Read(input,body->position);
        Read(input,body->velocity);
        PhysBAM::Read_Binary<float>(input,body->mass);
    }
    else if(so->id==data_exchange::ground_plane::fixed_id())
    {
        data_exchange::ground_plane * body = new data_exchange::ground_plane;
        so = body;
        Read(input,body->position);
        Read(input,body->normal);
    }
    else if(so->id==data_exchange::scripted_geometry::fixed_id())
    {
        data_exchange::scripted_geometry * body = new data_exchange::scripted_geometry;
        so = body;
        Read(input,body->mesh);
        Read(input,body->position);
    }
}

void Write(std::ostream& output,const data_exchange::simulation_object * so)
{
    PhysBAM::Write_Binary<float>(output,so->id);
    if(so->id==data_exchange::deformable_body::fixed_id())
    {
        const data_exchange::deformable_body * body = static_cast<const data_exchange::deformable_body*>(so);
        Write(output,body->mesh);
        Write(output,body->position);
        Write(output,body->velocity);
        PhysBAM::Write_Binary<float>(output,body->mass);
    }
    else if(so->id==data_exchange::ground_plane::fixed_id())
    {
        const data_exchange::ground_plane * body = static_cast<const data_exchange::ground_plane*>(so);
        Write(output,body->position);
        Write(output,body->normal);
    }
    else if(so->id==data_exchange::scripted_geometry::fixed_id())
    {
        const data_exchange::scripted_geometry * body = static_cast<const data_exchange::scripted_geometry*>(so);
        Write(output,body->mesh);
        Write(output,body->position);
    }
}

void Read(std::istream& input,data_exchange::force *& f)
{
    PhysBAM::Read_Binary<float>(input,f->id);
    if(f->id==data_exchange::gravity_force::fixed_id())
    {
        data_exchange::gravity_force * fo = new data_exchange::gravity_force;
        f = fo;
        PhysBAM::Read_Binary<float>(input,fo->magnitude);
        Read(input,fo->direction);
    }
    else if(f->id==data_exchange::volumetric_force::fixed_id())
    {
        data_exchange::volumetric_force * fo = new data_exchange::volumetric_force;
        f = fo;
        PhysBAM::Read_Binary<float>(input,fo->stiffness,fo->poissons_ratio,fo->damping);
    }
}

void Write(std::ostream& output,const data_exchange::force * f)
{
    PhysBAM::Write_Binary<float>(output,f->id);
    if(f->id==data_exchange::gravity_force::fixed_id())
    {
        const data_exchange::gravity_force * fo = static_cast<const data_exchange::gravity_force*>(f);
        PhysBAM::Write_Binary<float>(output,fo->magnitude);
        Write(output,fo->direction);
    }
    else if(f->id==data_exchange::volumetric_force::fixed_id())
    {
        const data_exchange::volumetric_force * fo = static_cast<const data_exchange::volumetric_force*>(f);
        PhysBAM::Write_Binary<float>(output,fo->stiffness,fo->poissons_ratio,fo->damping);
    }
}

#endif

