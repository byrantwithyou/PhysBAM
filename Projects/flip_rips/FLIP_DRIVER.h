//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLIP_DRIVER__
#define __FLIP_DRIVER__

#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <Tools/Log/LOG.h>
#include "FLIP_COLLIDABLE_OBJECT.h"
#include "FLIP_KRYLOV_VECTOR.h"
#include "FLIP_KRYLOV_SYSTEM.h"
#include "FLIP_ITERATOR.h"
#include <boost/format.hpp>
#include <fstream>

namespace PhysBAM{

template<class TV,int order>
class FLIP_DRIVER
{
public:

    enum WORKAROUND{w=order+1,ghost=3};
    enum CELL_TYPE{
        CELL_NEUMANN=0,
        CELL_DIRICHLET=1,
        CELL_INTERIOR=2};

    typedef typename TV::SCALAR T;
    typedef GRID<TV> T_GRID;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef FLIP_PARTICLE<TV,w> T_PARTICLE;
    typedef FLIP_ITERATOR<TV,w> T_ITERATOR;
    typedef ARRAY<T,TV_INT> T_CELL_ARRAY;
    typedef ARRAY<short int,TV_INT> INT_CELL_ARRAY;
    typedef FACE_INDEX<TV::m> T_FACE;
    typedef ARRAY<T,T_FACE> T_FACE_ARRAY;
    typedef FACE_ITERATOR<TV> T_FACE_ITERATOR;
    typedef CELL_ITERATOR<TV> T_CELL_ITERATOR;
    
    const T mass_threshold;
    const T dx;
    const T one_over_dx;
    const T_GRID grid;
    const int threads;

    std::string output_directory;
    const TV gravity;
    
    T cfl;
    T flip;

    int frames;
    T frame_dt;
    T time;
    T min_dt;
    T max_dt;

    bool frame_done;
    T new_time; // end of step time

    int output_number;

    T_FACE_ARRAY velocity,velocity_new,mass;
    INT_CELL_ARRAY cell_type; // 0 - neumann, 1 - dirichlet, 2 - interior
    T_CELL_ARRAY density;
    T_CELL_ARRAY cell_mass;

    ARRAY<T_PARTICLE> particles;
    ARRAY<TV_INT> interior_cells;

    FLIP_KRYLOV_SYSTEM<TV> system;
    FLIP_KRYLOV_VECTOR<TV> pressure,rhs,solver_q,solver_s,solver_r,solver_k,solver_z;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> solver_av;

    ARRAY<FLIP_COLLIDABLE_OBJECT<TV>*> collidable_objects;

    FLIP_DRIVER(const T_GRID grid_input,const int threads_input)
        :mass_threshold(1e-4),dx(grid_input.dX.x),
        one_over_dx((T)1/dx),grid(grid_input),threads(threads_input),
        output_directory("output"),gravity(TV::Axis_Vector(1)*(T)(-9.8)),
        cfl((T).2),flip(0.95),frames(100),frame_dt((T)1/24),
        time(0),min_dt(1e-5),max_dt(1e-3),output_number(0),
        system(one_over_dx,threads,interior_cells,cell_type,mass),
        pressure(interior_cells,threads),rhs(interior_cells,threads),
        solver_q(interior_cells,threads),solver_s(interior_cells,threads),
        solver_r(interior_cells,threads),solver_k(interior_cells,threads),
        solver_z(interior_cells,threads)
    {
        if(grid.dX.Min()!=grid.dX.Max()) PHYSBAM_FATAL_ERROR();
        solver_av.Remove_All();
        solver_av.Append(&solver_q);
        solver_av.Append(&solver_s);
        solver_av.Append(&solver_r);
        solver_av.Append(&solver_k);
        solver_av.Append(&solver_z);
    }

    ~FLIP_DRIVER(){}

    void Run()
    {
        LOG::cout<<"Initializing..."<<std::endl;
        Initialize();
        Write("Frame 0");
        for(int k=1;k<=frames;k++){
            LOG::cout<<"Frame "<<k<<std::endl;
            Advance_Frame();
            Write((boost::format("Frame %d")%k).str().c_str());
        }
    }
    
    void Initialize()
    {
        mass.Resize(grid,ghost);
        velocity.Resize(grid,ghost);
        velocity_new.Resize(grid,ghost);
        cell_type.Resize(grid.Cell_Indices(ghost));
        density.Resize(grid.Cell_Indices(ghost));
        cell_mass.Resize(grid.Cell_Indices(ghost));

        rhs.p.Resize(grid.Cell_Indices(ghost));
        pressure.p.Resize(grid.Cell_Indices(ghost));
        solver_q.p.Resize(grid.Cell_Indices(ghost));
        solver_s.p.Resize(grid.Cell_Indices(ghost));
        solver_r.p.Resize(grid.Cell_Indices(ghost));
        solver_k.p.Resize(grid.Cell_Indices(ghost));
        solver_z.p.Resize(grid.Cell_Indices(ghost));
        
        Initialize_Grid_Based_Variables();
        Update_Particle_Weights_And_Cell_Info();
        Rasterize_Mass_And_Velocity_To_Grid();
    }

    void Initialize_Grid_Based_Variables()
    {
        for(int j=0;j<TV::m;j++)
#pragma omp parallel for
        for(int i=0;i<mass.data(j).array.m;i++){
            mass.data(j).array(i)=T();
            velocity.data(j).array(i)=T();
            velocity_new.data(j).array(i)=T();}

#pragma omp parallel for
        for(int i=0;i<pressure.p.array.m;i++){
            cell_type.array(i)=CELL_DIRICHLET;
            cell_mass.array(i)=T();
            pressure.p.array(i)=T();
            rhs.p.array(i)=T();}

        interior_cells.Remove_All();
    }

    void Update_Particle_Weights_And_Cell_Info()
    {
#pragma omp parallel for
        for(int i=0;i<particles.m;i++)
            particles(i).Update_Base_And_Weights(grid);
    }

    void Rasterize_Mass_And_Velocity_To_Grid()
    {
        // do mass and momentum rasterization
        for(int i=0;i<particles.m;i++){
            for(int j=0;j<TV::m;j++)
                for(T_ITERATOR it(particles(i),j);it.Valid();it.Next()){
                    const T weight=it.Weight();
                    const T_FACE face(j,it.Index());
                    velocity(face)+=particles(i).V(j)*weight;
                    mass(face)+=weight;}
            for(T_ITERATOR it(particles(i),TV::m);it.Valid();it.Next()){
                const T weight=it.Weight();
                const TV_INT index=it.Index();
                cell_mass(index)+=weight;}}

        // divide off the mass
        for(int j=0;j<TV::m;j++)
#pragma omp parallel for
        for(int i=0;i<mass.data(j).array.m;i++)
            if(mass.data(j).array(i)>mass_threshold)
                velocity.data(j).array(i)/=mass.data(j).array(i);
            else{
                velocity.data(j).array(i)=T();
                mass.data(j).array(i)=T();}

#pragma omp parallel for
        for(int i=0;i<cell_mass.array.m;i++)
            if(cell_mass.array(i)<=mass_threshold)
                cell_mass.array(i)=T();
    }

    void Mark_Grid_Cells_And_Perform_Grid_Based_Collisions()
    {
        for(T_CELL_ITERATOR it(grid,ghost);it.Valid();it.Next()){
            const TV_INT index=it.Cell_Index();
            bool interior=true;
            bool neumann=true;
            if(!cell_mass(index)) interior=false;
            for(int j=0;j<TV::m;j++)
            for(int s=0;s<=1;s++){
                T_FACE face(j,index+TV_INT::Axis_Vector(j)*s);
                if(!mass(face)) interior=false;
                bool collides=false;
                for(int k=0;k<collidable_objects.m;k++)
                    if(collidable_objects(k)->Detect_Collision(grid.Face(T_FACE(face.axis,face.index)),time)){
                        collides=true;velocity(face)=T();break;}
                if(!collides) neumann=false;}
            if(neumann){cell_type(index)=CELL_NEUMANN;continue;}
            if(interior){
                cell_type(index)=CELL_INTERIOR;
                interior_cells.Append(index);}}
    }

    void Advance_Frame()
    {
        frame_done=false;
        const T initial_time=time;
        const T final_time=time+frame_dt;

        while(!frame_done){

            // compute dt
            T max_v=Max_Particle_Velocity();
            T dt=max(min_dt,min(max_dt,cfl*dx/max(max_v,(T)1e-2)));
            LOG::cout<<"Max Particle Velocity:  "<<max_v<<std::endl;

            // update dt and new time
            if(final_time-time<dt*1.001){
                dt=final_time-time;
                new_time=final_time;
                frame_done=true;}
            else{
                if(final_time-time<2*dt)
                    dt=(final_time-time)*0.5;
                new_time=time+dt;}

            // advance step
            Advance_Step(dt);

            // end timestep
            time=new_time;
            LOG::cout<<(int)((time-initial_time)*1000/frame_dt)/(T)10<<"%\t  dt="<<dt<<std::endl;

            // Write((boost::format("step")).str().c_str());
        }
    }

    void Advance_Step(const T dt)
    {
        Apply_Forces(dt);
        Initialize_Grid_Based_Variables();
        Update_Particle_Weights_And_Cell_Info();
        Rasterize_Mass_And_Velocity_To_Grid();
        Mark_Grid_Cells_And_Perform_Grid_Based_Collisions();
        Projection(dt);
        Update_Particle_Positions_And_Velocities(dt);
    }
    
    void Apply_Forces(const T dt)
    {
#pragma omp parallel for
        for(int i=0;i<particles.m;i++)
            particles(i).V+=gravity*dt;
    }

    void Dump_System_Matrix(const char* filename)
    {
        std::ofstream fout(filename);
        const int size=solver_k.Raw_Size();
        for(int j=0;j<size;j++)
        {
            for(int i=0;i<size;i++)
                solver_k.Raw_Get(i)=(i==j);
            system.Multiply(solver_k,solver_z);
            for(int i=0;i<size;i++)
                if(solver_z.Raw_Get(i))
                    fout<<i<<" "<<j<<" "<<solver_z.Raw_Get(i)<<std::endl;
        }
    }

    void Dump_Preconditioner(const char* filename)
    {
        std::ofstream fout(filename);
        const int size=solver_k.Raw_Size();
        for(int j=0;j<size;j++)
        {
            for(int i=0;i<size;i++)
                solver_k.Raw_Get(i)=(i==j);
            system.Precondition(solver_k,solver_z);
            for(int i=0;i<size;i++)
                if(solver_z.Raw_Get(i))
                    fout<<i<<" "<<j<<" "<<solver_z.Raw_Get(i)<<std::endl;
        }
    }

    void Projection(const T dt)
    {
        // fill in the rhs
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            T divergence=0;
            for(int j=0;j<TV::m;j++)
                divergence+=one_over_dx*(
                    velocity(T_FACE(j,index+TV_INT::Axis_Vector(j)))-
                    velocity(T_FACE(j,index)));
            rhs.p(index)=divergence;}

        // solve for pressure
        CONJUGATE_RESIDUAL<T> cr;
        system.dt=dt;
        // Dump_System_Matrix("A.txt");
        // Dump_Preconditioner("M.txt");
        cr.Solve(system,pressure,rhs,solver_av,1e-7,0,10000);
    }

    void Perform_Particle_Collision(T_PARTICLE& particle)
    {
        for(int i=0;i<collidable_objects.m;i++)
            if(collidable_objects(i)->Detect_Collision(particle.X,time)) particle.V=TV();
    }

    void Update_Particle_Positions_And_Velocities(const T dt)
    {
#pragma omp parallel for
        for(int i=0;i<particles.m;i++){
            TV V_pic,V_flip=particles(i).V;
            T div=0;
            for(int j=0;j<TV::m;j++)
            for(T_ITERATOR it(particles(i),j);it.Valid();it.Next())
            {
                const T weight=it.Weight();
                const T_FACE face(j,it.Index());
                const TV_INT cell0=face.index-TV_INT::Axis_Vector(j);
                const TV_INT cell1=face.index; 
                const short int type0=cell_type(cell0);
                const short int type1=cell_type(cell1);

                T grad_p=0;
                if(type0==CELL_INTERIOR && type1==CELL_INTERIOR)
                    grad_p=(pressure.p(cell1)-pressure.p(cell0))*one_over_dx;
                else if(type0==CELL_INTERIOR && type1==CELL_DIRICHLET)
                    grad_p=(-pressure.p(cell0))*one_over_dx;
                else if(type0==CELL_DIRICHLET && type1==CELL_INTERIOR)
                    grad_p=(pressure.p(cell1))*one_over_dx;

                if(mass(face))
                {
                    T mass_scale=1/mass(face);
                    V_pic(j)+=(velocity(face)+grad_p*mass_scale*dt)*weight;
                    V_flip(j)+=grad_p*mass_scale*dt*weight;
                    div+=(velocity(face)+grad_p*mass_scale*dt)*it.Weight_Gradient()(j);}
                }

            particles(i).V=V_flip*flip+V_pic*(1-flip);
            Perform_Particle_Collision(particles(i));
            particles(i).X+=V_pic*dt;
        }
    }
 
    T Max_Particle_Velocity() const
    {
        ARRAY<T> result_per_thread(threads);
#pragma omp parallel for
        for(int i=0;i<particles.m;i++){
            const int tid=omp_get_thread_num();
            T& r=result_per_thread(tid);
            r=max(r,particles(i).V.Magnitude_Squared());}
        T result=0;
        for(int tid=0;tid<threads;tid++)
            result=max(result,result_per_thread(tid));
        return sqrt(result);
    }

    void Write(const char* text)
    {
        FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        FILE_UTILITIES::Create_Directory((boost::format("%s/%d")%output_directory%output_number).str());
        
        FILE_UTILITIES::Write_To_File<T>((boost::format("%s/%d/mac_velocities")%output_directory%output_number).str(),velocity);

        FILE_UTILITIES::Write_To_File<T>((boost::format("%s/%d/time")%output_directory%output_number).str(),time);
        FILE_UTILITIES::Write_To_File<T>(output_directory+"/common/grid",grid);
        FILE_UTILITIES::Write_To_File<T>((boost::format("%s/%d/grid")%output_directory%output_number).str(),grid);

        for(T_CELL_ITERATOR it(grid,ghost);it.Valid();it.Next()){
            switch(cell_type(it.Cell_Index())){
            case CELL_DIRICHLET: density(it.Cell_Index())=0; break;
            case CELL_NEUMANN: density(it.Cell_Index())=0.333; break;
            case CELL_INTERIOR: density(it.Cell_Index())=0.666; break;
            default: PHYSBAM_FATAL_ERROR();}}
        FILE_UTILITIES::Write_To_File<T>((boost::format("%s/%d/density.gz")%output_directory%output_number).str(),density);

        GEOMETRY_PARTICLES<TV> geometry_particles;
        geometry_particles.Store_Velocity();
        for(int i=0;i<particles.m;i++){
            const int out_index=geometry_particles.Add_Element();
            geometry_particles.X(out_index)=particles(i).X;
            geometry_particles.V(out_index)=particles(i).V;}

        FILE_UTILITIES::Write_To_File<T>((boost::format("%s/%d/debug_particles")%output_directory%output_number).str(),geometry_particles);
        
        FILE_UTILITIES::Write_To_Text_File((boost::format("%s/%d/frame_title")%output_directory%output_number).str(),text);
        FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",output_number);
        output_number++;
    }
};
}

#endif
