//#####################################################################
// Copyright 2009, Eftychios Sifakis,Aleka McAdams
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIGRID_POISSON_TESTS
//#####################################################################
//
// Solving Lx=b
//
// WARNING: writing substeps is VERY slow!
//
// Test numbers:
// 1. Full square domain. Sinusoidal RHS
// 2. Test 1 with substep visualization output
// 3. Rectangle domain with Dirichlet boundaries and Neumann boundary sphere (e.g. flow past sphere). random RHS.
// 4. Test 3 with substep visualization output


#ifndef __MULTIGRID_POISSON_TESTS__
#define __MULTIGRID_POISSON_TESTS__

#include <Core/Log/SCOPE.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <memory>
#include "MG_PRECONDITIONED_CONJUGATE_GRADIENT.h"
#include "MULTIGRID_POISSON_SOLVER.h"
namespace PhysBAM{

template<class T,int d>
class MULTIGRID_POISSON_TESTS
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<int,d> TV_INT;

public:
    int test_number;
    std::string output_dir;
    int frames;
    bool write_substeps;
    ARRAY<T,TV_INT> x;
    ARRAY<T,TV_INT> b;
    ARRAY<T,TV_INT> tmp;
    std::unique_ptr<MULTIGRID_POISSON_SOLVER<T,d> > multigrid_poisson_solver;

    MULTIGRID_POISSON_TESTS(int test_number_input,int number_of_threads,int resolution=256)
        :test_number(test_number_input)
    {
        T_INDEX size=resolution*T_INDEX::All_Ones_Vector();
        int levels=1;
        switch(test_number){
            case 1:
            case 2:
                levels=(int)log2(resolution)-2;
                frames=1;
                break;
            case 3:
            case 4:
                size(2)*=1.5;
                levels=(int)log2(resolution)-3;
                frames=1;
                break;
            case 5:
            case 6:
                levels=(int)log2(resolution)-2;
                frames=2;
                break;
            default:
                LOG::cout<<"Unsupported test number"<<std::endl;
                exit(1);
                break;
        }
        
        if(test_number%2)
            write_substeps=false;
        else
            write_substeps=true;

        output_dir=LOG::sprintf("Test_%d_Resolution_%d",test_number,resolution);

        FILE_UTILITIES::Create_Directory(output_dir);
        FILE_UTILITIES::Create_Directory(output_dir+"/common");

        LOG::Initialize_Logging();
        LOG::Instance()->Copy_Log_To_File(output_dir+"/common/log.txt",false);
        LOG::cout<<"Running test number "<<test_number<<" at resolution "<<resolution<<" with "<<levels<<" levels"<<std::endl;
        multigrid_poisson_solver.reset(new MULTIGRID_POISSON_SOLVER<T,d>(size,(T)1/resolution,levels,number_of_threads));

    }

    ~MULTIGRID_POISSON_TESTS()
    {}

    void Initialize()
    {
        LOG::SCOPE scope("Test Initialization");
        RANDOM_NUMBERS<T> random_numbers;
        random_numbers.Set_Seed(1);
        // initialize cell types
        Set_Cell_Type(0);

        // initialize right hand side
        MULTIGRID_POISSON<T,d>& multigrid_poisson=multigrid_poisson_solver->Discretization();
        b.Resize(multigrid_poisson.grid.Domain_Indices(),false,false);
        switch(test_number){
            case 1:
            case 2:
                for(RANGE_ITERATOR<d> iterator(multigrid_poisson.padded_domain);iterator.Valid();iterator.Next()){
                    TV X=multigrid_poisson.grid.Node(iterator.index);
                    if(multigrid_poisson.cell_type(iterator.index)==MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE){
                        if(d==2)
                            b(iterator.index)=sin(X(1)*2*pi)*cos(X(2)*2*pi);
                        else if(d==3)
                            b(iterator.index)=sin(X(1)*2*pi)*cos(X(2)*2*pi)*sin(X(3)*2*pi);
                    }
                    else
                        b(iterator.index)=0;
                }
                break;        
            case 3:
            case 4:
                for(RANGE_ITERATOR<d> iterator(multigrid_poisson.padded_domain);iterator.Valid();iterator.Next()){
                    TV X=multigrid_poisson.grid.Node(iterator.index);
                    if(multigrid_poisson.cell_type(iterator.index)==MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE){ 
                        b(iterator.index)=random_numbers.Get_Number();
                    }
                    else
                        b(iterator.index)=0;
                }
                break;
            case 5:
            case 6:
                b.Fill(0);
                break;
        }

        // Initialize initial guess
        x.Resize(multigrid_poisson.grid.Domain_Indices(),false,false);
        tmp.Resize(multigrid_poisson.grid.Domain_Indices(),false,false);
        switch(test_number){
            case 1:
            case 2:
            case 3:
            case 4:
                for(RANGE_ITERATOR<d> iterator(multigrid_poisson.padded_domain);iterator.Valid();iterator.Next()){
                    const T_INDEX& index=iterator.index;
                    if(multigrid_poisson.cell_type(index)==MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE)
                        x(index)=random_numbers.Get_Number();
                    else
                        x(index)=0;
                }
                break;
            case 5:
            case 6:
                x.Fill(0);
                break;
        }

        // Initialize multigrid hierarchy
        multigrid_poisson_solver->Initialize_Multigrid_Hierarchy();                
    }

    void Run()
    {
        Initialize();
        Write_Frame(0);
        for(int frame=0;frame<frames;frame++){
            if(frame>1){
                Set_Cell_Type(frame);
                multigrid_poisson_solver->Reinitialize_Multigrid_Hierarchy();
            }
            Set_Boundary_Values(frame);

            multigrid_poisson_solver->Discretization().b=b;
            MULTIGRID_SYSTEM<T,d> multigrid_system(*multigrid_poisson_solver);
            MG_PRECONDITIONED_CONJUGATE_GRADIENT<T,d> cg;
            cg.print_residuals=true;
            
            // this is really slow. get substep data by actually restarting and solving with i iterations
            if(write_substeps){
                Write_Substep(0,frame);
                ARRAY<T,TV_INT> x_save(x);
                for(int i=0;i<100;i++){
                    x=x_save;
                    multigrid_poisson_solver->Discretization().b=b;
                    bool converged=cg.Solve(multigrid_system,x,multigrid_poisson_solver->Discretization().b,multigrid_poisson_solver->Discretization().u,tmp,1e-7,1,i);
                    Write_Substep(i,frame);
                    if(converged)
                        break;
                }        
            }else{

                multigrid_poisson_solver->Discretization().b=b;
                cg.Solve(multigrid_system,x,multigrid_poisson_solver->Discretization().b,multigrid_poisson_solver->Discretization().u,tmp,1e-7,1,50);
            }
            Write_Frame(frame);
            
        }

        LOG::Finish_Logging();

    }
    
    void Write_Frame(int frame)
    {
        LOG::SCOPE scope("Write frame","Write frame %d",frame);
        MULTIGRID_POISSON<T,d>& multigrid_poisson=multigrid_poisson_solver->Discretization();
        std::string f=LOG::sprintf("%i",frame);

        FILE_UTILITIES::Write_To_File<float>(output_dir+"/common/grid",multigrid_poisson.grid);
        ARRAY<T,TV_INT> x_as_density(x);
        T x_min=std::numeric_limits<T>::max();
        T x_max=-x_min;
        for(RANGE_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain);iterator.Valid();iterator.Next()){
            if(multigrid_poisson.cell_type(iterator.index)==MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE){
                const T& xi=x(iterator.index);
                x_min=min(x_min,xi);
                x_max=max(x_max,xi);
            }
        }
        x_as_density-=x_min;
        T x_range=x_max-x_min;
        if(x_range)
            x_as_density/=x_range;

        FILE_UTILITIES::Create_Directory(output_dir+"/"+f);
        FILE_UTILITIES::Write_To_File<float>(output_dir+"/"+f+"/density",x_as_density);
        FILE_UTILITIES::Write_To_Text_File(output_dir+"/common/last_frame",frame);
    }

    void Write_Substep(int substep,int frame)
    {
        LOG::SCOPE scope("Write substep");
        MULTIGRID_POISSON<T,d>& multigrid_poisson=multigrid_poisson_solver->Discretization();
        std::string f=LOG::sprintf("%i",frame);
        std::string s=LOG::sprintf("%i",substep);

        FILE_UTILITIES::Create_Directory(output_dir+"/Frame_"+f+"_x");        
        FILE_UTILITIES::Create_Directory(output_dir+"/Frame_"+f+"_residual");
        ARRAY<T,TV_INT> x_as_density(x);
        ARRAY<T,TV_INT> r_as_density(b);
        T x_min=std::numeric_limits<T>::max();
        T x_max=-x_min;
        T r_max=x_max;
        for(RANGE_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain);iterator.Valid();iterator.Next()){
            const T_INDEX& index=iterator.index;

            if(multigrid_poisson.cell_type(index)==MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE){
                const T& xi=x(index);
                x_min=min(x_min,xi);
                x_max=max(x_max,xi);
            }
            T& ri=r_as_density(index);
            ri=abs(ri-multigrid_poisson.Apply_System_Matrix(index,x))/sqr(multigrid_poisson.h);
            r_max=max(r_max,ri);
        }
        x_as_density-=x_min;
        T x_range=x_max-x_min;
        if(x_range)
            x_as_density/=x_range;
        if(r_max)
            r_as_density/=r_max;
        
        if(substep==0){
            FILE_UTILITIES::Write_To_File<float>(output_dir+"/Frame_"+f+"_x/grid",multigrid_poisson.grid);
            FILE_UTILITIES::Write_To_File<float>(output_dir+"/Frame_"+f+"_residual/grid",multigrid_poisson.grid);
        }
        FILE_UTILITIES::Write_To_File<float>(output_dir+"/Frame_"+f+"_x/density."+s,x_as_density);
        FILE_UTILITIES::Write_To_File<float>(output_dir+"/Frame_"+f+"_residual/density."+s,r_as_density);

        FILE_UTILITIES::Write_To_Text_File(output_dir+"/Frame_"+f+"_x/last_frame",s);
        FILE_UTILITIES::Write_To_Text_File(output_dir+"/Frame_"+f+"_residual/last_frame",s);


    }

    void Set_Cell_Type(int frame)
    {
        LOG::SCOPE scope("Set cell type");
        MULTIGRID_POISSON<T,d>& multigrid_poisson=multigrid_poisson_solver->Discretization();
        switch(test_number){
            case 1:
            case 2:                
                if(frame>0) return;
                for(RANGE_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain);iterator.Valid();iterator.Next())
                    multigrid_poisson.cell_type(iterator.index)=MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE;

                for(BOUNDARY_ITERATOR<d> iterator(multigrid_poisson.padded_domain);iterator.Valid();iterator.Next())
                    multigrid_poisson.cell_type(iterator.Index())=MULTIGRID_POISSON<T,d>::NEUMANN_CELL_TYPE;
                break;
            case 3:
            case 4:{
                if(frame>0) return;
                T radius=.125;

                TV circle_center=TV::All_Ones_Vector()*.5;
                circle_center(2)=2*radius;
                T radius_squared=radius*radius;

                for(BOUNDARY_ITERATOR<d> iterator(multigrid_poisson.padded_domain);iterator.Valid();iterator.Next())
                    multigrid_poisson.cell_type(iterator.Index())=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;

                for(RANGE_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain);iterator.Valid();iterator.Next()){
                    TV X=multigrid_poisson.grid.Node(iterator.index);
                    if((X-circle_center).Magnitude_Squared()<radius_squared)
                        multigrid_poisson.cell_type(iterator.index)=MULTIGRID_POISSON<T,d>::NEUMANN_CELL_TYPE;
                    else
                        multigrid_poisson.cell_type(iterator.index)=MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE;
                }
            }break;
            case 5:
            case 6:

                for(BOUNDARY_ITERATOR<d> iterator(multigrid_poisson.padded_domain);iterator.Valid();iterator.Next())
                    multigrid_poisson.cell_type(iterator.Index())=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;

                TV center;
                center(1)=(T).3*sin((T)frame*.1)+.5;
                center(2)=(T).3*cos((T)frame*.1)+.5;
                if(d==3)
                    center(3)=(T).3*cos((T)frame*.03)+.5;
                T radius=.125;
                T radius_squared=radius*radius;
                for(RANGE_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain);iterator.Valid();iterator.Next()){
                    TV X=multigrid_poisson.grid.Node(iterator.index);
                    if((X-center).Magnitude_Squared()<radius_squared)
                        multigrid_poisson.cell_type(iterator.index)=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;
                    else
                        multigrid_poisson.cell_type(iterator.index)=MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE;
                }
                break;
                   
                    
        }

    }

    void Set_Boundary_Values(int frame)
    {
        LOG::SCOPE scope("Set boundary values");
        if(test_number==5 ||test_number==6){
            MULTIGRID_POISSON<T,d>& multigrid_poisson=multigrid_poisson_solver->Discretization();
            for(RANGE_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain.Thickened(-1));iterator.Valid();iterator.Next()){
                const T_INDEX& index=iterator.index;

                b(index)=0;
                if(multigrid_poisson.cell_type(index)==MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE){
                    for(int v=0;v<d;v++){
                        const T_INDEX& axis_vector=T_INDEX::Axis_Vector(v);
                        if(multigrid_poisson.cell_type(index+axis_vector)==MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE)
                            b(index)-=(T)1/sqr(multigrid_poisson.h);

                        if(multigrid_poisson.cell_type(index-axis_vector)==MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE)
                            b(index)-=(T)1/sqr(multigrid_poisson.h);
                    }
                        
                }
                if(multigrid_poisson.cell_type(index)==MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE)
                      x(index)=1;
            }    
        }
    }

};

}
#endif
