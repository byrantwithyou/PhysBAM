//#####################################################################
// Copyright 2009, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.
//#####################################################################

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>

#include "../multigrid_poisson_3d_optimized_kernels/Thread_Queueing/PTHREAD_QUEUE.h"
#include "MG_PRECONDITIONED_CONJUGATE_GRADIENT.h"
#include "MULTIGRID_POISSON.h"
#include "MULTIGRID_POISSON_REFINEMENT.h"
#include "MULTIGRID_POISSON_SOLVER.h"


using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
int main(int argc,char* argv[])
{ 
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    static const int d=2;

    if(argc!=4){
        std::cout<<"Usage: "<<argv[0]<<"<test_number> <number_of_threads> <ic or mg>"<<std::endl; return 1;
    }

    const int test_number=atoi(argv[1]);
    const int number_of_threads=atoi(argv[2]);
    std::string preconditioner=argv[3];

    if(preconditioner!="mg" && preconditioner!="ic"){
        std::cout<<"4th argument must be ic or mg"<<std::endl; return 1;
    }
        

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;

    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);


    int m,n;
    int levels;
    std::string output_dir;

 
    pthread_queue=new PhysBAM::PTHREAD_QUEUE(number_of_threads);

    switch(test_number){
        case 1: // Water
            m=512;
            n=512;
            levels=7;
            output_dir="Water_512";
            break;
        case 2: // car
            m=768;
            n=768;
            levels=6;
            output_dir="Car_768";
            break;
        case 3: // sphere 64
            m=64;
            n=64;
            levels=4;
            output_dir="sphere_64";
            break;
        case 4:
            m=96;
            n=96;
            levels=3;
            output_dir="sphere_96";
            break;
        case 5:
            m=128;
            n=128;
            levels=5;
            output_dir="sphere_128";
            break;
        case 6:
            m=192;
            n=192;
            levels=4;
            output_dir="sphere_192";
            break;
        case 7:
            m=256;
            n=256;
            levels=6;
            output_dir="sphere_256";
            break;
        case 8:
            m=384;
            n=384;
            levels=5;
            output_dir="sphere_384";
            break;
        case 9:
            m=512;
            n=512;
            levels=7;
            output_dir="sphere_512";
            break;
        case 10:
            m=768;
            n=768;
            levels=6;
            output_dir="sphere_768";
            break;
        case 11:
            m=768;
            n=1152;
            levels=6;
            output_dir="sphere_768_1152";
            break;
    }

    T_INDEX size=m*T_INDEX::All_Ones_Vector();
    size(2)=n;
    output_dir=output_dir+"_"+preconditioner;

    const T h=(T)1/(T)m;                

    FILE_UTILITIES::Create_Directory(output_dir);

    LOG::Initialize_Logging();
    LOG::Instance()->Copy_Log_To_File(output_dir+"/log.txt",false);
    LOG::cout<<"m: "<<m<<" n : "<<n<<" levels : "<<levels<<std::endl;

    if(preconditioner=="mg"){        
        
        LOG::Time("MULTIGRID_POISSON_SOLVER::Construction");
        MULTIGRID_POISSON_SOLVER<T,d> multigrid_poisson_solver(size,h,levels,number_of_threads);
        MULTIGRID_POISSON<T,d>& multigrid_poisson=multigrid_poisson_solver.Discretization();
        
        switch(test_number){
            case 1: // Water
                FILE_UTILITIES::Read_From_File<T>("/u/gojira/s1/terangroup/amcadams/free_surface_output/medium_pour_512_deep_pool/cell_type.90",multigrid_poisson.cell_type);
                break;
            case 2: // car
                FILE_UTILITIES::Read_From_File<T>("/u/godzilla/s2/jteran/aleka/smoke_output/768_car/cell_type.10",multigrid_poisson.cell_type);
                break;
            case 3: // sphere 64
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
            case 9:
            case 10:
            case 11:
                multigrid_poisson.Initialize_Square_Minus_Circle_Domain();
                break;
        }
        
        multigrid_poisson.Initialize_Test_Right_Hand_Side(multigrid_poisson.b);
        
        
        {            
            LOG::SCOPE scope("Initialize_Multigrid_Hierarchy");
            multigrid_poisson_solver.Initialize_Multigrid_Hierarchy();
        }
        
        
        
#if 1 // no visualization

        {
            typedef ARRAY<T,VECTOR<int,d> > T_VARIABLE;
            T_VARIABLE x(multigrid_poisson.grid.Domain_Indices());
            T_VARIABLE p(multigrid_poisson.grid.Domain_Indices());
            
            MULTIGRID_SYSTEM<T,d> multigrid_system(multigrid_poisson_solver);
            MG_PRECONDITIONED_CONJUGATE_GRADIENT<T,d> cg;
            
            cg.print_residuals=true;
            
            cg.Solve(multigrid_system,
                x,
                /*b*/multigrid_poisson.b,
                /*q*/multigrid_poisson.u,
                p,
                1e-8,1,1000);
            
            
        }

#else

        {
            ARRAY<T,TV_INT> x(multigrid_poisson.grid);
            ARRAY<T,TV_INT> p(multigrid_poisson.grid);
            
            ARRAY<T,TV_INT> x_save(multigrid_poisson.grid); //for visualization only
            ARRAY<T,TV_INT> b_save(multigrid_poisson.grid); //for visualization only
            b_save=multigrid_poisson.b;
            x_save=x;
            
            MULTIGRID_SYSTEM<T,d> multigrid_system(multigrid_poisson_solver);
            MG_PRECONDITIONED_CONJUGATE_GRADIENT<T,d> cg;
            
            cg.print_residuals=true;
            
            FILE_UTILITIES::Create_Directory("output_x");
            FILE_UTILITIES::Create_Directory("output_r");
            FILE_UTILITIES::Create_Directory("output_rhs");
            FILE_UTILITIES::Write_To_File<float>("output_x/grid",multigrid_poisson.grid);
            FILE_UTILITIES::Write_To_File<float>("output_r/grid",multigrid_poisson.grid);
            FILE_UTILITIES::Write_To_File<float>("output_rhs/grid",multigrid_poisson.grid);
            ARRAYS_3D<T> r_as_density(multigrid_poisson.grid);
            ARRAYS_3D<T> rhs_as_density(multigrid_poisson.grid);
            
            T rhs_min=b_save.Min();
            T rhs_max=b_save.Max();
            
            rhs_as_density=b_save;
            rhs_as_density-=rhs_min;
            T rhs_range=(rhs_max-rhs_min);
            if(rhs_range)
                rhs_as_density/=rhs_range;
            
            LOG::cout<<"Max rhs : "<<rhs_max<<std::endl;
            
            FILE_UTILITIES::Write_To_File<RW>("output_rhs/density.0",rhs_as_density);
            
            for(int i=0;i<30;i++){
                x=x_save;
                multigrid_poisson.b=b_save;
                cg.Solve(multigrid_system,
                    x,
                    /*b*/multigrid_poisson.b,
                    /*q*/multigrid_poisson.u,
                    p,
                    ///*r*/multigrid_poisson.b,
                    ///*z*/multigrid_poisson.u,
                    1e-10,1,i);
                
                T x_min=x.Max();
                T x_max=x.Min();
                for(NODE_ITERATOR<TV> iterator(multigrid_poisson.grid);iterator.Valid();iterator.Next())
                    if(multigrid_poisson.cell_type(iterator.Node_Index())==MULTIGRID_POISSON<T,3>::INTERIOR_CELL_TYPE){
                        x_min=min(x_min,x(iterator.Node_Index()));
                        x_max=max(x_max,x(iterator.Node_Index()));
                        r_as_density(iterator.Node_Index())=abs(b_save(iterator.Node_Index())-multigrid_poisson.Apply_System_Matrix(iterator.Node_Index(),x)/(h*h));
                    }
                T r_max=r_as_density.Max();
                LOG::cout<<"Max residual : "<<r_max<<std::endl;
                LOG::cout<<"Max x : "<<x_max<<std::endl;
                if(r_max)
                    r_as_density/=r_max;
                x-=x_min;
                T x_range=x_max-x_min;
                LOG::cout<<"x range : "<<x_range<<std::endl;
                if(x_range)
                    x/=x_range;
                
                FILE_UTILITIES::Write_To_File<RW>("output_x/density."+FILE_UTILITIES::Number_To_String(i),x);
                FILE_UTILITIES::Write_To_File<RW>("output_r/density."+FILE_UTILITIES::Number_To_String(i),r_as_density);
                FILE_UTILITIES::Write_To_Text_File("output_x/last_frame",i,"\n");
                FILE_UTILITIES::Write_To_Text_File("output_r/last_frame",i,"\n");
            }
            
        }

#endif
    }
#if 0
    else if(preconditioner=="ic"){

        if(test_number==1||test_number==2){
            LOG::cout<<"Test number not implemented for icpcg"<<std::endl; return 1;
        }

        typedef POLICY_UNIFORM<VECTOR<MULTIGRID_POISSON<T,d>::CELL_TYPE,d> >::ARRAYS_SCALAR T_CELL_TYPE_FIELD;
        typedef POLICY_UNIFORM<T_INDEX>::ARRAYS_SCALAR T_INT_ARRAY;
        GRID<TV> grid(size+2,RANGE<TV>(-TV::All_Ones_Vector()*.5*h,size*h+TV::All_Ones_Vector()*.5*h));
        T_CELL_TYPE_FIELD cell_type(grid);
        PCG_SPARSE<T> icpcg;
        ARRAY<T> b_icpcg;
        ARRAY<T> x_icpcg;
        SPARSE_MATRIX_FLAT_MXN<T> icpcg_matrix;
        LIST_ARRAY<T_INDEX> interior_indices;
        T_INT_ARRAY index_ids(grid);

        icpcg.Enforce_Compatibility(false);
        icpcg.Remove_Null_Space_Solution_Component(false);
        icpcg.Show_Residuals(true);
        icpcg.Show_Results(true);
        icpcg.Use_Modified_Incomplete_Cholesky(.99);
        
        index_ids.Resize(grid);index_ids.Fill(0);
        
        T radius=.125;
        TV circle_center(.5,radius*(T)2,.5);
        T radius_squared=radius*radius;
        LOG::cout<<"M"<<std::endl;
        // Set cell types
        for(NODE_ITERATOR<TV> iterator(grid,0,GRID<TV>::INTERIOR_REGION);iterator.Valid();iterator.Next()){
            const T_INDEX& index=iterator.Node_Index();
            if((iterator.Location()-circle_center).Magnitude_Squared()<radius_squared)
                cell_type(index)=MULTIGRID_POISSON<T,d>::NEUMANN_CELL_TYPE;
            else
                cell_type(index)=MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE;
        }
        
        LOG::cout<<"M"<<std::endl;
        for(NODE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next())
            cell_type(iterator.Node_Index())=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;
        for(NODE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,1);iterator.Valid();iterator.Next())
            cell_type(iterator.Node_Index())=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;
        for(NODE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,2);iterator.Valid();iterator.Next())
            cell_type(iterator.Node_Index())=MULTIGRID_POISSON<T,d>::NEUMANN_CELL_TYPE;
        for(NODE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,3);iterator.Valid();iterator.Next())
            cell_type(iterator.Node_Index())=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;
        for(NODE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,4);iterator.Valid();iterator.Next())
            cell_type(iterator.Node_Index())=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;
        for(NODE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,5);iterator.Valid();iterator.Next())
            cell_type(iterator.Node_Index())=MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE;

        {LOG::SCOPE scope("ICPCG Initialize");
        // set up pcg matrix
        for(NODE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            const T_INDEX& index=iterator.Node_Index();
            if(cell_type(index)==MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE)
                index_ids(index)=interior_indices.Append(index);
        }

        ARRAY<int> row_lengths(interior_indices.m);
        for(int equation_id=0;equation_id<interior_indices.m;equation_id++){
            row_lengths(equation_id)=1; //diagonal_entry
            
            for(int axis=0;axis<d;axis++){ //off-diag
                if(cell_type(interior_indices(equation_id)+T_INDEX::Axis_Vector(axis))==MULTIGRID_POISSON<T,3>::INTERIOR_CELL_TYPE)
                    row_lengths(equation_id)++;
                if(cell_type(interior_indices(equation_id)-T_INDEX::Axis_Vector(axis))==MULTIGRID_POISSON<T,3>::INTERIOR_CELL_TYPE)
                    row_lengths(equation_id)++;
            }
        }
        
        icpcg_matrix.Set_Row_Lengths(row_lengths);
        for(int equation_id=0;equation_id<interior_indices.m;equation_id++){
            const T_INDEX& equation_index=interior_indices(equation_id);
            for(int axis=0;axis<3;axis++){
                const T_INDEX& neighbor_index=equation_index+T_INDEX::Axis_Vector(axis);
                const int& neighbor_id=index_ids(neighbor_index);
                if(cell_type(neighbor_index)!=MULTIGRID_POISSON<T,3>::NEUMANN_CELL_TYPE){
                    icpcg_matrix(equation_id,equation_id)+=(T)1/(h*h);
                    if(cell_type(neighbor_index)==MULTIGRID_POISSON<T,3>::INTERIOR_CELL_TYPE){
                        assert(neighbor_id>0 && neighbor_id<=interior_indices.m);
                        icpcg_matrix(equation_id,neighbor_id)=-(T)1/(h*h);}
                }
                const T_INDEX& neighbor_index1=equation_index-T_INDEX::Axis_Vector(axis);
                const int& neighbor_id1=index_ids(neighbor_index1);
                if(cell_type(neighbor_index1)!=MULTIGRID_POISSON<T,3>::NEUMANN_CELL_TYPE){
                    icpcg_matrix(equation_id,equation_id)+=(T)1/(h*h);
                    if(cell_type(neighbor_index1)==MULTIGRID_POISSON<T,3>::INTERIOR_CELL_TYPE){
                        assert(neighbor_id1>0 && neighbor_id1<=interior_indices.m);
                        icpcg_matrix(equation_id,neighbor_id1)=-(T)1/(h*h);}
                }
                
            }
        }
        }
        // set right hand side
        b_icpcg.Resize(interior_indices.m);
        for(int equation_id=0;equation_id<interior_indices.m;equation_id++){
            TV X=grid.X(interior_indices(equation_id));
            b_icpcg(equation_id)=-sin(X(1)*2*pi)*cos(X(2)*2*pi)*sin(X(3)*2*pi);
        }


        //SOLVE!

        x_icpcg.Resize(interior_indices.m);
        x_icpcg.Fill(0);
        icpcg.Set_Maximum_Iterations(1000);

        {LOG::SCOPE scope("ICPCG Solve");
            icpcg.Solve(icpcg_matrix,x_icpcg,b_icpcg,1e-8);}
        

    }
#endif
    LOG::Finish_Logging();
    return 0;
}
