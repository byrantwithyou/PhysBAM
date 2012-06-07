//#####################################################################
// Copyright 2010, Elliot English, Jon Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Standard_Tests/STANDARD_TESTS.h"

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>

using namespace PhysBAM;

template<class T>
VECTOR_ND<T> operator-(const VECTOR_ND<T>& v0,const VECTOR_ND<T>& v1)
{
    assert(v0.Size()==v1.Size());
    VECTOR_ND<T> r(v0.Size());
    for(int i=0;i<r.Size();i++)
        r(i)=v0(i)-v1(i);
    return r;
}

template<class T>
SPARSE_MATRIX_FLAT_MXN<T> Build_Smoothed_Aggregation_Interpolation(const SPARSE_MATRIX_FLAT_MXN<T>& A)
{
    T strong_threshold=0.0;

    int n=A.m;
    /*ARRAY<int> coarse_set(n);
    coarse_set.Fill(-1);
    ARRAY<T> row_magnitudes(n);
    for(int i=0;i<n;i++)
    {
        T row_magnitude=0;
        for(int j=A.offsets(i);j<A.offsets(i+1);j++)
        {
            const SPARSE_MATRIX_ENTRY<T>& entry=A.A(j);
            row_magnitude+=entry.a*entry.a;
        }
        row_magnitude=sqrt(row_magnitude);

        if(coarse_set(i)==-1)
        {
            coarse_set(i)=1;
            for(int j=A.offsets(i);j<A.offsets(i+1);j++)
            {
                const SPARSE_MATRIX_ENTRY<T>& entry=A.A(j);
                if(entry.j!=i && coarse_set(entry.j)==-1 && fabs(entry.a)>(row_magnitude*strong_threshold))
                    coarse_set(entry.j)=0;
            }
        }
    }

    std::cout << "coarse_set = " << coarse_set << std::endl;*/

    int n_aggregates=0;
    ARRAY<ARRAY<int> > aggregates;

    ARRAY<int> parent_set(n);
    parent_set.Fill(-1);
    for(int i=0;i<n;i++)
    {
        T row_magnitude=0;
        for(int j=A.offsets(i);j<A.offsets(i+1);j++)
        {
            const SPARSE_MATRIX_ENTRY<T>& entry=A.A(j);
            row_magnitude+=entry.a*entry.a;
        }
        row_magnitude=sqrt(row_magnitude);

        if(parent_set(i)==-1)
        {
            aggregates.Append(ARRAY<int>());
            n_aggregates++;
            aggregates(n_aggregates).Append(i);
            parent_set(i)=n_aggregates;
            for(int j=A.offsets(i);j<A.offsets(i+1);j++)
            {
                const SPARSE_MATRIX_ENTRY<T>& entry=A.A(j);
                if(entry.j!=i && parent_set(entry.j)==-1 && fabs(entry.a)>(row_magnitude*strong_threshold))
                {
                    parent_set(entry.j)=n_aggregates;
                    aggregates(n_aggregates).Append(j);
                }
            }
        }
    }

    //std::cout << "parent_set = " << parent_set << std::endl;

    SPARSE_MATRIX_FLAT_MXN<T> aggregator;
    ARRAY<int> aggregator_row_counts(n);
    aggregator_row_counts.Fill(1);
    aggregator.Set_Row_Lengths(aggregator_row_counts);
    aggregator.n=n_aggregates;

    for(int i=0;i<n;i++)
        aggregator.Set_Element(i,parent_set(i),1.0);

    //std::cout << "aggregator " << std::endl << aggregator << std::endl;

    T eigenvalue_maximum=4;

    SPARSE_MATRIX_FLAT_MXN<T> smoother;
    smoother=A;
    smoother*=(-4.0/(3.0*eigenvalue_maximum));
    for(int i=0;i<n;i++)
        smoother.Add_Element(i,i,1.0);

    //std::cout << "smoother " << std::endl << smoother << std::endl;

    SPARSE_MATRIX_FLAT_MXN<T> interpolation=smoother*aggregator;

    //std::cout << "interpolation " << std::endl << interpolation << std::endl;

    for(int i=0;i<n;i++)
    {
        T diagonal;
        for(int j=A.offsets(i);j<A.offsets(i+1);j++)
            if(A.A(j).j==i)
                diagonal=A.A(j).a;
        T scaling=1.0/diagonal;
        for(int j=interpolation.offsets(i);j<interpolation.offsets(i);j++)
            interpolation.A(j).a*=scaling;
    }
    
    return interpolation;
}

template<class T>
void Gauss_Seidel(const SPARSE_MATRIX_FLAT_MXN<T>& A,const VECTOR_ND<T>& b,VECTOR_ND<T>& x,int iterations=1)
{
    int n=A.m;
    for(int i=0;i<iterations;i++)
        for(int j=0;j<n;j++)
        {
            T row_sum=0;
            T diagonal=1;
            for(int k=A.offsets(j);k<A.offsets(j+1);k++)
            {
                const SPARSE_MATRIX_ENTRY<T>& entry=A.A(k);
                if(entry.j==j)
                    diagonal=entry.a;
                else
                    row_sum+=entry.a*x(entry.j);
            }
            x(j)=(b(j)-row_sum)/diagonal;
            //std::cout << "row " << j << " " << row_sum << " " << diagonal << std::endl;
        }
}

template<class T>
VECTOR_ND<T> Residual(const SPARSE_MATRIX_FLAT_MXN<T>& A,const VECTOR_ND<T>& b,const VECTOR_ND<T>& x)
{
    VECTOR_ND<T> Ax(A.m);
    A.Times(x,Ax);
    VECTOR_ND<T> r=Ax-b;
    return r;
}

template<class T>
void Solve_Linear_System(const SPARSE_MATRIX_FLAT_MXN<T>& A,const VECTOR_ND<T>& b,VECTOR_ND<T>& x)
{
    LOG::cout << "Solver_Linear_System " << A.m << std::endl;
    //LOG::cout << "A = " << A << std::endl;
    //sleep(1);

    if(A.m<5)
    {
        Gauss_Seidel(A,b,x,100);
    }
    else
    {
        //LOG::cout << "residual 0 = " << Residual(A,b,x).Maximum_Magnitude() << std::endl;
        //LOG::cout << "x = " << x << std::endl;
    Gauss_Seidel(A,b,x,10);

    //LOG::cout << "residual 1 = " << Residual(A,b,x).Maximum_Magnitude() << std::endl;

    SPARSE_MATRIX_FLAT_MXN<T> It,I=Build_Smoothed_Aggregation_Interpolation(A);
    I.Transpose(It);

    int n=I.m;
    int n2=I.n;

    SPARSE_MATRIX_FLAT_MXN<T> A2=It*A*I;

    VECTOR_ND<T> r=Residual(A,b,x);
    VECTOR_ND<T> e2(n2),r2(n2);
    It.Times(r,r2);
    e2.Fill(0);

    //Gauss_Seidel(A2,r2,e2,100);
    Solve_Linear_System(A2,r2,e2);

    //LOG::cout << "coarse_residual = " << Residual(A2,r2,e2).Maximum_Magnitude() << std::endl;

    VECTOR_ND<T> e(n);
    I.Times(e2,e);
    
    //x-=e;

    //LOG::cout << "residual 2 = " << Residual(A,b,x).Maximum_Magnitude() << std::endl;

    Gauss_Seidel(A,b,x,10);
    //LOG::cout << "residual 3 = " << Residual(A,b,x).Maximum_Magnitude() << std::endl;
    //LOG::cout << "x = " << x << std::endl;
    }
}

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,3> TV;

//#define MULTILEVEL
#ifdef MULTILEVEL
    LOG::cout << "smoothed aggregation test" << std::endl;

    SPARSE_MATRIX_FLAT_MXN<T> Mi,C;
        
    int dim=40;
    ARRAY<T> masses;
    ARRAY<TV> positions;
    ARRAY<VECTOR<int,2> > constraint_pairs;
    ARRAY<TV> constraint_normals;
    for(int i=0;i<dim;i++)
        for(int j=0;j<2;j++)
        {
            masses.Append(1);
            positions.Append(TV(i,j+(T)rand()/RAND_MAX,0));
        }
    
    int n=masses.Size();
    for(int i=0;i<n;i++)
        for(int j=i+1;j<=n;j++)
        {
            //LOG::cout << "pair distance " << i << " " << j << " " << (positions(i)-positions(j)).Magnitude() << std::endl;
            if((positions(i)-positions(j)).Magnitude()<1.2)
            {
                constraint_pairs.Append(VECTOR<int,2>(i,j));
                constraint_normals.Append((positions(i)-positions(j)).Normalized());
            }
        }
    
    LOG::cout << "positions " << positions << std::endl;
    LOG::cout << "constraint_pairs " << constraint_pairs << std::endl;
    LOG::cout << "constraint_normals " << constraint_normals << std::endl;
    
    ARRAY<int> Mi_row_counts(TV::dimension*n);
    Mi_row_counts.Fill(1);
    Mi.Set_Row_Lengths(Mi_row_counts);
    Mi.n=TV::dimension*n;
    for(int i=0;i<n;i++)
        for(int j=0;j<TV::dimension;j++)
        {
            int k=TV::dimension*(i-1)+j;
            Mi.Set_Element(k,k,1.0/masses(i));
        }
    
    ARRAY<int> C_row_counts(constraint_pairs.m);
    C_row_counts.Fill(TV::dimension*2);
    C.Set_Row_Lengths(C_row_counts);
    C.n=Mi.n;
    for(int i=0;i<constraint_pairs.m;i++)
    {
        for(int j=0;j<TV::dimension;j++)
        {
            C.Set_Element(i,(constraint_pairs(i)(1)-1)*TV::dimension+j,-constraint_normals(i)(j));
            C.Set_Element(i,(constraint_pairs(i)(2)-1)*TV::dimension+j,constraint_normals(i)(j));
        }
    }
    
    //std::cout << "Mi = " << std::endl << Mi << std::endl;
    //std::cout << "C = " << std::endl << C << std::endl;
    
    SPARSE_MATRIX_FLAT_MXN<T> Ct;
    C.Transpose(Ct);
    
    SPARSE_MATRIX_FLAT_MXN<T> A=C*Mi*Ct;
    
    //std::cout << "A = " << std::endl << A << std::endl;
    
    VECTOR_ND<T> b(A.m),x(A.m);
    
    b.Fill(1);
    x.Fill(0);
    
    LOG::cout << "before residual " << Residual(A,b,x).Maximum_Magnitude() << std::endl;
    
    Solve_Linear_System(A,b,x);
    Solve_Linear_System(A,b,x);
    Solve_Linear_System(A,b,x);
    
    LOG::cout << "after residual " << Residual(A,b,x).Maximum_Magnitude() << std::endl;
    
    /*SPARSE_MATRIX_FLAT_MXN<T> It,I=Build_Smoothed_Aggregation_Interpolation(A);
    I.Transpose(It);
      
    SPARSE_MATRIX_FLAT_MXN<T> A2=It*A*I;
    
    std::cout << "A2 " << std::endl << A2 << std::endl;
    
    SPARSE_MATRIX_FLAT_MXN<T> I2t,I2=Build_Smoothed_Aggregation_Interpolation(A2);
    I2.Transpose(I2t);
    
    SPARSE_MATRIX_FLAT_MXN<T> A3=I2t*A2*I2;
    
    std::cout << "I2 " << std::endl << I2 << std::endl;
    
    std::cout << "A3 " << std::endl << A3 << std::endl;*/
#else
    EXAMPLE<TV>* example=new STANDARD_TESTS<T>(stream_type);
    example->Parse(argc,argv);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* solid_fluid_example=dynamic_cast<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >*>(example);
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",solid_fluid_example->restart);
    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*solid_fluid_example);
    driver.Execute_Main_Program();
    delete example;
#endif

    return 0;
}
//#####################################################################
