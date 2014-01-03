//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/INTERVAL.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Parallel_Computation/INT_ITERATOR_THREADED.h>
#include <Tools/Parallel_Computation/PCG_SPARSE_THREADED.h>
using namespace PhysBAM;
//#####################################################################
// Function Parallel_Solve
//#####################################################################
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,const T tolerance)
{
    Init_Barriers();
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    assert(tid==domain_index(interior_domain.max_corner));
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);

    int global_n=A.n,interior_n=interior_indices.Size();
    T global_tolerance=tolerance;
    int desired_iterations=global_n;if(enforce_compatibility) desired_iterations--;if(maximum_iterations) desired_iterations=min(desired_iterations,maximum_iterations);

    ARRAY<T> z_interior(interior_n);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&sum_lock);
#endif
    temp.Resize(global_n);p.Resize(global_n);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&sum_lock);
#endif

    // build interior views of x,b,p,z,temp
    ARRAY_VIEW<T> x_interior(x.Array_View(interior_indices)),b_interior(b.Array_View(interior_indices)),p_interior(p.Array_View(interior_indices)),temp_interior(temp.Array_View(interior_indices));

    // adjust x for the null space
    if(enforce_compatibility && remove_null_space_solution_component) x_interior-=(T)(Global_Sum((T)x_interior.Sum_Double_Precision(),tid)/global_n);

    // find initial residual, r=b-Ax - reusing b for the residual
#ifdef USE_PTHREADS
    pthread_barrier_wait(&barr);
#endif
    A.Times(x,temp);
    b_interior-=temp_interior;
    if(enforce_compatibility) b_interior-=(T)(Global_Sum((T)b_interior.Sum_Double_Precision(),tid)/global_n);
    if(Global_Max(b_interior.Max_Abs())<=global_tolerance){
        if(show_results) LOG::cout<<"NO ITERATIONS NEEDED"<<std::endl;
        return;}

    SPARSE_MATRIX_FLAT_MXN<T>* C=0;
    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    if(incomplete_cholesky){
        C=A.Create_Submatrix(interior_indices);
        C->In_Place_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
            preconditioner_zero_tolerance,preconditioner_zero_replacement);}

    double rho=0,rho_old=0;
    for(int iteration=0;;iteration++){
        if(incomplete_cholesky){
            // solve Mz=r
            C->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
            C->Solve_Backward_Substitution(temp_interior,z_interior,false,true);} // diagonal is inverted to save on divides
        else z_interior=b_interior; // set z=r when there is no preconditioner

        // for Neumann boundary conditions only, make sure z sums to zero
        if(enforce_compatibility) z_interior-=(T)(Global_Sum((T)z_interior.Sum_Double_Precision(),tid)/global_n);

        // update search direction
        rho_old=rho;rho=Global_Sum((T)ARRAY<T>::Dot_Product_Double_Precision(z_interior,b_interior),tid);
        T beta=0;if(iteration==0) p_interior=z_interior;else{beta=(T)(rho/rho_old);for(int i=0;i<interior_n;i++) p_interior(i)=z_interior(i)+beta*p_interior(i);} // when iteration=1, beta=0

        // update solution and residual
#ifdef USE_PTHREADS
        pthread_barrier_wait(&barr);
#endif
        A.Times(p,temp);
        T alpha=(T)(rho/Global_Sum((T)p_interior.Dot_Product_Double_Precision(p_interior,temp_interior),tid));
        for(int i=0;i<interior_n;i++){x_interior(i)+=alpha*p_interior(i);b_interior(i)-=alpha*temp_interior(i);}

        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
        if(enforce_compatibility) b_interior-=(T)(Global_Sum((T)b_interior.Sum_Double_Precision(),tid)/global_n);

        T residual=Global_Max(b_interior.Max_Abs());

        // check for convergence
        if(show_residual) LOG::cout<<residual<<std::endl;
        if(residual<=global_tolerance){if(show_results) LOG::cout<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations-1){if(show_results) LOG::cout<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
    }
    delete C;
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_In_Parts(SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,const T tolerance)
{
    int global_n=A.n;
    T global_tolerance=tolerance;
    int desired_iterations=global_n;if(enforce_compatibility) desired_iterations--;if(maximum_iterations) desired_iterations=min(desired_iterations,maximum_iterations);
    ARRAY<T> z(global_n,false);
    temp.Resize(global_n);p.Resize(global_n);
    
    INT_ITERATOR_THREADED_ALPHA<PCG_SPARSE_THREADED<TV> > threaded_iterator(1,global_n,&thread_queue);
    int num_intervals=threaded_iterator.intervals.m;
    ARRAY<T> local_sum(num_intervals);
    
    // adjust x for the null space
    T sum=0;
    if(enforce_compatibility && remove_null_space_solution_component){
        threaded_iterator.template Run<ARRAY<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Sum,x,local_sum);
        for(int i=0;i<num_intervals;i++) sum+=local_sum(i);}
    
    // find initial residual, r=b-Ax - reusing b for the residual
    threaded_iterator.template Run<SPARSE_MATRIX_FLAT_MXN<T>&,ARRAY<T>&,ARRAY<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Part_One,A,x,b,sum/global_n);
    if(enforce_compatibility){sum=0;
        threaded_iterator.template Run<ARRAY<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Sum,b,local_sum);
        for(int i=0;i<num_intervals;i++) sum+=local_sum(i);
        threaded_iterator.template Run<ARRAY<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Subtract,b,sum/global_n);}

    threaded_iterator.template Run<ARRAY<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Max,b,local_sum);
    T max_val=0;for(int i=0;i<num_intervals;i++) max_val=max(max_val,local_sum(i));
    if(max_val<=global_tolerance){
        if(show_results) LOG::cout<<"NO ITERATIONS NEEDED"<<std::endl;
        return;}

    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    SPARSE_MATRIX_FLAT_MXN<T>* C=0;
    if(incomplete_cholesky){
        C=new SPARSE_MATRIX_FLAT_MXN<T>(A);
        C->In_Place_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
            preconditioner_zero_tolerance,preconditioner_zero_replacement);}

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
       if(incomplete_cholesky){
            // solve Mz=r
            C->Solve_Forward_Substitution(b,temp,true); // diagonal should be treated as the identity
            C->Solve_Backward_Substitution(temp,z,false,true);} // diagonal is inverted to save on divides
        else z=b; // set z=r when there is no preconditioner

        // for Neumann boundary conditions only, make sure z sums to zero
        if(enforce_compatibility){T sum=0;
            threaded_iterator.template Run<ARRAY<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Sum,z,local_sum);
            for(int i=0;i<num_intervals;i++) sum+=local_sum(i);
            threaded_iterator.template Run<ARRAY<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Subtract,z,sum/global_n);}
        
        // update search direction
        threaded_iterator.template Run<ARRAY<T>&,ARRAY<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Dot,z,b,local_sum);
        rho_old=rho;rho=0;for(int i=0;i<num_intervals;i++) rho+=local_sum(i);
        threaded_iterator.template Run<ARRAY<T>&,T,T,int>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Part_Two,z,(T)rho,(T)rho_old,iteration);
        
        // update solution and residual
        threaded_iterator.template Run<SPARSE_MATRIX_FLAT_MXN<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Part_Three,A);
        threaded_iterator.template Run<ARRAY<T>&,ARRAY<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Dot,p,temp,local_sum);
        T sum=0;for(int i=0;i<num_intervals;i++) sum+=local_sum(i);
        threaded_iterator.template Run<ARRAY<T>&,ARRAY<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Part_Four,x,b,(T)(rho/sum));

        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
        if(enforce_compatibility){T sum=0;
            threaded_iterator.template Run<ARRAY<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Sum,b,local_sum);
            for(int i=0;i<num_intervals;i++) sum+=local_sum(i);
            threaded_iterator.template Run<ARRAY<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Subtract,b,sum/global_n);}
        
        threaded_iterator.template Run<ARRAY<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Max,b,local_sum);
        T residual=0;for(int i=0;i<num_intervals;i++) residual=max(residual,local_sum(i));

        // check for convergence
        if(show_residual) LOG::cout<<residual<<std::endl;
        if(residual<=global_tolerance){if(show_results) LOG::cout<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations){if(show_results) LOG::cout<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
    }
    delete C;
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_In_Parts(DOMAIN_ITERATOR_THREADED_ALPHA<PCG_SPARSE_THREADED<TV>,TV>& threaded_iterator,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,const T tolerance)
{
    int global_n=A.n,num_domains=threaded_iterator.domains.m;
    T global_tolerance=tolerance;
    int desired_iterations=global_n;if(enforce_compatibility) desired_iterations--;if(maximum_iterations) desired_iterations=min(desired_iterations,maximum_iterations);
    ARRAY<ARRAY_VIEW<T> > z_interior(num_domains),x_interior(num_domains),b_interior(num_domains),p_interior(num_domains),temp_interior(num_domains);
    temp.Resize(global_n);
    p.Resize(global_n);
    z.Resize(global_n);
    ARRAY<T> local_sum(num_domains);
    
    // build interior views of x,b,p,z,temp
    threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<T>&,ARRAY<T>&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_One,domain_index,all_interior_indices,x,b,z_interior,x_interior,b_interior,p_interior,temp_interior);
    // adjust x for the null space
    if(enforce_compatibility && remove_null_space_solution_component){T sum=0;
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Sum,domain_index,all_interior_indices,x_interior,local_sum);
        for(int i=0;i<num_domains;i++) sum+=local_sum(i);
        //threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Distribute,domain_index,all_interior_indices,x_interior,sum/global_n);}
        for(int i=0;i<num_domains;i++) x_interior(i)-=sum/global_n;}
    
    // find initial residual, r=b-Ax - reusing b for the residual
    threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,const ARRAY<ARRAY<INTERVAL<int> > >&,SPARSE_MATRIX_FLAT_MXN<T>&,ARRAY<T>&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Two,domain_index,all_interior_indices,all_ghost_indices,A,x,b_interior,temp_interior);
    if(enforce_compatibility){T sum=0;
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Sum,domain_index,all_interior_indices,b_interior,local_sum);
        for(int i=0;i<num_domains;i++) sum+=local_sum(i);
        //threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Distribute,domain_index,all_interior_indices,b_interior,sum/global_n);}
        for(int i=0;i<num_domains;i++) b_interior(i)-=sum/global_n;}

    threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Max,domain_index,all_interior_indices,b_interior,local_sum);
    T max_val=0;for(int i=0;i<num_domains;i++) max_val=max(max_val,local_sum(i));
    if(max_val<=global_tolerance){
        if(show_results) LOG::cout<<"NO ITERATIONS NEEDED"<<std::endl;
        return;}

    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T>*> C(num_domains);
    threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,SPARSE_MATRIX_FLAT_MXN<T>&,ARRAY<SPARSE_MATRIX_FLAT_MXN<T>*>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Three,domain_index,all_interior_indices,A,C);

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<SPARSE_MATRIX_FLAT_MXN<T>*>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Four,domain_index,all_interior_indices,z_interior,b_interior,temp_interior,C);
        
        // for Neumann boundary conditions only, make sure z sums to zero
        if(enforce_compatibility){T sum=0;
            threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Sum,domain_index,all_interior_indices,z_interior,local_sum);
            for(int i=0;i<num_domains;i++) sum+=local_sum(i);
            //threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Distribute,domain_index,all_interior_indices,z_interior,sum/global_n);}
            for(int i=0;i<num_domains;i++) z_interior(i)-=sum/global_n;}
        
        // update search direction
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Dot,domain_index,all_interior_indices,z_interior,b_interior,local_sum);
        rho_old=rho;rho=0;for(int i=0;i<num_domains;i++) rho+=local_sum(i);
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,T,T,int>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Five,domain_index,all_interior_indices,z_interior,p_interior,(T)rho,(T)rho_old,iteration);
        
        // update solution and residual
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,const ARRAY<ARRAY<INTERVAL<int> > >&,const SPARSE_MATRIX_FLAT_MXN<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Six,domain_index,all_interior_indices,all_ghost_indices,A);
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Dot,domain_index,all_interior_indices,p_interior,temp_interior,local_sum);
        T sum=0;for(int i=0;i<num_domains;i++) sum+=local_sum(i);
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<ARRAY_VIEW<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Seven,domain_index,all_interior_indices,x_interior,b_interior,p_interior,temp_interior,(T)(rho/sum));

        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
        if(enforce_compatibility){T sum=0;
            threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Sum,domain_index,all_interior_indices,b_interior,local_sum);
            for(int i=0;i<num_domains;i++) sum+=local_sum(i);
            //threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Distribute,domain_index,all_interior_indices,b_interior,sum/global_n);}
            for(int i=0;i<num_domains;i++) b_interior(i)-=sum/global_n;}

        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<ARRAY_VIEW<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Max,domain_index,all_interior_indices,b_interior,local_sum);
        T residual=0;for(int i=0;i<num_domains;i++) residual=max(residual,local_sum(i));

        // check for convergence
        if(show_residual) LOG::cout<<residual<<std::endl;
        if(residual<=global_tolerance){if(show_results) LOG::cout<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations){if(show_results) LOG::cout<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
    }
    for(int i=0;i<num_domains;i++) delete C(i);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_One(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<T>& x,ARRAY<T>& b,ARRAY<ARRAY_VIEW<T> >& z_interior,ARRAY<ARRAY_VIEW<T> >& x_interior,ARRAY<ARRAY_VIEW<T> >& b_interior,ARRAY<ARRAY_VIEW<T> >& p_interior,ARRAY<ARRAY_VIEW<T> >& temp_interior)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    assert(tid==domain_index(interior_domain.max_corner));
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);

    z_interior(tid)=z.Array_View(interior_indices);
    x_interior(tid)=x.Array_View(interior_indices);
    b_interior(tid)=b.Array_View(interior_indices);
    p_interior(tid)=p.Array_View(interior_indices);
    temp_interior(tid)=temp.Array_View(interior_indices);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Two(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& x,ARRAY<ARRAY_VIEW<T> >& b_interior,ARRAY<ARRAY_VIEW<T> >& temp_interior)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
 
    A.Times(x,temp);
    b_interior(tid)-=temp_interior(tid);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Three(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<SPARSE_MATRIX_FLAT_MXN<T>*>& C)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    
    C(tid)=0;
    if(incomplete_cholesky){
        C(tid)=A.Create_Submatrix(interior_indices);
        C(tid)->In_Place_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
            preconditioner_zero_tolerance,preconditioner_zero_replacement);}
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Four(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<ARRAY_VIEW<T> >& z_interior,ARRAY<ARRAY_VIEW<T> >& b_interior,ARRAY<ARRAY_VIEW<T> >& temp_interior,ARRAY<SPARSE_MATRIX_FLAT_MXN<T>*>& C)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    
    if(incomplete_cholesky){
        // solve Mz=r
        C(tid)->Solve_Forward_Substitution(b_interior(tid),temp_interior(tid),true); // diagonal should be treated as the identity
        C(tid)->Solve_Backward_Substitution(temp_interior(tid),z_interior(tid),false,true);} // diagonal is inverted to save on divides
    else z_interior(tid)=b_interior(tid); // set z=r when there is no preconditioner
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Five(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<ARRAY_VIEW<T> >& z_interior,ARRAY<ARRAY_VIEW<T> >& p_interior,T rho,T rho_old,int iteration)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    int interior_n=interior_indices.Size();
    
    T beta=0;if(iteration==0) p_interior(tid)=z_interior(tid);else{beta=(T)(rho/rho_old);for(int i=0;i<interior_n;i++) p_interior(tid)(i)=z_interior(tid)(i)+beta*p_interior(tid)(i);}
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Six(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,const SPARSE_MATRIX_FLAT_MXN<T>& A)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();

    A.Times(p,temp);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Seven(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<ARRAY_VIEW<T> >& x_interior,ARRAY<ARRAY_VIEW<T> >& b_interior,ARRAY<ARRAY_VIEW<T> >& p_interior,ARRAY<ARRAY_VIEW<T> >& temp_interior,T alpha)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    int interior_n=interior_indices.Size();
   
    for(int i=0;i<interior_n;i++){x_interior(tid)(i)+=alpha*p_interior(tid)(i);b_interior(tid)(i)-=alpha*temp_interior(tid)(i);}
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Distribute(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<ARRAY_VIEW<T> >& interior,const T sum)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    interior(tid)-=sum;
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Sum(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<ARRAY_VIEW<T> >& interior,ARRAY<T>& sum)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    sum(tid)=(T)interior(tid).Sum_Double_Precision();
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Max(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<ARRAY_VIEW<T> >& interior,ARRAY<T>& sum)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    sum(tid)=interior(tid).Max_Abs();
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Dot(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<ARRAY_VIEW<T> >& interior_1,ARRAY<ARRAY_VIEW<T> >& interior_2,ARRAY<T>& sum)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    sum(tid)=(T)interior_1(tid).Dot_Double_Precision(interior_2(tid));
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Subtract(ARRAY<T>& vector,const T sum,int start_index,int end_index)
{
    for(int i=start_index;i<end_index;i++) vector(i)-=sum;
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Sum(ARRAY<T>& vector,ARRAY<T>& sum,int start_index,int end_index,int tid)
{
    sum(tid)=(T)vector.Array_View(start_index,end_index-start_index).Sum_Double_Precision();
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Dot(ARRAY<T>& vector1,ARRAY<T>& vector2,ARRAY<T>& sum,int start_index,int end_index,int tid)
{
    sum(tid)=(T)vector1.Array_View(start_index,end_index-start_index).Dot_Double_Precision(vector2.Array_View(start_index,end_index-start_index));
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Max(ARRAY<T>& vector,ARRAY<T>& sum,int start_index,int end_index,int tid)
{
    sum(tid)=vector.Array_View(start_index,end_index-start_index).Max_Abs();
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Part_One(SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,const T sum,int start_index,int end_index)
{
    for(int i=start_index;i<end_index;i++) x(i)-=sum;
    A.Times(p,temp);
    for(int i=start_index;i<end_index;i++) b(i)-=temp(i);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Part_Two(ARRAY<T>& z,T rho,T rho_old,int iteration,int start_index,int end_index)
{
    T beta=0;if(iteration==0){for(int i=start_index;i<end_index;i++) p(i)=z(i);}else{beta=(T)(rho/rho_old);for(int i=start_index;i<end_index;i++) p(i)=z(i)+beta*p(i);}
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Part_Three(SPARSE_MATRIX_FLAT_MXN<T>& A,int start_index,int end_index)
{
    A.Times(p,temp);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Part_Four(ARRAY<T>& x,ARRAY<T>& b,T alpha,int start_index,int end_index)
{
    for(int i=start_index;i<end_index;i++){x(i)+=alpha*p(i);b(i)-=alpha*temp(i);}
}
//#####################################################################
namespace PhysBAM{
template class PCG_SPARSE_THREADED<VECTOR<float,1> >;
template class PCG_SPARSE_THREADED<VECTOR<float,2> >;
template class PCG_SPARSE_THREADED<VECTOR<float,3> >;
template class PCG_SPARSE_THREADED<VECTOR<double,1> >;
template class PCG_SPARSE_THREADED<VECTOR<double,2> >;
template class PCG_SPARSE_THREADED<VECTOR<double,3> >;
}
