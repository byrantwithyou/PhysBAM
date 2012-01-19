//#####################################################################
// Copyright 2003-2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_PROCESSOR
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include "LEVELSET_PROCESSOR.h"
using namespace PhysBAM;
//#####################################################################
// Function Smooth_By_Curvature
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Smooth_By_Curvature(const int iters,const T sigma)
{
    levelset.Set_Curvature_Motion(sigma);
    FACE_ARRAYS<VECTOR<T,3> > V(grid);

    T step=levelset.CFL(V);
    std::cout<<"Beginning smoothing iterations with step "<<step<<"..."<<std::endl;
    for(int i=0;i<iters;i++){
        std::cout<<".";
        levelset.Euler_Step(V,step);}
    if(verbose) std::cout<<" done."<<std::endl;
}
//#####################################################################
// Function Window_Sum_3D
//#####################################################################
template<class T> template<class T2> void LEVELSET_PROCESSOR<T>::
Window_Sum_3D(ARRAYS<VECTOR<T2,3> >& phi,const int xradius,const int yradius,const int zradius) 
{
    if(xradius<0 || yradius<0 || zradius<0){std::cerr<<"invalid kernel: "<<xradius<<" "<<yradius<<" "<<zradius<<std::endl;exit(1);}
    T2 *vector;
    //for phi
    int m=phi.m,m_start=phi.domain.min_corner.x,m_end=phi.domain.max_corner.x;
    int n=phi.n,n_start=phi.domain.min_corner.y,n_end=phi.domain.max_corner.y;
    int mn=phi.mn,mn_start=phi.domain.min_corner.z,mn_end=phi.domain.max_corner.z;
    std::cout<<"x-axis sums:"<<std::endl;
    vector = new T2[m];
    for(int j=n_start;j<=n_end;j++)for(int k=mn_start;k<=mn_end;k++){
        // Initial summation
        vector[0]=T2(0);
        for(int x=-xradius;x<=xradius;x++) vector[0]+=phi(clamp(m_start+x,m_start,m_end),j,k);            
        // Iterative update
        for(int x=1;x<m;x++) vector[x] = vector[x-1] -
                phi(clamp(m_start+x-xradius-1,m_start,m_end),j,k) +
                phi(clamp(m_start+x+xradius,m_start,m_end),j,k);
        // Copy back
        for(int x=0;x<m;x++) phi(m_start+x,j,k)=vector[x];}
    delete[] vector;
    //------------------------------------------------------------------
    std::cout<<"y-axis sums:"<<std::endl;
    vector = new T2[n];
    for(int i=m_start;i<=m_end;i++)for(int k=mn_start;k<=mn_end;k++){
        // Initial summation
        vector[0] = T2(0);
        for(int y=-yradius;y<=yradius;y++) vector[0]+=phi(i,clamp(n_start+y,n_start,n_end),k);
        // Iterative update
        for(int y=1;y<n;y++) vector[y] = vector[y-1] -
                    phi(i,clamp(n_start+y-yradius-1,n_start,n_end),k) +
                    phi(i,clamp(n_start+y+yradius,n_start,n_end),k);
        // Copy back
        for(int y=0;y<n;y++) phi(i,n_start+y,k)=vector[y];}
    delete[] vector;
    //------------------------------------------------------------------
    std::cout<<"z-axis sums:"<<std::endl;
    vector = new T2[mn];
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++){
        // Initial summation
        vector[0] = T2(0);
        for(int z=-zradius;z<=zradius;z++) vector[0]+=phi(i,j,clamp(mn_start+z,mn_start,mn_end));
        // Iterative update
        for(int z=1;z<mn;z++) vector[z] = vector[z-1] -
                    phi(i,j,clamp(mn_start+z-zradius-1,mn_start,mn_end)) +
                    phi(i,j,clamp(mn_start+z+zradius,mn_start,mn_end));
        // Copy back
        for(int z=0;z<mn;z++) phi(i,j,mn_start+z)=vector[z];}
    delete[] vector;
}
//#####################################################################
// Function Create_Erosion_Kernel
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Create_Erosion_Kernel(ARRAYS<VECTOR<int,3> >& kernel,const int xradius,const int yradius,const int zradius)
{
    if(xradius<0 || yradius<0 || zradius<0){std::cerr<<"invalid kernel: "<<xradius<<" "<<yradius<<" "<<zradius<<std::endl;exit(1);}
    kernel.Resize(-xradius,xradius,-yradius,yradius,-zradius,zradius);
    for(int i=-xradius;i<=xradius;i++)for(int j=-yradius;j<=yradius;j++)for(int ij=-zradius;ij<=zradius;ij++){
        T i1=xradius>0?T(i)/xradius:0,j1=yradius>0?T(j)/yradius:0,ij1=zradius>0?T(ij)/zradius:0;
        if(xradius>1 || yradius>1 || zradius>1) kernel(i,j,ij)=sqr(i1)+sqr(j1)+sqr(ij1)<1;
        else kernel(i,j,ij)=sqr(i1)+sqr(j1)+sqr(ij1)<=1;}
}
//#####################################################################
// Function Dilation
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Dilation(const int xradius,const int yradius,const int zradius)
{
    if(xradius<0 || yradius<0 || zradius<0){std::cerr<<"invalid kernel: "<<xradius<<" "<<yradius<<" "<<zradius<<std::endl;exit(1);}
    int m=phi.m,m_start=phi.domain.min_corner.x,m_end=phi.domain.max_corner.x;
    int n=phi.n,n_start=phi.domain.min_corner.y,n_end=phi.domain.max_corner.y;
    int mn=phi.mn,mn_start=phi.domain.min_corner.z,mn_end=phi.domain.max_corner.z;
    
    ARRAYS<VECTOR<int,3> > mask(m_start,m_end,n_start,n_end,mn_start,mn_end);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++)mask(i,j,ij)=phi(i,j,ij)>0; // indicator of outside
    Window_Sum_3D(mask,1,1,1); // spread "outside" to the immediate neighbors
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++)mask(i,j,ij)=phi(i,j,ij)<=0 && mask(i,j,ij)>0; // inside + neighbor of outside
    
    ARRAYS<VECTOR<T,3> > new_phi=phi;
    ARRAYS<VECTOR<int,3> > kernel;Create_Erosion_Kernel(kernel,xradius,yradius,zradius);
    int num_gridpoints=m*n*mn;long int count=0;
    printf("Dilation: [%d %d]x[%d %d]x[%d %d]=>%d gridpoints\n",m_start,m_end,n_start,n_end,mn_start,mn_end,num_gridpoints);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++){
        if(++count%100000==0) std::cout<<int(100.0*count/num_gridpoints)<<"%"<<std::endl;
        if(mask(i,j,ij)) //turn all neighbors to "inside"
            for(int x=-xradius;x<=xradius;++x)for(int y=-yradius;y<=yradius;++y)for(int z=-zradius;z<=zradius;++z)if(kernel(x,y,z)){
                T& new_phi_value=new_phi(clamp(i+x,m_start,m_end),clamp(j+y,n_start,n_end),clamp(ij+z,mn_start,mn_end));
                new_phi_value=-fabs(new_phi_value);}}
    phi=new_phi;std::cout<<"done"<<std::endl;
}

//#####################################################################
// Function Erosion
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Erosion(const int xradius,const int yradius,const int zradius)
{
    if(xradius<0 || yradius<0 || zradius<0){printf("invalid kernel: %d %d %d!\n",xradius,yradius,zradius);exit(1);}
    int m=phi.m,m_start=phi.domain.min_corner.x,m_end=phi.domain.max_corner.x;
    int n=phi.n,n_start=phi.domain.min_corner.y,n_end=phi.domain.max_corner.y;
    int mn=phi.mn,mn_start=phi.domain.min_corner.z,mn_end=phi.domain.max_corner.z;
    
    ARRAYS<VECTOR<int,3> > mask(m_start,m_end,n_start,n_end,mn_start,mn_end);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++)mask(i,j,ij)=phi(i,j,ij)<0; // indicator of inside
    Window_Sum_3D(mask,1,1,1); // spread to the immediate neighbors
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++)mask(i,j,ij)=phi(i,j,ij)>0 && mask(i,j,ij)>0; // outside + neighbor of inside
    
    ARRAYS<VECTOR<T,3> > new_phi=phi;
    ARRAYS<VECTOR<int,3> > kernel;Create_Erosion_Kernel(kernel,xradius,yradius,zradius);
    int num_gridpoints=m*n*mn;long int count=0;  
    printf("Erosion: [%d %d]x[%d %d]x[%d %d]=>%d gridpoints\n",m_start,m_end,n_start,n_end,mn_start,mn_end,num_gridpoints);
    for(int i=m_start;i<=m_end;i++)for(int j=n_start;j<=n_end;j++)for(int ij=mn_start;ij<=mn_end;ij++){
        if(++count%100000==0)printf("%d%% ",int(100.0*count/num_gridpoints));
        if(mask(i,j,ij)) // turn all neighbors to "outside"
            for(int x=-xradius;x<=xradius;++x)for(int y=-yradius;y<=yradius;++y)for(int z=-zradius;z<=zradius;++z)if(kernel(x,y,z)){
                T& new_phi_value=new_phi(clamp(i+x,m_start,m_end),clamp(j+y,n_start,n_end),clamp(ij+z,mn_start,mn_end));
                new_phi_value=fabs(new_phi_value);}}
    phi=new_phi;std::cout<<"done"<<std::endl;
}

//#####################################################################
// Function Closing
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Closing(const int xradius,const int yradius,const int zradius)
{
    std::cout<<"Closing: "<<std::endl;
    Dilation(xradius,yradius,zradius);
    Erosion(xradius,yradius,zradius);
    std::cout<<"Closing done: "<<std::endl;
}
//#####################################################################
// Function Opening
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Opening(const int xradius,const int yradius,const int zradius)
{
    std::cout<<"Opening: "<<std::endl;
    Erosion(xradius,yradius,zradius);    
    Dilation(xradius,yradius,zradius);
    std::cout<<"Opening done: "<<std::endl;
}
//#####################################################################
// Function Smooth
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Smooth(const int xradius,const int yradius,const int zradius)
{
    Window_Sum_3D(levelset.phi,xradius,yradius,zradius);
}
//#####################################################################
// Function Subtract_Set_Same_Size
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Subtract_Set_Same_Size(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset_input)
{
    std::cout<<"Subtract_Set_Same_Size()"<<std::endl;
    if(levelset.grid != levelset_input.grid){printf("unequal grids!\n");exit(1);}
    if(!ARRAYS<VECTOR<T,3> >::Equal_Dimensions(phi,levelset_input.phi)){printf("unequal phi sizes!\n");exit(1);}
    for(int i=phi.domain.min_corner.x;i<=phi.domain.max_corner.x;i++)for(int j=phi.domain.min_corner.y;j<=phi.domain.max_corner.y;j++)for(int ij=phi.domain.min_corner.z;ij<=phi.domain.max_corner.z;ij++)
        phi(i,j,ij)=max(phi(i,j,ij),-levelset_input.phi(i,j,ij));
    std::cout<<"done"<<std::endl;
}
//#####################################################################
// Function Subtract_Set
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Subtract_Set(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset_input)
{   
    std::cout<<"Subtract_Set()"<<std::endl;
    std::cout<<"grid: "<<grid<<std::endl;
    std::cout<<"input grid: "<<levelset_input.grid<<std::endl;
    for(int i=phi.domain.min_corner.x;i<=phi.domain.max_corner.x;i++)for(int j=phi.domain.min_corner.y;j<=phi.domain.max_corner.y;j++)for(int ij=phi.domain.min_corner.z;ij<=phi.domain.max_corner.z;ij++)
        phi(i,j,ij)=max(phi(i,j,ij),-levelset_input.Extended_Phi(grid.X(i,j,ij)));
    std::cout<<"done"<<std::endl;
}
//#####################################################################
// Function Add_Set_Same_Size
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Add_Set_Same_Size(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset_input)
{
    std::cout<<"Add_Set_Same_Size()"<<std::endl;
    if(levelset.grid != levelset_input.grid){printf("unequal grids!\n");exit(1);}
    if(!ARRAYS<VECTOR<T,3> >::Equal_Dimensions(phi,levelset_input.phi)){printf("unequal phi sizes!\n");exit(1);}
    for(int i=phi.domain.min_corner.x;i<=phi.domain.max_corner.x;i++)for(int j=phi.domain.min_corner.y;j<=phi.domain.max_corner.y;j++)for(int ij=phi.domain.min_corner.z;ij<=phi.domain.max_corner.z;ij++)
        phi(i,j,ij)=max(phi(i,j,ij),levelset_input.phi(i,j,ij));
    std::cout<<"done"<<std::endl;
}
//#####################################################################
// Function Add_Set
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Add_Set(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset_input)
{
    std::cout<<"Add_Set()"<<std::endl;
    std::cout<<"grid: "<<levelset.grid<<std::endl;
    std::cout<<"input grid: "<<levelset_input.grid<<std::endl;
    for(int i=phi.domain.min_corner.x;i<=phi.domain.max_corner.x;i++)for(int j=phi.domain.min_corner.y;j<=phi.domain.max_corner.y;j++)for(int ij=phi.domain.min_corner.z;ij<=phi.domain.max_corner.z;ij++)
        phi(i,j,ij)=max(phi(i,j,ij),levelset_input.Extended_Phi(grid.X(i,j,ij)));
    std::cout<<"done"<<std::endl;
}
//#####################################################################
// Function Output_Subset
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Output_Subset(LEVELSET_3D<GRID<VECTOR<T,3> > >& subset,const int m0,const int m1,const int n0,const int n1,const int mn0,const int mn1) const
{
    VECTOR<int,3> M0(m0,n0,mn0),M1(m1,n1,mn1);
    if(&levelset.phi==&subset.phi){std::cerr<<"output needs to be different from input!"<<std::endl;exit(1);}
    bool arguments_constitute_a_subgrid=(1<=m0 && m0<=m1 && m1<=grid.counts.x && 1<=n0 && n0<=n1 && n1<=grid.counts.y && 1<=mn0 && mn0<=mn1 && mn1<=grid.counts.z);
    if(!arguments_constitute_a_subgrid){std::cerr<<"invalid subgrid!"<<std::endl;exit(1);}
    subset.grid=GRID<VECTOR<T,3> >(M1-M0+1,RANGE<VECTOR<T,3> >(grid.X(M0),grid.X(M1)));
    subset.phi.Resize(m0,m1,n0,n1,mn0,mn1); // should be fine to have phi dimensions not start from 1's
    ARRAYS<VECTOR<T,3> >::Get(subset.phi,levelset.phi);
}
//#####################################################################
// Function Num_Points
//#####################################################################
template<class T> int LEVELSET_PROCESSOR<T>::
Num_Points(const int start,const int end,const int stride,const bool include_end) const
{
    int m=end-start+1;if(stride==1) return m;
    int num_full_strides=m/stride,num_left_over=m%stride;
    int num_points_in_range=num_full_strides+((num_left_over>0)?1:0);
    if(!include_end || num_left_over==1) return num_points_in_range;
    else return num_points_in_range+1;
}
//#####################################################################
// Function Output_Subsample
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Output_Subsample(LEVELSET_3D<GRID<VECTOR<T,3> > >& subsample,const int x_stride,const int y_stride,const int z_stride,bool include_ends) const
{   
    if(&levelset.phi==&subsample.phi){std::cerr<<"output needs to be different from input!"<<std::endl;exit(1);}
    if(x_stride<=0 || y_stride<=0 || z_stride<=0){printf("invalid strides: %d %d %d\n",x_stride,y_stride,z_stride);exit(1);}
    int m=Num_Points(1,grid.counts.x,x_stride,include_ends);
    int n=Num_Points(1,grid.counts.y,y_stride,include_ends);
    int mn=Num_Points(1,grid.counts.z,z_stride,include_ends);
    printf("m %d, n %d, mn %d\n", m,n,mn);
    T new_dx=grid.dX.x*x_stride,new_dy=grid.dX.y*y_stride,new_dz=grid.dX.z*z_stride;
    T new_xmax=grid.domain.min_corner.x+new_dx*(m-1),new_ymax=grid.domain.min_corner.y+new_dy*(n-1),new_zmax=grid.domain.min_corner.z+new_dz*(mn-1);
    printf("new_xmax = grid.xmax + %g cells\n",(new_xmax-grid.domain.max_corner.x)/grid.dX.x);
    printf("new_ymax = grid.ymax + %g cells\n",(new_ymax-grid.domain.max_corner.y)/grid.dX.y);
    printf("new_zmax = grid.zmax + %g cells\n",(new_zmax-grid.domain.max_corner.z)/grid.dX.z);
    subsample.grid=GRID<VECTOR<T,3> >(m,n,mn,grid.domain.min_corner.x,new_xmax,grid.domain.min_corner.y,new_ymax,grid.domain.min_corner.z,new_zmax);
    subsample.phi.Resize(1,m,1,n,1,mn);
    for(int i=1;i<=m;++i)for(int j=1;j<=n;++j)for(int k=1;k<=mn;++k){
        int i1=1+(i-1)*x_stride,j1=1+(j-1)*y_stride,k1=1+(k-1)*z_stride;
        if(i1<=phi.domain.max_corner.x && j1<=phi.domain.max_corner.y && k1<=phi.domain.max_corner.z) subsample.phi(i,j,k)=phi(i1,j1,k1);//a necessary check because output.grid may well be larger
        else subsample.phi(i,j,k)=levelset.Extended_Phi(subsample.grid.X(i,j,k));}
}
//#####################################################################
// Function Output_Resample
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Output_Resample(LEVELSET_3D<GRID<VECTOR<T,3> > >& supersample,const int m,const int n,const int mn) const
{
    supersample.grid=GRID<VECTOR<T,3> >(m,n,mn,grid.Domain());
    Output_Resample(supersample);
}
//#####################################################################
// Function Smooth_By_Curvature
//#####################################################################
template<class T> void LEVELSET_PROCESSOR<T>::
Output_Resample(LEVELSET_3D<GRID<VECTOR<T,3> > >& resample) const
{
    for(int i=0;i<resample.grid.counts.x;i++)for(int j=0;j<resample.grid.counts.y;j++)for(int ij=0;ij<resample.grid.counts.z;ij++)
        resample.phi(i,j,ij)=levelset.Phi(resample.grid.X(i,j,ij));
}
//#####################################################################
template class LEVELSET_PROCESSOR<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_PROCESSOR<double>;
#endif
