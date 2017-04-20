//#####################################################################
// Class GROUND_ELEMENT
//##################################################################### 
//
//#####################################################################
// Fedkiw - February 14, 2002
// Nguyen - March 25, 2003
//#####################################################################
#ifndef __RISING_SMOKE__
#define __RISING_SMOKE__

#include <Grid_Tools/Grids/GRID.h>
#include "../SMOKE_3D_EXAMPLE.h"
namespace PhysBAM{

template <class T>
class RISING_SMOKE:public SMOKE_3D_EXAMPLE<T>
{
public:
    T rho,rho_bottom,rho_top;
    
    RISING_SMOKE()
    {
        start_frame=0;
        end_frame=200;
        frame_rate=24;
        m=20;n=50;mn=50;
        xmin=(T)0;xmax=(T).8;ymin=(T)0;ymax=(T)2;zmin=(T)0;zmax=(T)2;
        use_vorticity_confinement=true;confinement_parameter=(T).3;
        komolgorov=(T)0;
        incompressible_enforce_compatibility=false;
        rho=(T)1.;T_air=T(0.);T_burnt=(T)2.;
        gravity=TV();
        cooling_constant=(T)0;
        write_matlab_files=true;
        write_output_files=true;
        matlab_directory="Rising_Smoke/matlab";
        output_directory="Rising_Smoke/output";
    }

    ~RISING_SMOKE()
    {}

//#####################################################################
// Function Update_Source
//#####################################################################
void Update_Sources(const GRID<TV>& grid,ARRAY<T,VECTOR<int,3> >& u,ARRAY<T,VECTOR<int,3> >& v,ARRAY<T,VECTOR<int,3> >& w,ARRAY<T,VECTOR<int,3> >& density, 
                                  ARRAY<T,VECTOR<int,3> >& temperature, const T time)
{
    int i,j,ij;
    int m=grid.m,n=grid.n,mn=grid.mn;
    T source_xmax = (T).1;
    float source_zmax = (T).3;

    for(i=0;i<m;i++) for(j=0;j<n;j++) for(ij=0;ij<mn;ij++){
        T x=grid.x(i);
        T y=grid.y(j);
        T z=grid.z(ij);
        if(z <= source_zmax && sqr(x-.4)+sqr(y-.5)<=sqr(source_xmax)){
            w(i,j,ij)=(T)2.0;
            density(i,j,ij)=rho;temperature(i,j,ij)=T_burnt;}}

    // keep density >= 0 and T >=0
    for(i=0;i<m;i++) for(j=0;j<n;j++) for(ij=0;ij<mn;ij++){
        if(density(i,j,ij) < 0) density(i,j,ij)=0;
        if(temperature(i,j,ij) < T_air) temperature(i,j,ij)=T_air;}
}
//#####################################################################
// Function Set_Domain_Boundary_Conditions
//#####################################################################
void Set_Domain_Boundary_Conditions(ARRAY<bool,VECTOR<int,3> >& psi_N,ARRAY<bool,VECTOR<int,3> >& psi_D,ARRAY<T,VECTOR<int,3> >& p)
{
    int i,j,ij;
    // set up Dirichlet boundary conditions on the outside edges
    for(i=0;i<m;i++) for(j=0;j<n;j++){psi_D(i,j,0)=psi_D(i,j,mn+1)=1;p(i,j,0)=p(i,j,mn+1)=0;}
    for(i=0;i<m;i++) for(ij=0;ij<mn;ij++){psi_D(i,0,ij)=psi_D(i,n+1,ij)=1;p(i,0,ij)=p(i,n+1,ij)=0;}
    for(j=0;j<n;j++) for(ij=0;ij<mn;ij++){psi_D(0,j,ij)=psi_D(m+1,j,ij)=1;p(0,j,ij)=p(m+1,j,ij)=0;}
}
//#####################################################################
// Function Update_Forces
//#####################################################################
void Update_Forces(const GRID<TV>& grid,const ARRAY<T,VECTOR<int,3> >& u,const ARRAY<T,VECTOR<int,3> >& v,const ARRAY<T,VECTOR<int,3> >& w,
                                const ARRAY<T,VECTOR<int,3> >& density_ghost,const ARRAY<T,VECTOR<int,3> >& T_ghost,
                                ARRAY<T,VECTOR<int,3> >& force_x,ARRAY<T,VECTOR<int,3> >& force_y,ARRAY<T,VECTOR<int,3> >& force_z,const T time)
{
    int i,j,ij;
    // control buoyancy force
    for(i=0;i<m;i++) for(j=0;j<n+1;j++) for(ij=0;ij<mn;ij++){ // y-direction forces only
        T rho_atm=rho_bottom+(rho_top-rho_bottom)*(grid.y(j)-ymin)/(ymax-ymin);
        T den=(density_ghost(i,j-1,ij)+density_ghost(i,j,ij))/2,temp=(T_ghost(i,j-1,ij)+T_ghost(i,j,ij))/2;
        force_y(i,j,ij)=-(rho_atm-den);}
}    
//#####################################################################
// Function Write_Data_Files
//#####################################################################
void Write_Data_File(const GRID<TV>& grid,const ARRAY<T,VECTOR<int,3> >& u,const ARRAY<T,VECTOR<int,3> >& v,const ARRAY<T,VECTOR<int,3> >& w,
                                  const ARRAY<T,VECTOR<int,3> >& density,const ARRAY<T,VECTOR<int,3> >& temperature,const int frame)
{    
    int j,ij;char filename[256];
    GRID<TV> grid_2d(n,mn,ymin,ymax,zmin,zmax);grid_2d.Set_MAC_Grid();

    if(write_matlab_files){
        if(!Directory_Exists(matlab_directory.c_str())) Create_Directory(matlab_directory.c_str());
        MATLAB_OUTPUT matlab_output;ARRAY<T,VECTOR<int,2> > output(1,n,1,mn);
        sprintf(filename,"%s/header",matlab_directory.c_str());matlab_output.Write_Header_File(filename,grid_2d,frame);
        for(j=0;j<n;j++) for(ij=0;ij<mn;ij++)output(j,ij)=density(m/2,j,ij);
        sprintf(filename,"%s/density",matlab_directory.c_str());matlab_output.Write_Output_File(filename,output,frame);
        for(j=0;j<n;j++) for(ij=0;ij<mn;ij++)output(j,ij)=temperature(m/2,j,ij);
        sprintf(filename,"%s/temperature",matlab_directory.c_str());matlab_output.Write_Output_File(filename,output,frame);
        for(j=0;j<n;j++) for(ij=0;ij<mn;ij++)output(j,ij)=(T).5*(v(m/2,j,ij)+v(m/2,j+1,ij));
        sprintf(filename,"%s/velocity0",matlab_directory.c_str());matlab_output.Write_Output_File(filename,output,frame);
        for(j=0;j<n;j++) for(ij=0;ij<mn;ij++)output(j,ij)=(T).5*(w(m/2,j,ij)+w(m/2,j,ij+1));
        sprintf(filename,"%s/velocity1",matlab_directory.c_str());matlab_output.Write_Output_File(filename,output,frame);}

    if(write_output_files){
        std::ofstream output;
        if(!Directory_Exists(output_directory.c_str())) Create_Directory(output_directory.c_str());
        if(frame == start_frame){
            sprintf(filename,"%s/common/grid",output_directory.c_str());output.open(filename,std::ios::binary);grid.Write(output);output.close();}
        sprintf(filename,"%s/%d/density",output_directory.c_str(),frame);output.open(filename,std::ios::binary);density.Write(output);output.close();
        sprintf(filename,"%s/%d/temperature",output_directory.c_str(),frame);output.open(filename,std::ios::binary);temperature.Write(output);output.close();
        sprintf(filename,"%s/%d/velocities",output_directory.c_str(),frame);output.open(filename,std::ios::binary);u.Write(output);v.Write(output);output.close();
        sprintf(filename,"%s/common/last_frame",output_directory.c_str());output.open(filename,std::ios::out|std::ios::trunc);output<<frame<<std::endl;output.close();}

}
//#####################################################################
};      
}
#endif


