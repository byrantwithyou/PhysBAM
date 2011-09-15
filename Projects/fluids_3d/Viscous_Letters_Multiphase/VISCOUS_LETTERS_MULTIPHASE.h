//#####################################################################
// Copyright 2005, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VISCOUS_LETTERS_MULTIPHASE  
//##################################################################### 
#ifndef __VISCOUS_LETTERS_MULTIPHASE__
#define __VISCOUS_LETTERS_MULTIPHASE__


#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

using namespace PhysBAM;

template<class T,class RW=T>
class VISCOUS_LETTERS_MULTIPHASE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;
    using BASE::output_directory;using BASE::restart;using BASE::restart_frame;using BASE::verbose_dt;using BASE::data_directory;

    RIGID_BODY_LIST<T,TV> letters;

    VISCOUS_LETTERS_MULTIPHASE(const int resolution,const bool restart,const int restart_frame)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(10,fluids_parameters.WATER)
    {
        fluids_parameters.grid.Initialize(35*resolution+1,20*resolution+1,35*resolution+1,-.3,.3,-.22,-.22+(.6*20/35),-.3,.3);
        LOG::cout<<"\ngrid: dx,dy,dz "<<fluids_parameters.grid.dx<<", "<<fluids_parameters.grid.dy<<", "<<fluids_parameters.grid.dz<<std::endl;
        fluids_parameters.domain_walls[3][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.domain_walls[2][2]=false;
        frame_rate=480;
        last_frame=int(T(20)*frame_rate);
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_reacting_flow=false;
        fluids_parameters.use_old_velocities_for_boundary_conditions=true;
        fluids_parameters.incompressible_iterations=100;
        fluids_parameters.implicit_viscosity_iterations=50;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.second_order_cut_cell_method=true;
        write_output_files=true;fluids_parameters.write_debug_data=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        output_directory="Viscous_Letters_Multiphase/output";
        fluids_parameters.reseeding_frame_rate=10;

//        if(restart&&restart_frame>=181){
            LOG::cout<<"RESTARTING at "<<restart_frame<<std::endl;
            // densities
            fluids_parameters.densities(1)=(T)10;
            fluids_parameters.densities(2)=(T)50;
            fluids_parameters.densities(3)=(T)20;
            fluids_parameters.densities(4)=(T)50;
            fluids_parameters.densities(5)=(T)30;
            fluids_parameters.densities(6)=(T)100;
            fluids_parameters.densities(7)=(T)10;
            fluids_parameters.densities(8)=(T)50;
            fluids_parameters.densities(9)=(T)1000;
            fluids_parameters.densities(10)=(T)1;
            // dirichlet regions
            fluids_parameters.dirichlet_regions(10)=true;
            // surface_tensions
            fluids_parameters.surface_tensions(1,9)=fluids_parameters.surface_tensions(9,1)=.00001;
            fluids_parameters.surface_tensions(2,9)=fluids_parameters.surface_tensions(9,2)=.00002;
            fluids_parameters.surface_tensions(3,9)=fluids_parameters.surface_tensions(9,3)=.00001;
            fluids_parameters.surface_tensions(4,9)=fluids_parameters.surface_tensions(9,4)=.00005;
            fluids_parameters.surface_tensions(5,9)=fluids_parameters.surface_tensions(9,5)=.00001;
            fluids_parameters.surface_tensions(6,9)=fluids_parameters.surface_tensions(9,6)=.00006;
            fluids_parameters.surface_tensions(7,9)=fluids_parameters.surface_tensions(9,7)=.00002;
            fluids_parameters.surface_tensions(8,9)=fluids_parameters.surface_tensions(9,8)=.00003;
//            fluids_parameters.surface_tensions(10,9)=fluids_parameters.surface_tensions(9,10)=.00001;//}
/*        else{
            // densities
            fluids_parameters.densities(1)=(T)3000;
            fluids_parameters.densities(2)=(T)2000;
            fluids_parameters.densities(3)=(T)2500;
            fluids_parameters.densities(4)=(T)1500;
            fluids_parameters.densities(5)=(T)2000;
            fluids_parameters.densities(6)=(T)1200;
            fluids_parameters.densities(7)=(T)2700;
            fluids_parameters.densities(8)=(T)1200;
            fluids_parameters.densities(9)=(T)1000;
            fluids_parameters.densities(10)=(T)1;
            // dirichlet regions
            fluids_parameters.dirichlet_regions(10)=true;
            // viscosities
            T letters_viscosity=5;
            fluids_parameters.viscosities(1)=letters_viscosity*fluids_parameters.densities(1);
            fluids_parameters.viscosities(2)=letters_viscosity*fluids_parameters.densities(2);
            fluids_parameters.viscosities(3)=letters_viscosity*fluids_parameters.densities(3);
            fluids_parameters.viscosities(4)=letters_viscosity*fluids_parameters.densities(4);
            fluids_parameters.viscosities(5)=letters_viscosity*fluids_parameters.densities(5);
            fluids_parameters.viscosities(6)=letters_viscosity*fluids_parameters.densities(6);
            fluids_parameters.viscosities(7)=letters_viscosity*fluids_parameters.densities(7);
            fluids_parameters.viscosities(8)=letters_viscosity*fluids_parameters.densities(8);
            fluids_parameters.viscosities(9)=0;
            fluids_parameters.viscosities(10)=0;}*/

        int air_region=10;
        for(int i=1;i<=fluids_parameters.number_of_regions;i++){
            BOUNDARY_PHI_WATER<T,GRID<TV> >* boundary=new BOUNDARY_PHI_WATER<T,GRID<TV> >();
            boundary->Set_Velocity_Pointer(fluids_parameters.incompressible->projection.face_velocities);
            if(i==air_region)boundary->sign=-1;
            fluids_parameters.phi_boundary_multiphase(i)=boundary;}

        LOG::cout<<"DENSITIES: "<<fluids_parameters.densities<<std::endl;
        LOG::cout<<"VISCOSITIES: "<<fluids_parameters.viscosities<<std::endl;
        LOG::cout<<"SURFACE TENSION: "<<std::endl;
        for(int i=1;i<=fluids_parameters.number_of_regions;i++){
            for(int j=1;j<=fluids_parameters.number_of_regions;j++) LOG::cout<<fluids_parameters.surface_tensions(i,j)<<" ";LOG::cout<<std::endl;}
    }

    virtual ~VISCOUS_LETTERS_MULTIPHASE()
    {}

//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    T letters_scale=.15;
    T letters_height=letters_scale*.2;
    T surface_height=-.08;
    T letters_position_y=-.04;
    BOX_1D<T> letters_bounding_box_y(letters_position_y-letters_height/2,letters_position_y+letters_height/2);

    // read in implicit surfaces for the letters
    const std::string siggraph("SIGGRAPH");T two_pi_over_360=(T)2*pi/(T)360;
    ARRAY<PAIR<VECTOR<T,3>,VECTOR<T,3> > > maya_coordinates;Get_Maya_Coordinates(maya_coordinates);
    std::cout<<"\nReading "<<std::flush;
    for(unsigned int i=0;i<siggraph.length();i++){
        std::cout<<siggraph[i]<<" "<<std::flush;
        int id=letters.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/Letters/"+siggraph[i],letters_scale,true,true,false,false);
        letters(id)->frame.t=letters_scale*maya_coordinates(i+1).x*VECTOR<T,3>(1,1,-1)+VECTOR<T,3>(0,letters_position_y,0);
        letters(id)->frame.r=QUATERNION<T>(MATRIX<T,4>::Rotation_Matrix_Z_Axis(two_pi_over_360*(maya_coordinates(i+1).y.z))*
                                           MATRIX<T,4>::Rotation_Matrix_Y_Axis(-two_pi_over_360*(maya_coordinates(i+1).y.y))*
                                           MATRIX<T,4>::Rotation_Matrix_X_Axis(-two_pi_over_360*(maya_coordinates(i+1).y.x)));}

    std::cout<<"\nInitializing phis..."<<std::flush;
    ARRAY<ARRAY<T,VECTOR<int,3> > >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(int i=1;i<=8;i++)ARRAY<T,VECTOR<int,3> >::copy(5*fluids_parameters.grid.dx,phis(i));
    for(CELL_ITERATOR iterator(fluids_parameters.grid);iterator.Valid();iterator.Next()){
        VECTOR<int,3> cell_index=iterator.Cell_Index();VECTOR<T,3> X=iterator.Location();
        phis(9)(cell_index)=X.y-surface_height;phis(10)(cell_index)=-phis(9)(cell_index);
        if(letters_bounding_box_y.Lazy_Inside(VECTOR<T,1>(X.y)))for(int i=1;i<=letters.rigid_bodies.m;i++){
            T phi=letters(i)->Implicit_Geometry_Extended_Value(X)-.25*letters_height;
            T bandwidth=3*fluids_parameters.grid.dx;
            if(phi<bandwidth){phis(i)(cell_index)=phi;for(int j=1;j<=10;j++)if(i!=j)phis(j)(cell_index)=max(phis(j)(cell_index),-phi);}}}
    std::cout<<"DONE"<<std::endl;
}
//#####################################################################
// Function Get_Maya_Coordinates
//#####################################################################
void Get_Maya_Coordinates(ARRAY<PAIR<VECTOR<T,3>,VECTOR<T,3> > >& maya_coordinates)
{
    // .x :   x,y,z translations wrt origin
    // .y :   rotations in degrees about x,y,z axes
    maya_coordinates.Resize(8);
    maya_coordinates(1).x=VECTOR<T,3>(.003,0,-1.34);
    maya_coordinates(1).y=VECTOR<T,3>(-90.454,-35.512,179.9);
    maya_coordinates(2).x=VECTOR<T,3>(0.598,0,-1.195);
    maya_coordinates(2).y=VECTOR<T,3>(91.278,153.014,1.349);
    maya_coordinates(3).x=VECTOR<T,3>(1.125,0,-0.669);
    maya_coordinates(3).y=VECTOR<T,3>(90,145.362,0);
    maya_coordinates(4).x=VECTOR<T,3>(1.291,0,0.273);
    maya_coordinates(4).y=VECTOR<T,3>(90,96.478,0);
    maya_coordinates(5).x=VECTOR<T,3>(0.824,0,1.075);
    maya_coordinates(5).y=VECTOR<T,3>(90,50.711,0);
    maya_coordinates(6).x=VECTOR<T,3>(-0.157,0,1.269);
    maya_coordinates(6).y=VECTOR<T,3>(90,-4.945,0);
    maya_coordinates(7).x=VECTOR<T,3>(-1.043,0,0.948);
    maya_coordinates(7).y=VECTOR<T,3>(89.164,-65.138,1.522);
    maya_coordinates(8).x=VECTOR<T,3>(-1.288,0,0.132);
    maya_coordinates(8).y=VECTOR<T,3>(88.598,97.584,-1.415);
}
//##################################################################### 
};      
#endif


