//#####################################################################
// Copyright 2003-2005, Jiayi Chong, Cynthia Lau, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include "COMMANDS.h"
#include "LEVELSET_PROCESSOR.h"
#include "MISC.h"
using namespace PhysBAM;
//#####################################################################

template<class T,class RW>
struct COMMAND_LEVELSET_OP:public COMMAND_BASE
{
    virtual int Min_Inputs() { return 1; }
    virtual int Max_Inputs() { return 1; }
    virtual int Num_Outputs() { return 1; }
    virtual int Do()
    {
        GRID<VECTOR<T,3> > grid, grid2;
        ARRAYS<VECTOR<T,3> > phi, phi2;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset(grid,phi), levelset2(grid2,phi2);
        FILE_UTILITIES::Read_From_File<RW>(Get_Input(1),levelset);
        DoOp(levelset,levelset2);
        if(GrowBoundingBox()) FILE_UTILITIES::Write_To_File<RW>(Get_Output(),levelset2);
        else FILE_UTILITIES::Write_To_File<RW>(Get_Output(),levelset);
        return 1;
    }

    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) = 0;
    virtual bool GrowBoundingBox() { return false; }
};

template<class T,class RW>
struct COMMAND_LEVELSET_BINARY_OP:public COMMAND_BASE
{
    virtual int Min_Inputs() { return 2; }
    virtual int Max_Inputs() { return 2; }
    virtual int Num_Outputs() { return 1; }
    virtual int Do()
    {
        GRID<VECTOR<T,3> > grid1,grid2;
        ARRAYS<VECTOR<T,3> > phi1,phi2;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset1(grid1,phi1);
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset2(grid2,phi2);
        FILE_UTILITIES::Read_From_File<RW>(Get_Input(1),levelset1);
        FILE_UTILITIES::Read_From_File<RW>(Get_Input(2),levelset2);
        DoOp(levelset1,levelset2);
        FILE_UTILITIES::Write_To_File<RW>(Get_Output(),levelset2);
        return 1;
    }

    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset1,const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) = 0;
};

template<class T,class RW>
struct COMMAND_FAST_MARCHING:public COMMAND_LEVELSET_OP<T,RW>
{
    virtual std::string Name() { return "fmm"; }
    virtual std::string Description() { return "does the fast marching method."; }
    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2)
    {
        std::cout<<"fmm: dx is "<<levelset.grid.dX.x<<std::endl;
        std::cout<<"fmm: dy is "<<levelset.grid.dX.y<<std::endl;
        std::cout<<"fmm: dz is "<<levelset.grid.dX.z<<std::endl;
        levelset.Fast_Marching_Method();
    }
};

template<class T,class RW>
struct COMMAND_ADD_CONSTANT:public COMMAND_LEVELSET_OP<T,RW>
{
    T dphi;

    COMMAND_ADD_CONSTANT() {
        Add_Float_Arg("dphi",dphi);
    }
    virtual std::string Name() { return "add_constant"; }
    virtual std::string Description() { return "adds a constant phi."; }

    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2)
    {
        levelset.phi+=dphi;
    }
};

template<class T,class RW>
struct COMMAND_MOTION_MEAN_CURVATURE:public COMMAND_LEVELSET_OP<T,RW>
{
    typedef COMMAND_LEVELSET_OP<T,RW> BASE;using BASE::Add_Int_Arg;

    int iters;
    T sigma;

    COMMAND_MOTION_MEAN_CURVATURE():iters(5),sigma(1) {
        Add_Int_Arg("iters",iters);
        Add_Float_Arg("sigma",sigma);
    }
    virtual std::string Name() { return "mmc"; }
    virtual std::string Description() { return "smooth with mmc, sigma = strength."; }

    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2)
    {
        LEVELSET_PROCESSOR<T> proc(levelset);
        proc.Smooth_By_Curvature(iters,sigma);
    }
};

template<class T,class RW>
struct COMMAND_KERNEL_OP:public COMMAND_LEVELSET_OP<T,RW>
{
    typedef COMMAND_LEVELSET_OP<T,RW> BASE;using BASE::Add_Int_Arg;

    int xradius,yradius,zradius;
    COMMAND_KERNEL_OP()
        :xradius(1),yradius(1),zradius(1)
    {
        Add_Int_Arg("xradius",xradius);
        Add_Int_Arg("yradius",yradius);
        Add_Int_Arg("zradius",zradius);
    }
};

template<class T,class RW>
struct COMMAND_EROSION:public COMMAND_KERNEL_OP<T,RW>
{
    typedef COMMAND_KERNEL_OP<T,RW> BASE;using BASE::xradius;using BASE::yradius;using BASE::zradius;

    virtual std::string Name() { return "erosion"; }
    virtual std::string Description() { return "erode the levelset."; }
    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        LEVELSET_PROCESSOR<T> proc(levelset);
        proc.Erosion(xradius,yradius,zradius);
    }
};

template<class T,class RW>
struct COMMAND_DILATION:public COMMAND_KERNEL_OP<T,RW>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef COMMAND_KERNEL_OP<T,RW> BASE;using BASE::xradius;using BASE::yradius;using BASE::zradius;
    //
    virtual std::string Name() { return "dilation"; }
    virtual std::string Description() { return "dilate the levelset."; }
    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        TV_INT radius(xradius,yradius,zradius);
        levelset2.grid.Initialize(levelset.grid.counts+2*radius,
            RANGE<TV>(levelset.grid.domain.min_corner-TV(radius)*levelset.grid.dX,levelset.grid.domain.max_corner+TV(radius)*levelset.grid.dX));
        levelset2.phi.Resize(levelset2.grid.Domain_Indices());
        bool use_x_default,use_y_default,use_z_default;
        T phi_value = levelset2.grid.dX.x/2;
        std::cout << "levelset2 phi ranges: " << levelset2.phi.m_start << " " << levelset2.phi.m_end << ", " 
            << levelset2.phi.n_start << " " << levelset2.phi.n_end << ", " 
            << levelset2.phi.mn_start << " " << levelset2.phi.mn_end << std::endl;
        for (int i=levelset2.phi.m_start;i<=levelset2.phi.m_end;i++){
            use_x_default=(i<levelset2.phi.m_start+xradius || i>levelset2.phi.m_end-xradius);
            for (int j=levelset2.phi.n_start;j<=levelset2.phi.n_end;j++){
                use_y_default=(use_x_default || j<levelset2.phi.n_start+yradius || j>levelset2.phi.n_end-yradius);
                for (int k=levelset2.phi.mn_start;k<=levelset2.phi.mn_end;k++){
                    use_z_default=(use_y_default || k<levelset2.phi.mn_start+zradius || k>levelset2.phi.mn_end-zradius);
                    levelset2.phi(i,j,k) = (use_z_default)?phi_value:levelset.phi(i-xradius,j-yradius,k-zradius);
                    //if (!use_default) std::cout << "using original levelset value " << levelset.phi(i-xradius,j-yradius,k-zradius) << std::endl;
                }
            }
        }
        LEVELSET_PROCESSOR<T> proc(levelset2);
        proc.Dilation(xradius,yradius,zradius);
    }
    virtual bool GrowBoundingBox() { return true; }
};

template<class T,class RW>
struct COMMAND_CLOSING:public COMMAND_KERNEL_OP<T,RW>
{
    typedef COMMAND_KERNEL_OP<T,RW> BASE;using BASE::xradius;using BASE::yradius;using BASE::zradius;

    virtual std::string Name() { return "closing"; }
    virtual std::string Description() { return "do closing."; }
    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        LEVELSET_PROCESSOR<T> proc(levelset);
        proc.Closing(xradius,yradius,zradius);
    }
};

template<class T,class RW>
struct COMMAND_OPENING:public COMMAND_KERNEL_OP<T,RW>
{
    typedef COMMAND_KERNEL_OP<T,RW> BASE;using BASE::xradius;using BASE::yradius;using BASE::zradius;

    virtual std::string Name() { return "opening"; }
    virtual std::string Description() { return "do opening."; }
    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        LEVELSET_PROCESSOR<T> proc(levelset);
        proc.Closing(xradius,yradius,zradius);
    }
};

template<class T,class RW>
struct COMMAND_SMOOTH:public COMMAND_KERNEL_OP<T,RW>
{
    typedef COMMAND_KERNEL_OP<T,RW> BASE;using BASE::xradius;using BASE::yradius;using BASE::zradius;

    virtual std::string Name() { return "smooth"; }
    virtual std::string Description() { return "smooth the mesh using erosion/dilation."; }
    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        LEVELSET_PROCESSOR<T> proc(levelset);
        proc.Smooth(xradius,yradius,zradius);
    }
};

template<class T,class RW>
struct COMMAND_ADD_SET:public COMMAND_LEVELSET_BINARY_OP<T,RW>
{
    virtual std::string Name() { return "add_set"; }
    virtual std::string Description() { return "add 2 sets together."; }
    void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset1,const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        LEVELSET_PROCESSOR<T> proc(levelset1);
        if(levelset1.grid == levelset2.grid)
            proc.Add_Set_Same_Size(levelset2);
        else proc.Add_Set(levelset2);
    }
};

template<class T,class RW>
struct COMMAND_SUBTRACT_SET:public COMMAND_LEVELSET_BINARY_OP<T,RW>
{
    virtual std::string Name() { return "subtract_set"; }
    virtual std::string Description() { return "add 2 sets together."; }
    void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset1,const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        LEVELSET_PROCESSOR<T> proc(levelset1);
        if(levelset1.grid == levelset2.grid)
            proc.Subtract_Set_Same_Size(levelset2);
        else proc.Subtract_Set(levelset2);
    }
};

template<class T,class RW>
struct COMMAND_SUBSET:public COMMAND_LEVELSET_OP<T,RW>
{
    typedef COMMAND_LEVELSET_OP<T,RW> BASE;using BASE::Add_Int_Arg;

    int m1,m2,n1,n2,mn1,mn2;
    COMMAND_SUBSET()
        :m1(1),m2(1),n1(1),n2(1),mn1(1),mn2(1)
    {
        Add_Int_Arg("m1",m1);
        Add_Int_Arg("m2",m2);
        Add_Int_Arg("n1",n1);
        Add_Int_Arg("n2",n2);
        Add_Int_Arg("mn1",mn1);
        Add_Int_Arg("mn2",mn2);
    }
    virtual std::string Name() { return "subset"; }
    virtual std::string Description() { return "extracts a sub-grid of the levelset"; }
    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        GRID<VECTOR<T,3> > grid;
        ARRAYS<VECTOR<T,3> > phi;
        LEVELSET_3D<GRID<VECTOR<T,3> > > temp(grid,phi);
        LEVELSET_PROCESSOR<T> proc(levelset);
        proc.Output_Subset(temp,m1,m2,n1,n2,mn1,mn2);
        levelset.grid = grid;
        levelset.phi = phi;
    }
};  

template<class T,class RW>
struct COMMAND_SUBSAMPLE:public COMMAND_LEVELSET_OP<T,RW>
{
    typedef COMMAND_LEVELSET_OP<T,RW> BASE;using BASE::Add_Int_Arg;

    int xstride,ystride,zstride;
    COMMAND_SUBSAMPLE()
        :xstride(2),ystride(2),zstride(2)
    {
        Add_Int_Arg("xstride",xstride);
        Add_Int_Arg("ystride",ystride);
        Add_Int_Arg("zstride",zstride);
    }
    virtual std::string Name() { return "subsample"; }
    virtual std::string Description() { return "subsamples a levelset by skipping cells."; }
    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        GRID<VECTOR<T,3> > grid;
        ARRAYS<VECTOR<T,3> > phi;
        LEVELSET_3D<GRID<VECTOR<T,3> > > temp(grid,phi);
        LEVELSET_PROCESSOR<T> proc(levelset);
        proc.Output_Subsample(temp,xstride,ystride,zstride,true);
        levelset.grid = grid;
        levelset.phi = phi;
    }
};  

template<class T,class RW>
struct COMMAND_RESAMPLE:public COMMAND_LEVELSET_OP<T,RW>
{
    typedef COMMAND_LEVELSET_OP<T,RW> BASE;using BASE::Add_Int_Arg;

    int m,n,mn;
    COMMAND_RESAMPLE()
        :m(1),n(1),mn(1)
    {
        Add_Int_Arg("m",m);
        Add_Int_Arg("n",n);
        Add_Int_Arg("mn",mn);
    }
    virtual std::string Name() { return "resample"; }
    virtual std::string Description() { return "resamples a levelset at the given resolution."; }
    virtual void DoOp(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset2) {
        GRID<VECTOR<T,3> > grid;
        ARRAYS<VECTOR<T,3> > phi;
        LEVELSET_3D<GRID<VECTOR<T,3> > > temp(grid,phi);
        LEVELSET_PROCESSOR<T> proc(levelset);
        proc.Output_Resample(temp,m,n,mn);
        levelset.grid = grid;
        levelset.phi = phi;
    }
};  

     // else if(!strcmp(control_unit.key_string,"-flat")){processor.Reinitialize_Levelset_Flat(processor.grid.dx/2);default_descriptive_output_filename+="_flat";}

template<class T,class RW>
struct COMMAND_PHI2TRI:public COMMAND_BASE
{
    int m,n,mn;

    COMMAND_PHI2TRI():m(-1),n(-1),mn(-1) {}
    virtual std::string Name() { return "phi2tri"; }
    virtual std::string Description() { return "create tri on default grid."; }
    virtual int Min_Inputs() { return 1; }
    virtual int Max_Inputs() { return 1; }
    virtual int Num_Outputs() { return 1; }

    virtual int Do()
    {
        GRID<VECTOR<T,3> > grid;
        ARRAYS<VECTOR<T,3> > phi;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset(grid,phi);
        FILE_UTILITIES::Read_From_File<RW>(Get_Input(1),levelset);
        TRIANGULATED_SURFACE<T>* surf = TRIANGULATED_SURFACE<T>::Create();
        if(m < 0) {
            levelset.Calculate_Triangulated_Surface_From_Marching_Tetrahedra(*surf);
        }
        else {
            GRID<VECTOR<T,3> > triGrid(m,n,mn,grid.Domain());
            levelset.Calculate_Triangulated_Surface_From_Marching_Tetrahedra(triGrid,*surf);
        }
        FILE_UTILITIES::Write_To_File<RW>(Get_Output(),*surf);
        delete surf;
        return 1;
    }
};

template<class T,class RW>
struct COMMAND_PHI2TRI_GRID:public COMMAND_PHI2TRI<T,RW>
{
    typedef COMMAND_PHI2TRI<T,RW> BASE;using BASE::Add_Int_Arg;using BASE::m;using BASE::n;using BASE::mn;

    COMMAND_PHI2TRI_GRID() {
        Add_Int_Arg("m",m);
        Add_Int_Arg("n",n);
        Add_Int_Arg("mm",mn);
    }
    virtual std::string Name() { return "phi2tri_grid"; }
    virtual std::string Description() { return "create tri on custom grid."; }
};

template<class T,class RW>
struct COMMAND_TRI2PHI:public COMMAND_BASE
{
    int m,n,mn;
    T dx,dy,dz;
    bool use_divs;
    bool use_heaviside;

    COMMAND_TRI2PHI(bool divs,bool heaviside)
        :m(10),n(10),mn(10),dx(1),dy(1),dz(1),use_divs(divs),use_heaviside(heaviside)
    {
        if(use_divs) {
            Add_Int_Arg("m",m);
            Add_Int_Arg("n",n);
            Add_Int_Arg("mn",mn);
        } else {
            Add_Float_Arg("dx",dx);
            Add_Float_Arg("dy",dy);
            Add_Float_Arg("dz",dz);
        }
    }
    virtual std::string Description() { return use_divs ? "grid divisions for triangulation (-1 for default)." : "grid cell size."; }
    virtual int Min_Inputs() { return 1; }
    virtual int Max_Inputs() { return 1; }
    virtual int Num_Outputs() { return 1; }

    virtual int Do()
    {
        GRID<VECTOR<T,3> > grid;
        ARRAYS<VECTOR<T,3> > phi;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset(grid,phi);

        TRIANGULATED_SURFACE<T>* surf = TRIANGULATED_SURFACE<T>::Create();
        TETRAHEDRON_MESH tet_mesh;

        DEFORMABLE_GEOMETRY_PARTICLES<VECTOR<T,3> > tet_particles;
        TETRAHEDRALIZED_VOLUME<T> tet_vol(tet_mesh,tet_particles);
        std::string fn = Get_Input(1);
        if(FILE_UTILITIES::File_Extension_Matches(fn,".tri"))
            FILE_UTILITIES::Read_From_File<RW>(fn,*surf);
        else if(FILE_UTILITIES::File_Extension_Matches(fn,".tet")){ // ooh we can handle tet files now
            std::cout<<"Loading a .tet file!"<<std::endl;
            FILE_UTILITIES::Read_From_File<RW>(fn,tet_vol);
            tet_vol.Initialize_Triangulated_Surface();
            surf = tet_vol.triangulated_surface;}

        surf->Update_Bounding_Box();
        //setup the grid
        if(use_divs)
            grid = Padded_Bounding_Grid(m,n,mn,*surf->bounding_box);
        else
            grid = Padded_Bounding_Grid_By_Size(dx,dy,dz,*surf->bounding_box);

        std::cout<<"Grid size: "<<grid.counts<<std::endl;
        std::cout<<"Grid resolution: "<<grid.dX<<std::endl;

        if(use_heaviside) {
            std::cout<<"Calculating heaviside function..."<<std::endl;
            surf->Calculate_Heaviside_Function(grid,phi,true);
            std::cout<<"Done calculating heaviside function."<<std::endl;
        }
        else {
            std::cout<<"Calculating signed distance function..."<<std::endl;
            surf->Calculate_Signed_Distance_Function(grid,phi,true);
            std::cout<<"done."<<std::endl;
        }
        delete surf;
        tet_vol.triangulated_surface = 0;

        FILE_UTILITIES::Write_To_File<RW>(Get_Output(),levelset);
        return 1;
    }
};

template<class T,class RW>
struct COMMAND_CREATE_HEAVISIDE:public COMMAND_TRI2PHI<T,RW>
{
    COMMAND_CREATE_HEAVISIDE():COMMAND_TRI2PHI<T,RW>(true,true) {}
    virtual std::string Name() { return "create_heaviside"; }
};

template<class T,class RW>
struct COMMAND_CREATE_DISTANCE:public COMMAND_TRI2PHI<T,RW>
{
    COMMAND_CREATE_DISTANCE():COMMAND_TRI2PHI<T,RW>(true,false) {}
    virtual std::string Name() { return "create_distance"; }
};

template<class T,class RW>
struct COMMAND_CREATE_HEAVISIDE_SIZE:public COMMAND_TRI2PHI<T,RW>
{
    COMMAND_CREATE_HEAVISIDE_SIZE():COMMAND_TRI2PHI<T,RW>(false,true) {}
    virtual std::string Name() { return "create_heaviside_size"; }
};

template<class T,class RW>
struct COMMAND_CREATE_DISTANCE_SIZE:public COMMAND_TRI2PHI<T,RW>
{
    COMMAND_CREATE_DISTANCE_SIZE():COMMAND_TRI2PHI<T,RW>(false,false) {}
    virtual std::string Name() { return "create_distance_size"; }
};


template<class T,class RW>
struct COMMAND_UNION:public COMMAND_BASE
{
    int m,n,mn;
    bool heaviside;

    COMMAND_UNION():m(-1),n(-1),mn(-1),heaviside(false) {}
    virtual std::string Name() { return "union"; }
    virtual std::string Description() { return "unions level sets together."; }
    virtual int Min_Inputs() { return 2; }
    virtual int Max_Inputs() { return -1; }
    virtual int Num_Outputs() { return 1; }
    virtual int Do()
    {
        std::cout<<"Begin union..."<<std::endl;
        LEVELSET_CSG<T> csg(heaviside);
        if(!csg.Union_Files(inputs)) return 0;
        std::cout<<"Done, writing to "<<Get_Output()<<std::endl;
        FILE_UTILITIES::Write_To_File<RW>(Get_Output(),csg.levelset_full);
        return 1;
    }
};

template<class T,class RW>
struct COMMAND_UNION_HEAVISIDE:public COMMAND_UNION<T,RW>
{
    using COMMAND_UNION<T,RW>::heaviside;

    COMMAND_UNION_HEAVISIDE() { heaviside = true; }
    virtual std::string Name() { return "union_heaviside"; }
    virtual std::string Description() { return "unions heaviside level sets together (faster)."; }
};

template<class T,class RW>
struct COMMAND_CSG_SUBTRACT:public COMMAND_BASE
{
    int m,n,mn;
    bool heaviside;

    COMMAND_CSG_SUBTRACT():m(-1),n(-1),mn(-1),heaviside(false) {}
    virtual std::string Name() { return "csg_subtract"; }
    virtual std::string Description() { return "performs csg subtraction on level sets."; }
    virtual int Min_Inputs() { return 2; }
    virtual int Max_Inputs() { return -1; }
    virtual int Num_Outputs() { return 1; }
    virtual int Do()
    {
        std::cout<<"Begin csg subtraction..."<<std::endl;
        LEVELSET_CSG<T> csg(heaviside);
        if(!csg.Subtract_Files(inputs)) return 0;
        std::cout<<"Done, writing to "<<Get_Output()<<std::endl;
        FILE_UTILITIES::Write_To_File<RW>(Get_Output(),csg.levelset_full);
        return 1;
    }
};

template<class T,class RW>
struct COMMAND_CSG_SUBTRACT_HEAVISIDE:public COMMAND_CSG_SUBTRACT<T,RW>
{
    using COMMAND_CSG_SUBTRACT<T,RW>::heaviside;

    COMMAND_CSG_SUBTRACT_HEAVISIDE() { heaviside = true; }
    virtual std::string Name() { return "csg_subtract_heaviside"; }
    virtual std::string Description() { return "performs csg subtraction on heaviside level sets (faster)."; }
};

template<class T,class RW>
struct COMMAND_INTERSECT:public COMMAND_BASE
{
    int m,n,mn;
    bool heaviside;

    COMMAND_INTERSECT():m(-1),n(-1),mn(-1),heaviside(false) {}
    virtual std::string Name() { return "intersect"; }
    virtual std::string Description() { return "intersects level sets together."; }
    virtual int Min_Inputs() { return 2; }
    virtual int Max_Inputs() { return -1; }
    virtual int Num_Outputs() { return 1; }
    virtual int Do()
    {
        std::cout<<"Begin intersection..."<<std::endl;
        LEVELSET_CSG<T> csg(heaviside);
        if(!csg.Intersect_Files(inputs)) return 0;
        std::cout<<"Done, writing to "<<Get_Output()<<std::endl;
        FILE_UTILITIES::Write_To_File<RW>(Get_Output(),csg.levelset_full);
        return 1;
    }
};

template<class T,class RW>
struct COMMAND_INTERSECT_HEAVISIDE:public COMMAND_INTERSECT<T,RW>
{
    using COMMAND_INTERSECT<T,RW>::heaviside;

    COMMAND_INTERSECT_HEAVISIDE() { heaviside = true; }
    virtual std::string Name() { return "intersect_heaviside"; }
    virtual std::string Description() { return "interesects heaviside level sets together (faster)."; }
};

template<class T,class RW>
struct COMMAND_FLAT:public COMMAND_BASE
{
    COMMAND_FLAT() {}
    virtual std::string Name() { return "flat"; }
    virtual std::string Description() { return "convert to heaviside function."; }
    virtual int Min_Inputs() { return 1; }
    virtual int Max_Inputs() { return 1; }
    virtual int Num_Outputs() { return 1; }
    virtual int Do()
    {
        std::cout << "Begin flattening..."<<std::endl;
        GRID<VECTOR<T,3> > grid;
        ARRAYS<VECTOR<T,3> > phi;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset_temp(grid, phi);
        FILE_UTILITIES::Read_From_File<RW>(inputs(1),levelset_temp);
        T phi_value = levelset_temp.grid.dX.x/2;
        for(int i=levelset_temp.phi.m_start;i<=levelset_temp.phi.m_end;i++) {
            for(int j=levelset_temp.phi.n_start;j<=levelset_temp.phi.n_end;j++) {
                for(int k=levelset_temp.phi.mn_start;k<=levelset_temp.phi.mn_end;k++) {
                    levelset_temp.phi(i,j,k) = (levelset_temp.phi(i,j,k) >= 0) ? phi_value : -phi_value;
                }}}
        std::cout<<"Done, writing to "<<Get_Output()<<std::endl;
        FILE_UTILITIES::Write_To_File<RW>(Get_Output(),levelset_temp);
        return 1;
    }
};

template<class T,class RW>
struct COMMAND_MAC_GRID:public COMMAND_BASE
{
    COMMAND_MAC_GRID() {}
    virtual int Min_Inputs() { return 1; }
    virtual int Max_Inputs() { return 1; }
    virtual int Num_Outputs() { return 1; }
    virtual std::string Name() { return "mac_grid"; }
    virtual std::string Description() { return "creates a mac grid from a uniform grid."; }
    virtual int Do()
    {
        std::cout << "Begin conversion..."<<std::endl;
        GRID<VECTOR<T,3> > grid, mac_grid;
        ARRAYS<VECTOR<T,3> > phi;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset_temp(grid,phi);
        FILE_UTILITIES::Read_From_File<RW>(inputs(1),levelset_temp);
        mac_grid=grid.Get_MAC_Grid();
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset_mac(mac_grid,phi);
        std::cout<<"Done, writing to "<<Get_Output()<<std::endl;
        FILE_UTILITIES::Write_To_File<RW>(Get_Output(),levelset_mac);
        return 1;
    }
};

template<class T,class RW>
struct COMMAND_THICKEN:public COMMAND_BASE
{
    int sm,sn,smn;
    COMMAND_THICKEN():sm(0),sn(0),smn(0) { Add_Int_Arg("sm",sm);Add_Int_Arg("sn",sn);Add_Int_Arg("smn",smn);}
    virtual int Min_Inputs() { return 1; }
    virtual int Max_Inputs() { return 1; }
    virtual int Num_Outputs() { return 1; }
    virtual std::string Name() { return "thicken"; }
    virtual std::string Description() { return "Thickens the grid by the input."; }
    virtual int Do()
    {
        std::cout<<"Begin thickening..."<<std::endl;
        GRID<VECTOR<T,3> > grid;
        ARRAYS<VECTOR<T,3> > phi;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset(grid,phi);
        FILE_UTILITIES::Read_From_File<RW>(inputs(1),levelset);
        GRID<VECTOR<T,3> > grid_big(phi.m_end+sm*2,phi.n_end+sn*2,phi.mn_end+smn,grid.domain.min_corner.x-(T)sm*grid.dX.x,grid.domain.max_corner.x+(T)sm*grid.dX.x,grid.domain.min_corner.y-(T)sn*grid.dX.y,grid.domain.max_corner.y+(T)sn*grid.dX.y,grid.domain.min_corner.z,grid.domain.max_corner.z+(T)smn*grid.dX.z);
        //GRID<VECTOR<T,3> > grid_big(phi.m_end+sm*2,phi.n_end+sn*2,phi.mn_end+smn*2,grid.domain.min_corner.x-(T)sm*grid.dX.x,grid.domain.max_corner.x+(T)sm*grid.dX.x,grid.domain.min_corner.y-(T)sn*grid.dX.y,grid.domain.max_corner.y+(T)sn*grid.dX.y,grid.domain.min_corner.z-(T)smn*grid.dX.z,grid.domain.max_corner.z+(T)smn*grid.dX.z);
        ARRAYS<VECTOR<T,3> > phi_big(grid_big.Domain_Indices());
        for(int i=phi.m_start;i<=phi.m_end;i++) for(int j=phi.n_start;j<=phi.n_end;j++) for(int ij=phi.mn_start;ij<=phi.mn_end;ij++) phi_big(i+sm,j+sn,ij+smn)=phi(i,j,ij);
        for(int s=0;s<sm;s++){
            for(int i=phi_big.n_start;i<=phi_big.n_end;i++) for(int j=phi_big.mn_start;j<=phi_big.mn_end;j++) {phi_big(phi_big.m_start+s,i,j)=grid.min_dX;phi_big(phi_big.m_end-s,i,j)=grid.min_dX;}}
        for(int s=0;s<sn;s++){
            for(int i=phi_big.m_start;i<=phi_big.m_end;i++) for(int j=phi_big.mn_start;j<=phi_big.mn_end;j++) {phi_big(i,phi_big.n_start+s,j)=grid.min_dX;phi_big(i,phi_big.n_end-s,j)=grid.min_dX;}}
        for(int s=0;s<smn;s++){
            for(int i=phi_big.n_start;i<=phi_big.n_end;i++) for(int j=phi_big.m_start;j<=phi_big.m_end;j++) {phi_big(j,i,phi_big.mn_start)=grid.min_dX;phi_big(j,i,phi_big.mn_end-s)=grid.min_dX;}}
            //for(int i=phi_big.n_start;i<=phi_big.n_end;i++) for(int j=phi_big.m_start;j<=phi_big.m_end;j++) {phi_big(j,i,phi_big.mn_start+s)=grid.min_dX;phi_big(j,i,phi_big.mn_end-s)=grid.min_dX;}}
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset_big(grid_big,phi_big);
        std::cout<<"Done, writing to "<<Get_Output()<<std::endl;
        FILE_UTILITIES::Write_To_File<RW>(Get_Output(),levelset_big);
        return 1;
    }
};

COMMAND_BASE* commands[] = {
    new COMMAND_PHI2TRI<float,float>,
    new COMMAND_PHI2TRI_GRID<float,float>,
    new COMMAND_CREATE_HEAVISIDE<float,float>,
    new COMMAND_CREATE_DISTANCE<float,float>,
    new COMMAND_CREATE_HEAVISIDE_SIZE<float,float>,
    new COMMAND_CREATE_DISTANCE_SIZE<float,float>,
    new COMMAND_EROSION<float,float>,
    new COMMAND_DILATION<float,float>,
    new COMMAND_CLOSING<float,float>,
    new COMMAND_OPENING<float,float>,
    new COMMAND_SMOOTH<float,float>,
    new COMMAND_ADD_SET<float,float>,
    new COMMAND_SUBTRACT_SET<float,float>,
    new COMMAND_SUBSET<float,float>,
    new COMMAND_SUBSAMPLE<float,float>,
    new COMMAND_RESAMPLE<float,float>,
    new COMMAND_MOTION_MEAN_CURVATURE<float,float>,
    new COMMAND_FAST_MARCHING<float,float>,
    new COMMAND_ADD_CONSTANT<float,float>,
    new COMMAND_UNION<float,float>,
    new COMMAND_UNION_HEAVISIDE<float,float>,
    new COMMAND_CSG_SUBTRACT<float,float>,
    new COMMAND_CSG_SUBTRACT_HEAVISIDE<float,float>,
    new COMMAND_INTERSECT<float,float>,
    new COMMAND_INTERSECT_HEAVISIDE<float,float>,
    new COMMAND_FLAT<float,float>,
    new COMMAND_MAC_GRID<float,float>,
    new COMMAND_THICKEN<float,float>};
#define ARRAYSIZE(x) (sizeof(x)/sizeof(x[0]));
int num_commands = ARRAYSIZE(commands);

int main(int argc,const char** argv)
{
    int res=Command_Reader(commands,num_commands,argc,argv);
    if(res < 0) {
        std::cout << "Press any key to continue."<<std::endl;
        getchar();
    }
    return res;
}
