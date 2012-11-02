//#####################################################################
// Copyright 2003-2005, Jiayi Chong, Eran Guendelman, Cynthia Lau, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
using namespace PhysBAM;
//#####################################################################

template <class T>
GRID<VECTOR<T,3> > Padded_Bounding_Grid_By_Size(T dx,T dy,T dz,const RANGE<VECTOR<T,3> >& surfBox,T pad=1);

template <class T>
T csgsub(T a, T b) {return (a >= 0 && b > 0) ? b : -a;}

template<class T>
T Min(T a,T b)
{return std::min(a,b);}

template<class T>
T Max(T a,T b)
{return std::max(a,b);}

template <class T>
struct LEVELSET_CSG
{
    GRID<VECTOR<T,3> > grid_full;
    ARRAYS<VECTOR<T,3> > phi_full;
    LEVELSET_3D<GRID<VECTOR<T,3> > > levelset_full;
    bool heaviside;

    const static T big;

    RANGE<VECTOR<T,3> > box;
    int m,n,mn;
    VECTOR<T,3> min_cell;

    LEVELSET_CSG(bool heaviside_input)
        :levelset_full(grid_full,phi_full),heaviside(heaviside_input),m(-1),n(-1),mn(-1),min_cell(big,big,big)
    {}

    bool Union_Files(const ARRAY<std::string>& files) {
        assert(files.m >= 1);
        GRID<VECTOR<T,3> > grid_temp;
        ARRAYS<VECTOR<T,3> > phi_temp;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset_temp(grid_temp,phi_temp);

        //load the grids -- this depends on having
        FILE_UTILITIES::Read_From_File<T>(files(1),grid_temp);
        Start_Domain(grid_temp);
        for(int i=2;i<=files.m;i++){
            FILE_UTILITIES::Read_From_File<T>(files(i),grid_temp);
            Add_Domain(grid_temp);}
        Finish_Domain();

        for(int i=0;i<files.m;i++){
            FILE_UTILITIES::Read_From_File<T>(files(i),levelset_temp);
            Union(levelset_temp);}
        return true;
    }

    bool Subtract_Files(const ARRAY<std::string>& files) {
        assert(files.m >= 1);
        GRID<VECTOR<T,3> > grid_temp;
        ARRAYS<VECTOR<T,3> > phi_temp;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset_temp(grid_temp,phi_temp);

        //load the grids -- this depends on having
        FILE_UTILITIES::Read_From_File<T>(files(1),grid_temp);
        Start_Domain(grid_temp);
        for(int i=2;i<=files.m;i++){
            FILE_UTILITIES::Read_From_File<T>(files(i),grid_temp);
            Add_Domain(grid_temp);}
        Finish_Domain();

        FILE_UTILITIES::Read_From_File<T>(files(1),levelset_temp);
        Union(levelset_temp);
        for(int i=2;i<=files.m;i++){
            FILE_UTILITIES::Read_From_File<T>(files(i),levelset_temp);
            Subtract(levelset_temp);}
        return true;
    }
    bool Intersect_Files(const ARRAY<std::string>& files) {
        assert(files.m >= 1);
        GRID<VECTOR<T,3> > grid_temp;
        ARRAYS<VECTOR<T,3> > phi_temp;
        LEVELSET_3D<GRID<VECTOR<T,3> > > levelset_temp(grid_temp,phi_temp);

        //load the grids -- this depends on having
        FILE_UTILITIES::Read_From_File<T>(files(1),grid_temp);
        Start_Domain(grid_temp);
        for(int i=2;i<=files.m;i++){
            FILE_UTILITIES::Read_From_File<T>(files(i),grid_temp);
            Add_Domain(grid_temp);}
        Finish_Domain();

        FILE_UTILITIES::Read_From_File<T>(files(1),levelset_temp);
        Union(levelset_temp);
        for(int i=2;i<=files.m;i++) {
            FILE_UTILITIES::Read_From_File<T>(files(i),levelset_temp);
            Intersect(levelset_temp);
        }
        return true;
    }
    //default resolution calculates it from the minimum cell
    void Set_Resolution(const int m_input, const int n_input, const int mn_input)
    { m=m_input; n=n_input; mn=mn_input; }

    //first the domain must be set, either by Set... or {Start...Add...Finish...}
    void Set_Domain(const GRID<VECTOR<T,3> >* grids, int num_grids) {
        //calculate bounding box of all level sets
        Start_Domain(grids[0]);
        for(int i=1;i<num_grids;i++)
            Add_Domain(grids[i]);
        Finish_Domain();
    }
    void Set_Domain(const LEVELSET_3D<GRID<VECTOR<T,3> > >* levelsets, int num_levelsets);

    void Start_Domain(const GRID<VECTOR<T,3> >& grid) {
        box = grid.domain;
        min_cell = grid.dX;
    }
    void Start_Domain(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset) { Start_Domain(levelset.grid); }

    void Add_Domain(const GRID<VECTOR<T,3> >& grid) {
        box.Enlarge_To_Include_Box(grid.Domain());
        min_cell = VECTOR<T,3>::Componentwise_Min(min_cell,grid.dX);
    }
    void Add_Domain(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset) { Add_Domain(levelset.grid); }

    void Finish_Domain() {
        if(m < 0 || n < 0 || mn < 0) {  //calculate default divisions
            grid_full = Padded_Bounding_Grid_By_Size<float>(min_cell.x,min_cell.y,min_cell.z,box,0.0);
            m=grid_full.counts.x;n=grid_full.counts.y;mn=grid_full.counts.z;
        }
        else grid_full.Initialize(m,n,mn,box);

        std::cout<<"Bounding box: "<<box<<", divisions: "<<m<<" "<<n<<" "<<mn<<std::endl;

        phi_full.Resize(1,m,1,n,1,mn);
        T init_value = heaviside ? 1 : big;
        phi_full.array.Fill(init_value);
    }

    void CSG_Op(const LEVELSET_3D<GRID<VECTOR<T,3> > >& unioned, T (*op) (T a, T b)) {
        if(heaviside) {
            //we can take some shortcuts - only consider the domain of the unioned
            int m0,m1,n0,n1,mn0,mn1;        //extents of the unioned object's box
            VECTOR<T,3> minpos(unioned.grid.domain.min_corner),maxpos(unioned.grid.domain.max_corner);
            VECTOR<int,3> minindex=grid_full.Clamped_Index(minpos),maxindex=grid_full.Clamped_Index(maxpos);
            m0=minindex.x; m1=maxindex.x; n0=minindex.y;n1=maxindex.y; mn0=minindex.z; mn1=maxindex.z;
            std::cout<<"range "<<m0<<" "<<m1<<" "<<n0<<" "<<n1<<" "<<mn0<<" "<<mn1<<std::endl;

            VECTOR<T,3> pos;
            for(int i=m0;i<=m1;i++) { for(int j=n0;j<=n1;j++) { for(int k=mn0;k<=mn1;k++) {
                pos=grid_full.X(i,j,k);
                T phi1 = unioned.Phi(pos);
                T phi2 = levelset_full.phi(i,j,k);
                levelset_full.phi(i,j,k) = op(phi1,phi2);
            }}}
        }
        else {
            VECTOR<T,3> pos;
            for(int i=0;i<grid_full.counts.x;i++) { for(int j=0;j<grid_full.counts.y;j++) { for(int k=0;k<grid_full.counts.z;k++) {
                pos=grid_full.X(i,j,k);
                T phi1 = unioned.Extended_Phi(pos);
                T phi2 = levelset_full.phi(i,j,k);
                levelset_full.phi(i,j,k) = op(phi1,phi2);
            }}}
        }
    }
    void Union(const LEVELSET_3D<GRID<VECTOR<T,3> > >& unioned) {
        std::cout<<"Union..."<<std::endl;
        CSG_Op(unioned,Min);
        std::cout << "Done." << std::endl;
    }
    void Subtract(const LEVELSET_3D<GRID<VECTOR<T,3> > >& unioned) {
        std::cout<<"Subtract..."<<std::endl;
        CSG_Op(unioned,csgsub);
        std::cout << "Done." << std::endl;
    }
    void Intersect(const LEVELSET_3D<GRID<VECTOR<T,3> > >& unioned) {
        std::cout<<"Intersect..."<<std::endl;
        CSG_Op(unioned,Max);
        std::cout << "Done." << std::endl;
    }
};

template <class T>
const T LEVELSET_CSG<T>::big = 2e6;



template <class T>
GRID<VECTOR<T,3> > Padded_Bounding_Grid(int m,int n,int mn, const RANGE<VECTOR<T,3> >& surfBox,T pad=1)
{   assert(m > 0 && n > 0 && mn > 0);
    RANGE<VECTOR<T,3> > bbox = surfBox;
    //figure out how much padding we should do
    T dx = 1.0/(T)m;
    T dy = 1.0/(T)n;
    T dz = 1.0/(T)mn;
    bbox.Scale_About_Center(1.0+2.0*pad*dx,1.0+2.0*pad*dy,1.0+2.0*pad*dz);
    int extra=(int)ceil(2*pad);
    return GRID<VECTOR<T,3> >(m+extra,n+extra,mn+extra,bbox);
}

template <class T>
GRID<VECTOR<T,3> > Padded_Bounding_Grid_By_Size(T dx,T dy,T dz,const RANGE<VECTOR<T,3> >& surfBox,T pad)
{   assert(dx > 0 && dy > 0 && dz > 0);
    //create a slightly larger bbox than necessary
    T divx,divy,divz;  //number of divisions required
    int m,n,mn;     //# of grid cells on each side
    divx = (T)ceil((surfBox.max_corner.x-surfBox.min_corner.x)/dx);
    divy = (T)ceil((surfBox.max_corner.y-surfBox.min_corner.y)/dy);
    divz = (T)ceil((surfBox.max_corner.z-surfBox.min_corner.z)/dz);
    std::cout << divx << " " << divy << " " << divz<< std::endl;
    m = (int)divx;
    n = (int)divy;
    mn = (int)divz;
    //add a padding per side
    RANGE<VECTOR<float,3> > bbox(-pad,divx+pad,-pad,divy+pad,-pad,divz+pad);
    bbox.min_corner.x*=dx; bbox.max_corner.x*=dx;
    bbox.min_corner.y*=dy; bbox.max_corner.y*=dy;
    bbox.min_corner.z*=dz; bbox.max_corner.z*=dz;
    //offset the box by the "lower corner"
    bbox.min_corner+=surfBox.min_corner;bbox.max_corner+=surfBox.min_corner;
    int extra = 1+int(2.0*pad); //+1 for the edges, +1 for the padding
    return GRID<VECTOR<T,3> >(m+extra,n+extra,mn+extra,bbox); 
}

