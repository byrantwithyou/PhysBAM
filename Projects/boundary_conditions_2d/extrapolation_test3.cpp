#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <iomanip>

using namespace PhysBAM;

template<class T> T f(const VECTOR<T,3>& r){return cos(r.x)*sin(r.y)*sin(r.z);}

int main(int argc,char* argv[])
{
    typedef double T;
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    LOG::cout<<std::setprecision(2*sizeof(T));

    int resolution=atoi(argv[1]);
    int iterations=atoi(argv[2]);
    int order=atoi(argv[3]);

    GRID<TV> grid(TV_INT(resolution,resolution,resolution),RANGE<TV>::Centered_Box()*(T)pi,true);
    ARRAY<T,FACE_INDEX<TV::m> > x(grid.Domain_Indices());
    ARRAY<T,TV_INT> p(grid.Domain_Indices());
    LEVELSET_3D<GRID<TV> > phi(grid,p);

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        TV X(it.Location());
        T r=X.Magnitude()-sqrt((T)2.0);
        p(it.index)=r;}

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        TV X(it.Location());
        T r=X.Magnitude()-sqrt((T)2.0);
        if(r<=0) x(it.Full_Index())=f(it.Location())+it.Axis();}

    EXTRAPOLATION_HIGHER_ORDER<TV,T>::Extrapolate_Face(grid,phi,0,x,iterations,order,3);

    T m=0;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        TV X(it.Location()),Y=X;
        T r=Y.Normalize(),c=sqrt((T)2.0);
        T e=x(it.Full_Index())-(f(X)+it.Axis());
        if(r-c<grid.dX.Max()*3) m=max(abs(e),m);}
    LOG::cout<<"L-inf "<<m<<std::endl;

    return 0;
}
