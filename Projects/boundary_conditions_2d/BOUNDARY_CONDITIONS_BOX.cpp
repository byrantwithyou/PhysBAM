#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <climits>
#include <cstdio>
#include "BOUNDARY_CONDITIONS_BOX.h"
#include "HEADER.h"

template<class TV> BOUNDARY_CONDITIONS_BOX<TV>::
BOUNDARY_CONDITIONS_BOX(GRID<TV>& grid_input,T offset,T angle,std::string& bc_types)
    :BOUNDARY_CONDITIONS<TV>(grid_input)
{
    boundary_gap=1e-6;
    PHYSBAM_ASSERT(!offset || !angle);
    base_domain=RANGE<TV>::Unit_Box();
    if(angle) Set_Planes_From_Angle_And_Box(base_domain,angle);
    else Set_Planes_From_Offset_And_Box(base_domain,offset);
    for(int i=0;i<4;i++) bounding_box.Enlarge_To_Include_Point(corners[i]);
    for(int i=0;i<4;i++) poly.X.Append(corners[i]);

    int bc_translate[CHAR_MAX]={0};
    bc_translate[(int)'s']=BOUNDARY_CONDITIONS_BOX<TV>::source;
    bc_translate[(int)'w']=BOUNDARY_CONDITIONS<TV>::noslip;
    bc_translate[(int)'f']=BOUNDARY_CONDITIONS<TV>::free;

    type[0]=bc_translate[(int)bc_types[0]];
    type[1]=bc_translate[(int)bc_types[1]];
    type[2]=bc_translate[(int)bc_types[2]];
    type[3]=bc_translate[(int)bc_types[3]];

}

template<class TV> BOUNDARY_CONDITIONS_BOX<TV>::
~BOUNDARY_CONDITIONS_BOX()
{
    int n=0;
    T L1=0,Li=0;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,2);it.Valid();it.Next()){
        TV X(grid.Axis_X_Face(it.Full_Index()));
        if(Theta(X)>grid.dX.Max()*2) continue;
        if(!Inside(X)) continue;
        T e=bounding_box.min_corner.y;
        TV Y(X.y*(6*e+(3-3*e)*X.y-2*X.y*X.y),X.x*X.x*(3-2*X.x));

        
        Y=(X-bounding_box.Center()).Orthogonal_Vector()+(X-bounding_box.Center())*.0625;
        Y=cos(X)+sin((T)4*X).Orthogonal_Vector();

        MATRIX<T,2> J(-sin(X.x),-4*cos(4*X.y),4*cos(4*X.x),-sin(X.y));
        Y-=((T)0.0625)*J.Transposed()*Y;

        n++;
        T q=abs((*saved_u)(it.Full_Index())-Y(it.Axis()));
        TV dX=X-bounding_box.Center();
        printf("%.15g %.15g -> %.15g %.15g %.15g\n",(*saved_u)(it.Full_Index()),Y(it.Axis()),q,Theta(X),atan2(dX.y,dX.x));
        L1+=q;
        Li=max(Li,q);}

    LOG::cout<<"W@Z-L1 "<<grid.counts.x<<" "<<L1/n<<std::endl;
    LOG::cout<<"W@Z-Linf "<<grid.counts.x<<" "<<Li<<std::endl;
}

template<class TV> LINE_2D<typename TV::SCALAR> BOUNDARY_CONDITIONS_BOX<TV>::
Bounding_Edge_From_Endpoints(const TV& p1,const TV& p2) const
{return LINE_2D<T>((p1-p2).Normalized().Orthogonal_Vector(),(p1+p2)/2);}

template<class TV> void BOUNDARY_CONDITIONS_BOX<TV>::
Set_Planes_From_Corners()
{for(int i=0;i<4;i++) planes[i]=Bounding_Edge_From_Endpoints(corners[i],corners[i+1]);}

template<class TV> void BOUNDARY_CONDITIONS_BOX<TV>::
Set_Corners_From_Box(const RANGE<TV>& box)
{
    corners[0]=TV(box.min_corner.x,box.max_corner.y);
    corners[1]=box.min_corner;
    corners[2]=TV(box.max_corner.x,box.min_corner.y);
    corners[3]=box.max_corner;
    corners[4]=corners[0];
}

template<class TV> void BOUNDARY_CONDITIONS_BOX<TV>::
Set_Planes_From_Angle_And_Box(const RANGE<TV>& box,T angle)
{
    Set_Corners_From_Box(box);
    T shift=tan(angle)*(box.max_corner.x-box.min_corner.x);
    corners[2].y+=shift;
    corners[3].y+=shift;
    Set_Planes_From_Corners();
}

template<class TV> void BOUNDARY_CONDITIONS_BOX<TV>::
Set_Planes_From_Offset_And_Box(const RANGE<TV>& box,T offset)
{
    Set_Corners_From_Box(box);
    corners[1].y-=offset;
    corners[2].y-=offset;
    Set_Planes_From_Corners();
}

template<class TV> typename TV::SCALAR BOUNDARY_CONDITIONS_BOX<TV>::
Theta(const TV& X) const
{
    return poly.Signed_Distance(X);
}

template<class TV> bool BOUNDARY_CONDITIONS_BOX<TV>::
Inside(const TV& X) const
{
    return poly.Signed_Distance(X)<=-boundary_gap;
}

template<class TV> int BOUNDARY_CONDITIONS_BOX<TV>::
Boundary_Condition(const TV& X,const TV& Y,T& theta,TV& value,T time) const
{
    for(int i=0;i<4;i++) if(planes[i].Signed_Distance(Y)>-boundary_gap){
        T py=planes[i].Signed_Distance(Y),px=planes[i].Signed_Distance(X);
        theta=px/(px-py);
        int tp=type[i];
        if(tp==source){
            int ns=BOUNDARY_CONDITIONS<TV>::noslip;
            tp=ns;
            TV p=X+theta*(Y-X);
            TV u=corners[i+1]-corners[i];
            T y=TV::Dot_Product(p-corners[i],u)/u.Magnitude_Squared();
            value=planes[i].normal*Source_Velocity_Curve(type[(i+3)%4]==ns,type[(i+1)%4]==ns,y);}
        else value=TV();
        return tp;}
    LOG::cout<<"BC FAIL ON "<<X<<"  "<<Y<<std::endl;
    return 0;
}

template<class TV> void BOUNDARY_CONDITIONS_BOX<TV>::
Initialize_Velocity_Field(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const
{
    saved_u=&u;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,2);it.Valid();it.Next()){
      if(this->check_leaks && Theta(grid.Axis_X_Face(it.Full_Index()))>boundary_gap){u(it.Full_Index())=1e10;continue;}
        TV X=it.Location();
        X-=bounding_box.Center();
//        X-=bounding_box.min_corner;
//        X/=bounding_box.Edge_Lengths();
//        T s=sin(X*pi).Product()*10;
//        X-=.5;
        TV Y=(X.Orthogonal_Vector()+X+(T).01*sin(X));

        T e=bounding_box.min_corner.y;
        X=it.Location();
        Y.x=X.y*(6*e+(3-3*e)*X.y-2*X.y*X.y);
        Y.y=X.x*X.x*(3-2*X.x);

        Y=cos(X)+sin((T)4*X).Orthogonal_Vector();

        u(it.Full_Index())=Y(it.Axis());}
}

template<class TV> TV BOUNDARY_CONDITIONS_BOX<TV>::
Analytic_Velocity(const TV& X,T time) const
{
    return TV();
}

template struct BOUNDARY_CONDITIONS_BOX<VECTOR<double,2> >;
