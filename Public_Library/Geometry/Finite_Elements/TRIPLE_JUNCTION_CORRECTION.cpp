#include <Tools/Arrays/ARRAY_VIEW.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND_VIEW.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Tools/Math_Tools/integer_log.h>
#include <Tools/Math_Tools/INTERVAL.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Finite_Elements/TRIPLE_JUNCTION_CORRECTION.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_TETRAHEDRA.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <climits>
#include <iomanip>
using namespace PhysBAM;
template<class T>
struct RAT
{
    T n0,n1,d0,d1;
    T x0,x1;
    bool in_x;

    T eval(T x) const {return (n1*x+n0)/(d1*x+d0);}
};

template<class T,class TV,class TV_INT>
void Emit_Rat(const GRID<TV>& grid,const RAT<T>& r,const TV_INT& cell,const VECTOR<T,3>& col)
{
    int N=50;
    T dx=(r.x1-r.x0)/N;
    for(int t=0;t<N;t++){
        T x0=r.x0+dx*t,x1=x0+dx,y0=r.eval(x0),y1=r.eval(x1);
        VECTOR<TV,2> pts=r.in_x?VECTOR<TV,2>(TV(x0,y0),TV(x1,y1)):VECTOR<TV,2>(TV(y0,x0),TV(y1,x1));
        pts=(pts+TV(cell))*grid.dX+grid.domain.min_corner;
        Add_Debug_Object(pts,col);}
}

template<class T,class TV,class TV_INT>
void Handle_Square(const GRID<TV>& grid,const VECTOR<T,8>& phi,const TV_INT& cell,const VECTOR<T,3>& col)
{
}

template<class T,class TV,class TV_INT>
void Handle_Square(const GRID<TV>& grid,const VECTOR<T,4>& phi,const TV_INT& cell,const VECTOR<T,3>& col)
{
    VECTOR<PAIR<T,T>,4> pts(PAIR<T,T>(2,2),PAIR<T,T>(2,2),PAIR<T,T>(2,2),PAIR<T,T>(2,2));
    if((phi(0)<0)!=(phi(1)<0)){pts(0).x=phi(0)/(phi(0)-phi(1));pts(0).y=0;}
    if((phi(2)<0)!=(phi(3)<0)){pts(1).x=phi(2)/(phi(2)-phi(3));pts(1).y=1;}
    if((phi(0)<0)!=(phi(2)<0)){pts(2).y=phi(0)/(phi(0)-phi(2));pts(2).x=0;}
    if((phi(1)<0)!=(phi(3)<0)){pts(3).y=phi(1)/(phi(1)-phi(3));pts(3).x=1;}
    pts.Sort();
    for(int st=0;st<4;st+=2){
        if(pts(st+1).x>1) break;
        if(pts(st+1).x-pts(st).x>abs(pts(st+1).y-pts(st).y)){
            RAT<T> rat={phi(0),phi(1)-phi(0),phi(0)-phi(2),phi(2)+phi(1)-phi(0)-phi(3),pts(st).x,pts(st+1).x,true};
            Emit_Rat(grid,rat,cell,col);}
        else{
            RAT<T> rat={phi(0),phi(2)-phi(0),phi(0)-phi(1),phi(1)+phi(2)-phi(0)-phi(3),pts(st).y,pts(st+1).y,false};
            Emit_Rat(grid,rat,cell,col);}}
}

template<class T,class TV,class TV_INT>
void Dump_Interface(const GRID<TV>& grid,ARRAY_VIEW<T,TV_INT> p,const VECTOR<T,3>& col)
{
    T default_phi=(T)1.12349871352389e-10;
    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());

    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        VECTOR<T,(1<<TV::m)> phi(p.Subset(bits+it.index));
        if(phi.Contains(default_phi)) continue;
        if(phi.Min()>=0 || phi.Max()<=0) continue;
        Handle_Square(grid,phi,it.index,col);}
}

//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIPLE_JUNCTION_CORRECTION<TV>::
TRIPLE_JUNCTION_CORRECTION(const GRID<TV>& grid,ARRAY<ARRAY<T,TV_INT> >& phi,int ghost)
    :grid(grid),phi(phi),ghost(ghost),default_phi((T)1.12349871352389e-10),trust_buffer(grid.dX.Max()*4),
    valid_width(ghost*grid.dX.Max()),extent(3*valid_width),extrap_width(2*ghost+3),bc_colors(3)
{
}
//#####################################################################
// Function Compute_Pairwise_Data
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Compute_Pairwise_Data()
{
    pairwise_data.Resize(grid.Node_Indices(ghost));

    for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        PAIRWISE_LEVEL_SET_DATA& data=pairwise_data(it.index);
        VECTOR<T,3> min(FLT_MAX,FLT_MAX,FLT_MAX);
        for(int i=0;i<phi.m;i++){
            T p=phi(i)(it.index);
            if(p<min(0)){min=min.Remove_Index(2).Insert(p,0);data.trust=data.trust.Remove_Index(1).Insert(i,0);}
            else if(p<min(1)){min=min.Remove_Index(2).Insert(p,1);data.trust(1)=i;}
            else if(p<min(2)) min(2)=p;
            if(p<=valid_width) data.valid_flags|=1<<i;}
        min-=(min(0)+min(1))/2; // In case of boundary conditions
        if(min(1)>valid_width || min(2)>extent){
            data=PAIRWISE_LEVEL_SET_DATA();
            continue;}
        if(min(1)>min(2)-trust_buffer){
            data.trust=VECTOR<short,2>(-1,-1);
            continue;}
        data.phi=min.x;
        if(data.trust.y<data.trust.x){
            exchange(data.trust.x,data.trust.y);
            data.phi=min.y;}}

    for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        PAIRWISE_LEVEL_SET_DATA& data=pairwise_data(it.index);
        if(data.trust.Contains(-1)) continue;
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,0));
    }
    Flush_Frame<TV>("trust");
}
//#####################################################################
// Function Initialize_Pairwise_Level_Set
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Initialize_Pairwise_Level_Set()
{
    pairwise_phi.Resize(phi.m);
    for(int i=0;i<phi.m;i++){
        pairwise_phi(i).Resize(phi.m);
        for(int j=i+1;j<phi.m;j++){
            pairwise_phi(i)(j).Resize(grid.Node_Indices(extrap_width),false);
            pairwise_phi(i)(j).Fill(default_phi);}}

    for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        PAIRWISE_LEVEL_SET_DATA& data=pairwise_data(it.index);
        if(data.trust.x>=0){
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(data.trust.y==1,data.trust.Sum()==2,data.trust.x==1)/(1+(data.phi<0)));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(data.phi));
            pairwise_phi(data.trust.x)(data.trust.y)(it.index)=data.phi;}}
    Flush_Frame<TV>(__FUNCTION__);
}
//#####################################################################
// Function Fill_Valid_Region_With_Exprapolation
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Fill_Valid_Region_With_Exprapolation()
{
    HASHTABLE<VECTOR<short,2> > trust_pairs;
    for(RANGE_ITERATOR<TV::m> it(pairwise_data.domain);it.Valid();it.Next()){
        VECTOR<short,2> trust=pairwise_data(it.index).trust;
        if(trust.x!=-1)
            trust_pairs.Set(trust);}
    for(HASHTABLE<VECTOR<short,2> >::ITERATOR it(trust_pairs);it.Valid();it.Next()){
        VECTOR<short,2> k=it.Key();
        EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::Extrapolate_Node(grid,
            [&](const TV_INT& index){return pairwise_data(index).trust==k;},
            extrap_width,pairwise_phi(k.x)(k.y),3,extrap_width);
        for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
            int mask=(1<<k.x)|(1<<k.y);
            T& p=pairwise_phi(k.x)(k.y)(it.index);
            if((pairwise_data(it.index).valid_flags&mask)==mask){
                PHYSBAM_ASSERT(p!=default_phi);}
            else p=default_phi;}
        Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(4),VECTOR<T,3>(1,0,0));
        Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(5),VECTOR<T,3>(0,1,0));
        Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(5),VECTOR<T,3>(0,0,1));
        Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(6),VECTOR<T,3>(0,1,1));
        Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(6),VECTOR<T,3>(1,0,1));
        Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(5)(6),VECTOR<T,3>(1,1,0));
        Flush_Frame<TV>("extrap");}
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(4),VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(5),VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(5),VECTOR<T,3>(0,0,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(6),VECTOR<T,3>(0,1,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(6),VECTOR<T,3>(1,0,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(5)(6),VECTOR<T,3>(1,1,0));
    Flush_Frame<TV>("after extrap");
}
//#####################################################################
// Function Meet_Phi
//#####################################################################
template<class T> static VECTOR<T,2>
Meet_Phi(const VECTOR<VECTOR<T,4> ,2>& phi)
{
    typedef VECTOR<T,2> TV;
    QUADRATIC<T> quad(-phi(0)(0)*phi(1)(3)+phi(0)(3)*phi(1)(0)+phi(0)(1)*phi(1)(2)-phi(0)(1)*phi(1)(0)+phi(0)(2)*phi(1)(3)-phi(0)(2)*phi(1)(1)+phi(0)(0)*phi(1)(1)-phi(0)(3)*phi(1)(2),-2*phi(0)(0)*phi(1)(1)+phi(0)(0)*phi(1)(3)-phi(0)(1)*phi(1)(2)+2*phi(0)(1)*phi(1)(0)+phi(0)(2)*phi(1)(1)-phi(0)(3)*phi(1)(0),-phi(0)(1)*phi(1)(0)+phi(0)(0)*phi(1)(1));
    quad.Compute_Roots();
    if(quad.roots==0){
        return TV()+1e20;
//        printf("XXX   %.15g %.15g %.15g\n", quad.a, quad.b, quad.c);
        quad.roots=1;
        quad.root1=-quad.b/(2*quad.a);}
    T ya[2]={quad.root1,quad.root2},ph[2]={FLT_MAX,FLT_MAX},ps[2]={FLT_MAX,FLT_MAX};(void)ps;
    TV X[2];
    for(int i=0;i<quad.roots;i++){
        if(abs(ya[i])>10) continue;
        T num=-phi(1)(0)-ya[i]*phi(1)(2)+ya[i]*phi(1)(0);
        T den=phi(1)(1)-phi(1)(0)+ya[i]*phi(1)(3)-ya[i]*phi(1)(2)-ya[i]*phi(1)(1)+ya[i]*phi(1)(0);
        X[i]=TV(num/den,ya[i]);
        ph[i]=phi(0)(0)+(phi(0)(1)-phi(0)(0))*X[i].x+(phi(0)(2)-phi(0)(0))*X[i].y+(phi(0)(3)-phi(0)(2)-phi(0)(1)+phi(0)(0))*X[i].x*X[i].y;
        ps[i]=phi(1)(0)+(phi(1)(1)-phi(1)(0))*X[i].x+(phi(1)(2)-phi(1)(0))*X[i].y+(phi(1)(3)-phi(1)(2)-phi(1)(1)+phi(1)(0))*X[i].x*X[i].y;
}

//    LOG::cout<<phi<<std::endl;
//    printf("xx %.15g %.15g    %.15g %.15g\n", ph[0],ps[0],ph[1],ps[1]);
    return X[abs(ph[1])<abs(ph[0])];
}
//#####################################################################
// Function Meet_Phi
//#####################################################################
template<class T> static VECTOR<T,3>
Meet_Phi(const VECTOR<VECTOR<T,8> ,2>& phi)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Zero_Phi
//#####################################################################
template<class T> VECTOR<T,2>
Zero_Phi(const VECTOR<VECTOR<T,4> ,3>& phi,VECTOR<T,3>& p)
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<T,(1<<TV::m)> PHI;
    p=VECTOR<T,3>()+FLT_MAX;
    TV bestX=TV()+1e20,X;
    T pp=0;
    X=Meet_Phi(VECTOR<PHI,2>(phi.x-phi.y,phi.x-phi.z));
    if(X.Magnitude_Squared()<1e10){
        pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
        if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(pp,pp,pp);bestX=X;}}
    X=Meet_Phi(VECTOR<PHI,2>(phi.x+phi.y,phi.x-phi.z));
    if(X.Magnitude_Squared()<1e10){
        pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
        if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(pp,-pp,pp);bestX=X;}}
    X=Meet_Phi(VECTOR<PHI,2>(phi.x-phi.y,phi.x+phi.z));
    if(X.Magnitude_Squared()<1e10){
        pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
        if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(pp,pp,-pp);bestX=X;}}
    X=Meet_Phi(VECTOR<PHI,2>(-phi.x+phi.y,-phi.x+phi.z));
    if(X.Magnitude_Squared()<1e10){
        pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
        if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(-pp,pp,pp);bestX=X;}}
    return bestX;
}
//#####################################################################
// Function Zero_Phi
//#####################################################################
template<class T> VECTOR<T,3>
Zero_Phi(const VECTOR<VECTOR<T,8> ,3>& phi,VECTOR<T,3>& p)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function One_Step_Triple_Junction_Correction
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
One_Step_Triple_Junction_Correction()
{
    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(4),VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(5),VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(5),VECTOR<T,3>(0,0,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(6),VECTOR<T,3>(0,1,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(6),VECTOR<T,3>(1,0,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(5)(6),VECTOR<T,3>(1,1,0));

    T max_move=(T).5*grid.dX.Max();
    HASHTABLE<TRIPLE<int,int,TV_INT>,T> total_size;
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        int mask=~0;
        for(int j=0;j<(1<<TV::m);j++)
            mask&=pairwise_data(bits(j)+it.index).valid_flags;
        for(int a=0;a<phi.m;a++)
            if(mask&(1<<a))
                for(int b=a+1;b<phi.m;b++)
                    if(mask&(1<<b))
                        for(int c=b+1;c<phi.m;c++)
                            if(mask&(1<<c)){
                                PHI ab(pairwise_phi(a)(b).Subset(bits+it.index));
                                PHI ac(pairwise_phi(a)(c).Subset(bits+it.index));
                                PHI bc(pairwise_phi(b)(c).Subset(bits+it.index));

                                bool b01=RANGE<TV>::Unit_Box().Lazy_Inside(Meet_Phi(VECTOR<PHI,2>(ab,ac)));
                                bool b02=RANGE<TV>::Unit_Box().Lazy_Inside(Meet_Phi(VECTOR<PHI,2>(ab,bc)));
                                bool b12=RANGE<TV>::Unit_Box().Lazy_Inside(Meet_Phi(VECTOR<PHI,2>(ac,bc)));
                                if(!b01 && !b02 && !b12) continue;
                                // if(b01){Add_Debug_Particle(Meet_Phi(VECTOR<PHI,2>(ab,ac))*grid.dX+grid.Node(it.index),VECTOR<T,3>(1,0,0));LOG::cout<<std::setprecision(16)<<it.index<<"   "<<ab<<"   "<<ac<<"   "<<Meet_Phi(VECTOR<PHI,2>(ab,ac))<<std::endl;}
                                // if(b02) Add_Debug_Particle(Meet_Phi(VECTOR<PHI,2>(ab,bc))*grid.dX+grid.Node(it.index),VECTOR<T,3>(0,1,0));
                                // if(b12) Add_Debug_Particle(Meet_Phi(VECTOR<PHI,2>(ac,bc))*grid.dX+grid.Node(it.index),VECTOR<T,3>(0,0,1));
//                                LOG::cout<<it.index<<"  "<<b01<<b02<<b12<<"   "<<ab<<"   "<<ac<<"   "<<bc<<std::endl;

                                VECTOR<T,3> pp;
                                TV X=Zero_Phi(VECTOR<PHI,3>(ab,ac,bc),pp),Y=X*grid.dX+grid.Node(it.index);
                                Add_Debug_Particle(Y,VECTOR<T,3>(pp(0)>0,pp(1)>0,pp(2)>0));
                                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(pp(0)));
//                                Add_Debug_Particle(it.Location(),VECTOR<T,3>(pp(0)>0,pp(1)>0,pp(2)>0));
                                if(abs(pp(0))>max_move) pp*=max_move/abs(pp(0));
                                T div=(b01 && b02 && b12)?1:3;
                                for(int j=0;j<(1<<TV::m);j++){
                                    T& p01=total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(a,b,bits(j)+it.index));
                                    T& p02=total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(a,c,bits(j)+it.index));
                                    T& p12=total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(b,c,bits(j)+it.index));
                                    if(abs(p01)<abs(pp(0))) p01=pp(0)/div;
                                    if(abs(p02)<abs(pp(1))) p02=pp(1)/div;
                                    if(abs(p12)<abs(pp(2))) p12=pp(2)/div;}}}
    Flush_Frame<TV>("start");

    for(typename HASHTABLE<TRIPLE<int,int,TV_INT>,T>::ITERATOR it(total_size);it.Valid();it.Next())
        pairwise_phi(it.Key().x)(it.Key().y)(it.Key().z)-=it.Data();
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(4),VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(5),VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(5),VECTOR<T,3>(0,0,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(6),VECTOR<T,3>(0,1,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(6),VECTOR<T,3>(1,0,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(5)(6),VECTOR<T,3>(1,1,0));
    Flush_Frame<TV>(__FUNCTION__);
}
//#####################################################################
// Function Fill_Combined_Level_Set_At_Index
//#####################################################################
template<class TV> int TRIPLE_JUNCTION_CORRECTION<TV>::
Fill_Combined_Level_Set_At_Index(const TV_INT& node)
{
    int min_index=-1;
    T min1=FLT_MAX,min2=FLT_MAX;
    for(int i=0;i<phi.m;i++){
        T p=phi(i)(node);
        if(p<min1){min2=min1;min1=p;min_index=i;}
        else if(p<min2) min2=p;}
    T shift=(T).5*(min2+min1);
    for(int i=0;i<phi.m;i++)
        phi(i)(node)-=shift;
    combined_color(node)=min_index;
    combined_phi(node)=min1-shift;
    return min_index;
}
//#####################################################################
// Function Fill_Combined_Level_Set_At_Index
//#####################################################################
template<class TV> int TRIPLE_JUNCTION_CORRECTION<TV>::
Fill_Phi_From_Pairwise_Level_Set_At_Index(const TV_INT& node,int color)
{
    T min_phi=FLT_MAX;
    int mask=pairwise_data(node).valid_flags;
    for(int a=0;a<phi.m;a++){
        if(a==color || !(mask&(1<<a))) continue;
        T p=a<color?pairwise_phi(a)(color)(node):-pairwise_phi(color)(a)(node);
        min_phi=min(p,min_phi);
        phi(a)(node)=p;}
    phi(color)(node)=-min_phi;
    return Fill_Combined_Level_Set_At_Index(node);
}
//#####################################################################
// Function Update_Color_Level_Sets
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Update_Color_Level_Sets()
{
    combined_phi.Resize(grid.Node_Indices(ghost));
    combined_color.Resize(grid.Node_Indices(ghost));

    ARRAY<TV_INT> todo;
    for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        int mask=pairwise_data(it.index).valid_flags;
        int num_bits=count_bits(mask);
        if(num_bits<2 || !pairwise_data(it.index).trust.Contains(-1)){
            todo.Append(it.index);
            int color=Fill_Combined_Level_Set_At_Index(it.index);
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(color==0,color==1,color==2));
            continue;}
        combined_color(it.index)=-1;}

    while(todo.m){
        TV_INT index(todo.Pop());
        int mask=pairwise_data(index).valid_flags;
        int color=combined_color(index);

        for(int i=0;i<GRID<TV>::number_of_one_ring_neighbors_per_cell;i++){
            TV_INT neighbor=GRID<TV>::One_Ring_Neighbor(index,i);
            if(!combined_color.domain.Lazy_Inside_Half_Open(neighbor)) continue;
            if(combined_color(neighbor)>=0) continue;
            int neighbor_mask=pairwise_data(neighbor).valid_flags&mask;
            bool separated=false;
            for(int a=0;a<phi.m;a++)
                if(neighbor_mask&(1<<a))
                    for(int b=a+1;b<phi.m;b++)
                        if(neighbor_mask&(1<<b))
                            if((pairwise_phi(a)(b)(index)>0)!=(pairwise_phi(a)(b)(neighbor)>0)){
                                separated=true;
                                break;}
            if(separated) continue;
            Fill_Phi_From_Pairwise_Level_Set_At_Index(neighbor,color);
            combined_color(neighbor)=color;
            todo.Append(neighbor);
            Add_Debug_Object(VECTOR<TV,2>(grid.Node(index),grid.Node(neighbor)),VECTOR<T,3>(color==0,color==1,color==2)/3);
            Add_Debug_Particle(grid.Node(neighbor),VECTOR<T,3>(color==0,color==1,color==2)/3);}}

    for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next())
        PHYSBAM_ASSERT(combined_color(it.index)>=0);
    Flush_Frame<TV>(__FUNCTION__);
}
//#####################################################################
// Function Cut_Interface
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Cut_Interface(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_data)
{
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(4),VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(5),VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(5),VECTOR<T,3>(0,0,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(3)(6),VECTOR<T,3>(0,1,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(4)(6),VECTOR<T,3>(1,0,1));
    Dump_Interface<T,TV,TV_INT>(grid,pairwise_phi(5)(6),VECTOR<T,3>(1,1,0));
    Flush_Frame<TV>("pairwise ls");

    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();
    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        int full_mask=~0;
        for(int j=0;j<bits.m;j++){
            full_mask&=pairwise_data(bits(j)+it.index).valid_flags;}
        VECTOR<int,(1<<TV::m)> colors(combined_color.Subset(bits+it.index));
        if(count_bits(full_mask)<3){
            if(colors.Count_Matches(colors(0))<(1<<TV::m)){
                VECTOR<T,(1<<TV::m)> phis(combined_phi.Subset(bits+it.index));
                CELL_ELEMENTS& ce=index_to_cell_data.Get_Or_Insert(it.index);
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,1));
                MARCHING_CUBES_COLOR<TV>::Get_Elements_For_Cell(ce.interface,ce.boundary,colors-bc_colors,phis,
                    grid.Cell_Domain(it.index));}
            continue;}

        int lone_a=-1,lone_b=-1;
        for(int a=0;a<phi.m;a++)
            if(full_mask&(1<<a))
                for(int b=a+1;b<phi.m;b++)
                    if(full_mask&(1<<b)){
                        if(!colors.Contains(a) && !colors.Contains(b)) continue;
                        VECTOR<T,(1<<TV::m)> phis(pairwise_phi(a)(b).Subset(bits+it.index));
                        PHYSBAM_ASSERT(phis.Max() && phis.Min());
                        if(phis.Max()>0 && phis.Min()<0){
                            if(lone_a==-1){
                                lone_a=a;
                                lone_b=b;}
                            else{
                                lone_a=-2;
                                a=phi.m;
                                break;}}}
        if(lone_a==-1) continue;
        if(lone_a>=0 && colors.Contains(lone_a) && colors.Contains(lone_b)){
            VECTOR<T,(1<<TV::m)> phis(pairwise_phi(lone_a)(lone_b).Subset(bits+it.index));
            CELL_ELEMENTS& ce=index_to_cell_data.Get_Or_Insert(it.index);
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,0,1));
//            LOG::cout<<phis<<"  "<<(colors-bc_colors)<<std::endl;
            for(int i=0;i<phis.m;i++) phis(i)=abs(phis(i));
            MARCHING_CUBES_COLOR<TV>::Get_Elements_For_Cell(ce.interface,ce.boundary,colors-bc_colors,phis,
                grid.Cell_Domain(it.index));
            continue;}

        Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
        Cut_Cell_With_Pairwise_Phi(index_to_cell_data,it.index);}
    Flush_Frame<TV>(__FUNCTION__);
    combined_color.array-=bc_colors;
}
//#####################################################################
// Function Cut_Stencil_With_Pairwise_Phi
//#####################################################################
template<class TV,class TV_INT,class CELL_ELEMENTS> static void
Cut_Cell_With_Pairwise_Phi_Helper(TRIPLE_JUNCTION_CORRECTION<TV>& self,HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_data,const VECTOR<int,2>& cell)
{
    typedef typename TV::SCALAR T;
    typedef typename MARCHING_CUBES_COLOR<TV>::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    typedef typename MARCHING_CUBES_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    typedef VECTOR<T,(1<<TV::m)> PHI;
    struct CROSSING
    {
        T theta;
        int c0;
        int c1;
        TV X;
    };

    int vertices_of_side[4][2]={{0,1},{1,3},{3,2},{2,0}};

    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());
    int full_mask=~0;
    for(int j=0;j<bits.m;j++)
        full_mask&=self.pairwise_data(bits(j)+cell).valid_flags;

    VECTOR<int,(1<<TV::m)> colors(self.combined_color.Subset(bits+cell));
    VECTOR<ARRAY<CROSSING>,4> crossings;
    for(int a=0;a<self.phi.m;a++)
        if(full_mask&(1<<a))
            for(int b=a+1;b<self.phi.m;b++)
                if(full_mask&(1<<b)){
                    VECTOR<T,(1<<TV::m)> phis(self.pairwise_phi(a)(b).Subset(bits+cell));
                    if(phis.Max()<=0 || phis.Min()>=0) continue;
                    for(int e=0;e<4;e++){
                        int v0=vertices_of_side[e][0],v1=vertices_of_side[e][1];
                        T p0=phis(v0),p1=phis(v1);
                        if((p0>0)==(p1>0)) continue;
                        CROSSING t={p0/(p0-p1),a,b};
                        TV X0=self.grid.Node(cell+bits(v0)),X1=self.grid.Node(cell+bits(v1));
                        t.X=X0+t.theta*(X1-X0);
                        if(p0>0) exchange(t.c0,t.c1);
                        crossings(e).Append(t);}}

    ARRAY<int> best(self.phi.m),index(self.phi.m);
    for(int e=0;e<4;e++){
        ARRAY<CROSSING>& ar=crossings(e);
        ar.Sort([](const CROSSING& a,const CROSSING& b){return a.theta<b.theta;});
        int c0=colors(vertices_of_side[e][0]),c1=colors(vertices_of_side[e][1]),k=0;
        ARRAY<int> pa(ar.m);
        best.Fill(-INT_MAX);
        best(c1)=0;
        index.Fill(-1);
        for(int i=ar.m-1;i>=0;i--){
            pa(i)=index(ar(i).c1);
            if(best(ar(i).c1)+1>best(ar(i).c0)){
                best(ar(i).c0)=best(ar(i).c1)+1;
                index(ar(i).c0)=i;}}
        for(int i=index(c0);i>=0;i=pa(i))
            ar(k++)=ar(i);
        ar.Resize(k);}
    // if(cell.x==50 && cell.y==23){
    //     LOG::cout<<colors<<std::endl;
    // for(int e=0;e<4;e++){
    //     LOG::cout<<e<<" : ";
    //     ARRAY<CROSSING>& ar=crossings(e);
    //     for(int i=0;i<ar.m;i++){
    //         LOG::cout<<ar(i).c0<<"  "<<ar(i).c1<<"     ";
    //     }
    //     LOG::cout<<std::endl;
    // }}

    bool empty=true;
    for(int e=0;e<4;e++)
        if(crossings(e).m)
            empty=false;
    if(empty) return;

    CELL_ELEMENTS& ce=index_to_cell_data.Get_Or_Insert(cell);
    for(int s=0;s<2;s++){
        int c0=colors(3*s);
        TV X0=self.grid.Node(cell+bits(3*s)),X1=self.grid.Node(cell+bits(3*s^1)),X=X0;
        for(int i=0;i<crossings(2*s).m;i++){
            BOUNDARY_ELEMENT be={SEGMENT_2D<T>(X,crossings(2*s)(i).X),c0-self.bc_colors};
            ce.boundary.Append(be);
            X=crossings(2*s)(i).X;
            c0=crossings(2*s)(i).c1;}
        BOUNDARY_ELEMENT be={SEGMENT_2D<T>(X,X1),c0-self.bc_colors};
        ce.boundary.Append(be);}

    HASHTABLE<VECTOR<int,2>,TV> pts;
    for(int e=0;e<4;e++)
        for(int i=0;i<crossings(e).m;i++){
            VECTOR<int,2> key0(crossings(e)(i).c0,crossings(e)(i).c1),key1(crossings(e)(i).c1,crossings(e)(i).c0);
            TV X0=crossings(e)(i).X,X1;
            if(pts.Get(key1,X1)){
                pts.Delete(key1);
                if(key0.x>key0.y){
                    exchange(X0,X1);
                    key1=key0;}
                INTERFACE_ELEMENT ie={SEGMENT_2D<T>(X0,X1),key1-self.bc_colors};
                ce.interface.Append(ie);}
            pts.Insert(key0,X0);}

    TV centroid;
    RANGE<TV> cell_range(self.grid.Cell_Domain(cell));
    HASHTABLE<VECTOR<int,3> > found;
    for(typename HASHTABLE<VECTOR<int,2>,TV>::ITERATOR it(pts);it.Valid();it.Next()){
        for(int c=0;c<self.phi.m;c++)
            if(c!=it.Key().x && c!=it.Key().y && full_mask&(1<<c)){
                VECTOR<int,2> key(it.Key().y,c);
                if(!pts.Contains(key)) continue;
                VECTOR<int,3> tkey(key.Append(it.Key().x).Sorted());
                if(found.Contains(tkey)) continue;
                found.Insert(tkey);
                VECTOR<T,3> p;
                PHI phi0(self.pairwise_phi(tkey.x)(tkey.y).Subset(bits+cell));
                PHI phi1(self.pairwise_phi(tkey.x)(tkey.z).Subset(bits+cell));
                PHI phi2(self.pairwise_phi(tkey.y)(tkey.z).Subset(bits+cell));
                centroid+=Zero_Phi(VECTOR<PHI,3>(phi0,phi1,phi2),p)*cell_range.Edge_Lengths()+cell_range.min_corner;}}

    if(found.Size()){
        centroid/=found.Size();
        for(typename HASHTABLE<VECTOR<int,2>,TV>::ITERATOR it(pts);it.Valid();it.Next()){
            INTERFACE_ELEMENT ie={SEGMENT_2D<T>(centroid,it.Data()),it.Key()-self.bc_colors};
            if(ie.color_pair.x>ie.color_pair.y){
                exchange(ie.color_pair.x,ie.color_pair.y);
                exchange(ie.face.X.x,ie.face.X.y);}
            ce.interface.Append(ie);}}
}
//#####################################################################
// Function Cut_Stencil_With_Pairwise_Phi
//#####################################################################
template<class TV,class TV_INT,class CELL_ELEMENTS> static void
Cut_Cell_With_Pairwise_Phi_Helper(TRIPLE_JUNCTION_CORRECTION<TV>& self,HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_data,const VECTOR<int,3>& cell)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Cut_Cell_With_Pairwise_Phi
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Cut_Cell_With_Pairwise_Phi(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_data,const TV_INT& cell)
{
    Cut_Cell_With_Pairwise_Phi_Helper(*this,index_to_cell_data,cell);
}
//#####################################################################
// Function Compute_Pairwise_Level_Set_Data
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Compute_Pairwise_Level_Set_Data()
{
    Compute_Pairwise_Data();
    Initialize_Pairwise_Level_Set();

    Fill_Valid_Region_With_Exprapolation();

    for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        T p=pairwise_phi(0)(1)(it.index);
        if((pairwise_data(it.index).valid_flags&3)==3)
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(pairwise_data(it.index).trust==VECTOR<short,2>(0,1),0,1));
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(p<0,p>=0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}

    for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        T p=pairwise_phi(0)(1)(it.index);
        if((pairwise_data(it.index).valid_flags&3)==3)
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(pairwise_data(it.index).trust==VECTOR<short,2>(0,1),0,1));
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(p<0,p>=0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
    Flush_Frame<TV>(__FUNCTION__);

    for(int t=0;t<20;t++)
        One_Step_Triple_Junction_Correction();

    Update_Color_Level_Sets();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIPLE_JUNCTION_CORRECTION<TV>::PAIRWISE_LEVEL_SET_DATA::
PAIRWISE_LEVEL_SET_DATA()
    :phi(FLT_MAX),valid_flags(0),trust(-1,-1)
{
}
namespace PhysBAM{
template class TRIPLE_JUNCTION_CORRECTION<VECTOR<float,2> >;
template class TRIPLE_JUNCTION_CORRECTION<VECTOR<double,2> >;
template class TRIPLE_JUNCTION_CORRECTION<VECTOR<float,3> >;
template class TRIPLE_JUNCTION_CORRECTION<VECTOR<double,3> >;
}
