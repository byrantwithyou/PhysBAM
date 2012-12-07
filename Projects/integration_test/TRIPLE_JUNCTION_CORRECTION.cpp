#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_VIEW.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_TETRAHEDRA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "TRIPLE_JUNCTION_CORRECTION.h"
using namespace PhysBAM;
template<class T,class TV>
extern void Flush_Frame(const char* title);
template<class TV> GRID<TV>* Global_Grid(GRID<TV>* grid_in=0);

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

template<class T,class TV_INT>
void Dump_Interface(ARRAY_VIEW<T,TV_INT> p,const VECTOR<T,3>& col)
{
    typedef VECTOR<T,TV_INT::m> TV;
    const GRID<TV>& grid=*Global_Grid<TV>();
    T default_phi=(T)1.12349871352389e-10;
    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        VECTOR<T,4> phi(p.Subset(bits+it.index));
        if(phi.Contains(default_phi)) continue;
        if(phi.Min()>=0 || phi.Max()<=0) continue;
        Handle_Square(grid,phi,it.index,col);}
}


//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIPLE_JUNCTION_CORRECTION<TV>::
TRIPLE_JUNCTION_CORRECTION(const GRID<TV>& grid,ARRAY<ARRAY<T,TV_INT> >& phi,int ghost)
    :grid(grid),phi(phi),ghost(ghost),default_phi((T)1.12349871352389e-10),trust_buffer(grid.dX.Max()),
    valid_width(ghost*grid.dX.Max()),extent(3*valid_width),extrap_width(2*ghost+3)
{
}
//#####################################################################
// Function Compute_Pairwise_Data
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Compute_Pairwise_Data()
{
    pairwise_data.Resize(grid.Node_Indices(ghost));

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
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

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
        PAIRWISE_LEVEL_SET_DATA& data=pairwise_data(it.index);
        if(data.trust.x>=0){
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(data.trust.y==1,data.trust.Sum()==2,data.trust.x==1)/(1+(data.phi<0)));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(data.phi));
            pairwise_phi(data.trust.x)(data.trust.y)(it.index)=data.phi;}}
}
//#####################################################################
// Function Fill_Valid_Region_With_Exprapolation
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Fill_Valid_Region_With_Exprapolation()
{
    for(int i=0;i<phi.m;i++)
        for(int j=i+1;j<phi.m;j++){
            struct LOCAL_MASK:public EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::MASK
            {
                const ARRAY<PAIRWISE_LEVEL_SET_DATA,TV_INT>& pairwise_data;
                VECTOR<short,2> current_pair;
                LOCAL_MASK(const ARRAY<PAIRWISE_LEVEL_SET_DATA,TV_INT>& pairwise_data,const VECTOR<short,2>& current_pair)
                    :pairwise_data(pairwise_data),current_pair(current_pair)
                {}
                bool Inside(const TV_INT& index) PHYSBAM_OVERRIDE
                {return pairwise_data(index).trust==current_pair;}
            } inside_mask(pairwise_data,VECTOR<short,2>(i,j));
            EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T>::Extrapolate_Node(grid,inside_mask,extrap_width,pairwise_phi(i)(j),3,extrap_width);
            for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
                int mask=(1<<i)|(1<<j);
                T& p=pairwise_phi(i)(j)(it.index);
                if((pairwise_data(it.index).valid_flags&mask)==mask){
                    PHYSBAM_ASSERT(p!=default_phi);}
                else p=default_phi;}}
}
//#####################################################################
// Function One_Step_Triple_Junction_Correction
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
One_Step_Triple_Junction_Correction()
{
    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());

    T max_move=(T).5*grid.dX.Max();
    HASHTABLE<TRIPLE<int,int,TV_INT>,T> total_size;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
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

                                VECTOR<T,3> pp;
                                TV X=Zero_Phi(VECTOR<PHI,3>(ab,ac,bc),pp),Y=grid.domain.min_corner+grid.dX*(X+TV(it.index));
                                Add_Debug_Particle(Y,VECTOR<T,3>(pp(0)>0,pp(1)>0,pp(2)>0));
                                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(pp(0)));
                                Add_Debug_Particle(it.Location(),VECTOR<T,3>(pp(0)>0,pp(1)>0,pp(2)>0));
                                if(abs(pp(0))>max_move) pp*=max_move/abs(pp(0));
                                for(int j=0;j<(1<<TV::m);j++){
                                    T& p01=total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(a,b,bits(j)+it.index));
                                    T& p02=total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(a,c,bits(j)+it.index));
                                    T& p12=total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(b,c,bits(j)+it.index));
                                    if(abs(p01)<abs(pp(0))) p01=pp(0);
                                    if(abs(p02)<abs(pp(1))) p02=pp(1);
                                    if(abs(p12)<abs(pp(2))) p12=pp(2);}}}

    Dump_Interface<T,TV_INT>(pairwise_phi(0)(1),VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(0)(2),VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(1)(2),VECTOR<T,3>(0,0,1));
    Flush_Frame<T,TV>("meet points");
    for(typename HASHTABLE<TRIPLE<int,int,TV_INT>,T>::ITERATOR it(total_size);it.Valid();it.Next())
        pairwise_phi(it.Key().x)(it.Key().y)(it.Key().z)-=it.Data();
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
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
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

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next())
        PHYSBAM_ASSERT(combined_color(it.index)>=0);

    Dump_Interface<T,TV_INT>(pairwise_phi(0)(1),VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(0)(2),VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(1)(2),VECTOR<T,3>(0,0,1));
    Flush_Frame<T,TV>("flood fill");
}
//#####################################################################
// Function Cut_Interface
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Cut_Interface(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_data)
{
#if 0
    const VECTOR<TV_INT,(1<<TV::m)>& bits=GRID<TV>::Binary_Counts(TV_INT());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        CELL_ELEMENTS ce;
        const VECTOR<VECTOR<TV_INT,TV::m+1>,num_stencils>& stencils_parity=stencils(it.index.Sum()&1);
        int full_mask=0;
        for(int j=0;j<bits.m;j++)
            full_mask|=pairwise_data(bits(j)).valid_flags;
        if(count_bits(full_mask)<3){
            VECTOR<int,(1<<TV::m)> colors(combined_color.Subset(bits+it.index));
            if(colors.Count_Matches(colors(0))<=(1<<TV::m))
                for(int i=0;i<stencils_parity.m;i++)
                    Cut_Stencil_With_Phi(ce,it.index,VECTOR<TV_INT,TV::m+1>(stencils_parity(i)+it.index));
            continue;}

        for(int i=0;i<stencils_parity.m;i++){
            VECTOR<TV_INT,TV::m+1> st(stencils_parity(i)+it.index);
            int mask=~0;
            for(int j=0;j<st.m;j++)
                mask&=pairwise_data(st(j)).valid_flags;
            if(count_bits(mask)<3) Cut_Stencil_With_Phi(ce,it.index,st);
            else Cut_Stencil_With_Pairwise_Phi(ce,it.index,st);}
        if(ce.interface.m)
            index_to_cell_data.Get_Or_Insert(it.index)=ce;}
#endif
}
//#####################################################################
// Function Cut_Stencil_With_Phi
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Cut_Stencil_With_Phi(CELL_ELEMENTS& ce,const TV_INT& cell,const VECTOR<TV_INT,TV::m+1>& st)
{
    VECTOR<TV,TV::m+1> X;
    VECTOR<T,TV::m+1> element_phi;
    VECTOR<int,2> color_pair;
    for(int i=0;i<TV::m+1;i++){
        X(i)=grid.Node(st(i));
        T mn=combined_phi(st(i));
        int col=combined_color(st(i));
        if(color_pair.x==-1 || color_pair.x==col) color_pair.x=col;
        else{
            PHYSBAM_ASSERT(color_pair.y==-1 || color_pair.y==col);
            if(color_pair.x>col){
                element_phi=-element_phi;
                color_pair.y=color_pair.x;
                color_pair.x=col;}
            else{
                mn=-mn;
                color_pair.y=col;}}
        element_phi(i)=mn;}

    ARRAY<T_FACE> interface,boundary_faces[2];
    VECTOR<VECTOR<ARRAY<T_FACE>*,2>,TV::m+1> boundary;
    int spec=0,num_diff=0;
    for(int i=1;i<TV::m+1;i++)
        if(st(i).x!=st(0).x){
            num_diff++;
            spec=i;}
    if(num_diff==TV::m) spec=0;
    if(TV::m!=3 || num_diff!=2)
        boundary(spec)=VECTOR<ARRAY<T_FACE>*,2>(&boundary_faces[1],&boundary_faces[0]);
    MARCHING_TETRAHEDRA<TV>::Get_Elements_For_Tetrahedron(interface,boundary,element_phi,X);

    for(int i=0;i<interface.m;i++){
        typename MARCHING_CUBES_COLOR<TV>::INTERFACE_ELEMENT e={interface(i),color_pair};
        ce.interface.Append(e);}

    for(int s=0;s<2;s++)
        for(int i=0;i<boundary_faces[s].m;i++){
            typename MARCHING_CUBES_COLOR<TV>::BOUNDARY_ELEMENT e={boundary_faces[s](i),color_pair(s)};
            ce.boundary.Append(e);}
}
//#####################################################################
// Function Cut_Stencil_With_Pairwise_Phi
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Cut_Stencil_With_Pairwise_Phi(CELL_ELEMENTS& ce,const TV_INT& cell,const VECTOR<TV_INT,TV::m+1>& st)
{
    
}
//#####################################################################
// Function Fill_Combined_Level_Set
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Fill_Combined_Level_Set()
{
}
//#####################################################################
// Function Meet_Phi
//#####################################################################
template<class TV> TV TRIPLE_JUNCTION_CORRECTION<TV>::
Meet_Phi(const VECTOR<PHI,2>& phi)
{
    QUADRATIC<T> quad(-phi(0)(0)*phi(1)(3)+phi(0)(3)*phi(1)(0)+phi(0)(1)*phi(1)(2)-phi(0)(1)*phi(1)(0)+phi(0)(2)*phi(1)(3)-phi(0)(2)*phi(1)(1)+phi(0)(0)*phi(1)(1)-phi(0)(3)*phi(1)(2),-2*phi(0)(0)*phi(1)(1)+phi(0)(0)*phi(1)(3)-phi(0)(1)*phi(1)(2)+2*phi(0)(1)*phi(1)(0)+phi(0)(2)*phi(1)(1)-phi(0)(3)*phi(1)(0),-phi(0)(1)*phi(1)(0)+phi(0)(0)*phi(1)(1));
    quad.Compute_Roots();
    if(quad.roots==0){
        quad.roots=1;
        quad.root1=-quad.b/(2*quad.a);}
    T ya[2]={quad.root1,quad.root2},ph[2]={FLT_MAX,FLT_MAX};
    TV X[2];
    for(int i=0;i<quad.roots;i++){
        if(abs(ya[i])>10) continue;
        T num=-phi(1)(0)-ya[i]*phi(1)(2)+ya[i]*phi(1)(0);
        T den=phi(1)(1)-phi(1)(0)+ya[i]*phi(1)(3)-ya[i]*phi(1)(2)-ya[i]*phi(1)(1)+ya[i]*phi(1)(0);
        X[i]=TV(num/den,ya[i]);
        ph[i]=phi(0)(0)+(phi(0)(1)-phi(0)(0))*X[i].x+(phi(0)(2)-phi(0)(0))*X[i].y+(phi(0)(3)-phi(0)(2)-phi(0)(1)+phi(0)(0))*X[i].x*X[i].y;}

    return X[abs(ph[1])<abs(ph[0])];
}
//#####################################################################
// Function Zero_Phi
//#####################################################################
template<class TV> TV TRIPLE_JUNCTION_CORRECTION<TV>::
Zero_Phi(const VECTOR<PHI,3>& phi,VECTOR<T,3>& p)
{
    p=VECTOR<T,3>()+FLT_MAX;
    TV bestX,X;
    T pp=0;
    X=Meet_Phi(VECTOR<PHI,2>(phi.x-phi.y,phi.x-phi.z));
    pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
    if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(pp,pp,pp);bestX=X;}
    X=Meet_Phi(VECTOR<PHI,2>(phi.x+phi.y,phi.x-phi.z));
    pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
    if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(pp,-pp,pp);bestX=X;}
    X=Meet_Phi(VECTOR<PHI,2>(phi.x-phi.y,phi.x+phi.z));
    pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
    if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(pp,pp,-pp);bestX=X;}
    X=Meet_Phi(VECTOR<PHI,2>(-phi.x+phi.y,-phi.x+phi.z));
    pp=LINEAR_INTERPOLATION<T,T>::Linear(&phi.x(0),X);
    if(abs(pp)<abs(p.x)){p=VECTOR<T,3>(-pp,pp,pp);bestX=X;}
    return bestX;
}
//#####################################################################
// Function Compute_Pairwise_Level_Set_Data
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Compute_Pairwise_Level_Set_Data()
{
    Compute_Pairwise_Data();
    Initialize_Pairwise_Level_Set();

    Dump_Interface<T,TV_INT>(pairwise_phi(0)(1),VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(0)(2),VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(1)(2),VECTOR<T,3>(0,0,1));
    Flush_Frame<T,TV>("Initial pairwise level sets");

    Fill_Valid_Region_With_Exprapolation();

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
        T p=pairwise_phi(0)(1)(it.index);
        if((pairwise_data(it.index).valid_flags&3)==3)
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(pairwise_data(it.index).trust==VECTOR<short,2>(0,1),0,1));
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(p<0,p>=0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
    Flush_Frame<T,TV>("level set 01");

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
        T p=pairwise_phi(0)(1)(it.index);
        if((pairwise_data(it.index).valid_flags&3)==3)
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(pairwise_data(it.index).trust==VECTOR<short,2>(0,1),0,1));
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(p<0,p>=0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(p));}
    Flush_Frame<T,TV>("level set 01");

    Dump_Interface<T,TV_INT>(pairwise_phi(0)(1),VECTOR<T,3>(1,0,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(0)(2),VECTOR<T,3>(0,1,0));
    Dump_Interface<T,TV_INT>(pairwise_phi(1)(2),VECTOR<T,3>(0,0,1));
    Flush_Frame<T,TV>("Extrapolated pairwise level sets");

    for(int t=0;t<3;t++){
        One_Step_Triple_Junction_Correction();
        Dump_Interface<T,TV_INT>(pairwise_phi(0)(1),VECTOR<T,3>(1,0,0));
        Dump_Interface<T,TV_INT>(pairwise_phi(0)(2),VECTOR<T,3>(0,1,0));
        Dump_Interface<T,TV_INT>(pairwise_phi(1)(2),VECTOR<T,3>(0,0,1));
        Flush_Frame<T,TV>("After triple junction correction");}

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
}
