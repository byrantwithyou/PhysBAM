//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
EXTRAPOLATION_HIGHER_ORDER()
//    :grid(grid_input),ghost(ghost_input),phi(phi_input),dt((T).5,(T).33333,(T).25)
{
//    dt*=grid.dX.Min();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
~EXTRAPOLATION_HIGHER_ORDER()
{
}
//#####################################################################
// Function Add_Neighbors
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Add_Neighbors(MAPPING& m,ARRAY<TV_INT>& next,const ARRAY<TV_INT>& neighbors,const TV_INT& index,int unregistered,int registered)
{
    for(int i=0;i<neighbors.m;i++){
        TV_INT neighbor_index=index+neighbors(i);
        int& n=m.node_to_index(neighbor_index);
        if(n!=unregistered) continue;
        next.Append(neighbor_index);
        n=registered;}
}
//#####################################################################
// Function Register_Nodes
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Register_Nodes(const GRID<TV>& grid,const LEVELSET<TV>& phi,boost::function<bool(const TV_INT& index)> inside_mask,int ghost,MAPPING& m,ARRAY<TV>& normal,
    ARRAY<VECTOR<STENCIL,TV::m> >& stencil,int order,int fill_width)
{
    m.node_to_index.Resize(grid.Domain_Indices(ghost+1)); // Need an extra ring for the sentinals
    m.index_to_node.Append(TV_INT()+INT_MAX); // First index is the "outside" index.
    normal.Append(TV()+FLT_MAX);

    ARRAY<TV_INT> next_inside,next_outside,current;
    ARRAY<TV_INT> neighbors_outside,neighbors_inside;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT()-1,TV_INT()+2));it.Valid();it.Next())
        neighbors_outside.Append(it.index);
    for(int d=0;d<TV::m;d++)
        for(int s=-1;s<=1;s+=2){
            TV_INT index;
            index(d)=s;
            neighbors_inside.Append(index);}

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        m.node_to_index(it.index)=-1;

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        bool b=inside_mask(it.index);
        m.node_to_index(it.index)=-1-b;}

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,1,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        Add_Neighbors(m,next_inside,neighbors_inside,it.index,-2,-4);

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        bool b=inside_mask(it.index);
        if(b) Add_Neighbors(m,next_outside,neighbors_outside,it.index,-1,-3);
        else Add_Neighbors(m,next_inside,neighbors_inside,it.index,-2,-4);}

    VECTOR<T,3> col(1,0,0);

    for(int j=0;j<current.m;j++)
        Add_Debug_Particle(grid.Node(current(j)),VECTOR<T,3>(0,0,1));

    for(int o=0;o<fill_width;o++){
        current.Exchange(next_outside);
        for(int j=0;j<current.m;j++){
            m.node_to_index(current(j))=m.index_to_node.Append(current(j));
            Add_Neighbors(m,next_outside,neighbors_outside,current(j),-1,-3);}

        for(int j=0;j<next_outside.m;j++)
            Add_Debug_Particle(grid.Node(next_outside(j)),col);

        current.Remove_All();
        exchange(col.x,col.y);}
    next_outside.Remove_All();
    Add_Debug_Particle(grid.Node(m.node_to_index.domain.min_corner),VECTOR<T,3>(1,1,1));
    Add_Debug_Particle(grid.Node(m.node_to_index.domain.max_corner),VECTOR<T,3>(1,1,1));
    Flush_Frame<TV>("register");

    for(int i=1;i<m.index_to_node.m;i++){
        TV N=phi.Normal(grid.X(m.index_to_node(i)));
        normal.Append(N);
        for(int d=0;d<TV::m;d++){
            TV_INT index=m.index_to_node(i);
            int s=N(d)<0?1:-1;
            for(int j=0;j<3;j++){
                int& n=m.node_to_index(index);
                if((n&-3)==-3) n=m.index_to_node.Append(index);
                index(d)+=s;}}}

    for(int o=0;o<order;o++){
        m.solve_indices(o)=INTERVAL<int>(1,m.index_to_node.m);
        m.fill_indices(o).min_corner=m.index_to_node.m;
        current.Exchange(next_inside);
        for(int j=0;j<current.m;j++){
            m.node_to_index(current(j))=m.index_to_node.Append(current(j));
            Add_Neighbors(m,next_inside,neighbors_inside,current(j),-2,-4);}
        current.Remove_All();}

    for(int o=order-1;o>=0;o--){
        current.Exchange(next_inside);
        for(int j=0;j<current.m;j++){
            m.node_to_index(current(j))=m.index_to_node.Append(current(j));
            Add_Neighbors(m,next_inside,neighbors_inside,current(j),-2,-4);}
        current.Remove_All();
        m.fill_indices(o).max_corner=m.index_to_node.m;}
    next_inside.Remove_All();

    for(int i=normal.m;i<m.fill_indices(1).max_corner;i++)
        normal.Append(phi.Normal(grid.X(m.index_to_node(i))));

    VECTOR<STENCIL,TV::m> st;
    stencil.Append(st);
    for(int i=m.solve_indices(order-1).min_corner;i<m.solve_indices(order-1).max_corner;i++){
        TV N=normal(i);
        for(int d=0;d<TV::m;d++){
            TV_INT index=m.index_to_node(i);
            int s=N(d)<0?1:-1;
            for(int j=0;j<3;j++){
                st(d).nodes(j+1)=m.node_to_index(index);
                st(d).scale=s*N(d)*grid.one_over_dX(d);
                index(d)+=s;}
            index(d)-=4*s;
            int n=m.node_to_index(index);
            st(d).nodes(0)=n==-3?0:n;}
        stencil.Append(st);}
}
//#####################################################################
// Function Fill_un
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Fill_un(const MAPPING& m,const TV& one_over_dx,const ARRAY<TV>& normal,const ARRAY<T2>& x,ARRAY<T2>& xn,int o,int mo)
{
    xn.Resize(m.fill_indices(o).max_corner);
    for(int i=m.fill_indices(o).min_corner;i<m.fill_indices(o).max_corner;i++){
        const TV_INT& index=m.index_to_node(i);
        T2 v=0;
        for(int d=0;d<TV::m;d++){
            TV_INT a=index,b=index;
            a(d)--;
            b(d)++;
            v+=(x(m.node_to_index(b))-x(m.node_to_index(a)))*(T).5*one_over_dx(d)*normal(i)(d);}
        xn(i)=v;}
}
//#####################################################################
// Function Constant_Extrapolate_FE
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_FE(const MAPPING& m,const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,const ARRAY<T2>& u,ARRAY<T2>& y,const ARRAY<T2>* z,int o,T dt,T alpha)
{
    for(int i=m.solve_indices(o).min_corner;i<m.solve_indices(o).max_corner;i++){
        T2 dot=0,zi=z?(*z)(i):0,a=0;
        for(int d=0;d<TV::m;d++){
            // Second order ENO.
            const STENCIL& s=stencil(i)(d);
            VECTOR<T2,4> f(u(s.nodes(0)),u(s.nodes(1)),u(s.nodes(2)),u(s.nodes(3)));
            VECTOR<T2,3> df(f(1)-f(0),f(2)-f(1),f(3)-f(2));
            VECTOR<T2,2> ddf(df(1)-df(0),df(2)-df(1));
            if(abs(ddf(0))<abs(ddf(1))) a=(T).5*(f(2)-f(0));
            else a=df(1)-(T).5*ddf(1);
            dot+=a*s.scale;}
        y(i)+=alpha*(u(i)-y(i)-dt*(dot-zi));}
}
//#####################################################################
// Function Constant_Extrapolate_RK2
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_RK2(const MAPPING& m,const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,ARRAY<T2>& u,const ARRAY<T2>* z,ARRAY<T2>& tmp,int o,T dt)
{
    tmp.Resize(u.m);
    for(int i=m.fill_indices(o).min_corner;i<m.fill_indices(o).max_corner;i++) tmp(i)=u(i);
    Extrapolate_FE(m,stencil,u,tmp,z,o,dt,1);
    Extrapolate_FE(m,stencil,tmp,u,z,o,dt,(T).5);
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Node(const GRID<TV>& grid,const LEVELSET<TV>& phi,boost::function<bool(const TV_INT& index)> inside_mask,
    int ghost,ARRAYS_ND_BASE<T2,TV_INT>& u,int iterations,int order,int fill_width)
{
    PHYSBAM_ASSERT(order>=1 && order<=3);
    PHYSBAM_ASSERT(!grid.Is_MAC_Grid());
    T dt=grid.dX.Max()/(TV::m+1);
    MAPPING m;
    ARRAY<TV> normal;
    ARRAY<VECTOR<STENCIL,TV::m> > stencil;
    Register_Nodes(grid,phi,inside_mask,ghost,m,normal,stencil,order,fill_width);
    int num_indices=m.fill_indices(0).max_corner;
    ARRAY<T2> du[3],tmp(num_indices);
    du[0].Resize(num_indices);
    for(int i=m.fill_indices(0).min_corner;i<m.fill_indices(0).max_corner;i++) du[0](i)=u(m.index_to_node(i));
    for(int o=1;o<order;o++) Fill_un(m,grid.one_over_dX,normal,du[o-1],du[o],o,order);
    for(int i=0;i<order;i++) du[i](0)=FLT_MAX/100; // Sentinal values for ENO.
    tmp(0)=FLT_MAX/100;
    for(int o=order-1;o>=0;o--) for(int i=0;i<iterations;i++) Extrapolate_RK2(m,stencil,du[o],(o!=order-1?&du[o+1]:0),tmp,o,dt);
    for(int i=m.solve_indices(0).min_corner;i<m.solve_indices(0).max_corner;i++) u(m.index_to_node(i))=du[0](i);
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Cell(const GRID<TV>& grid,const LEVELSET<TV>& phi,boost::function<bool(const TV_INT& index)> inside_mask,
    int ghost,ARRAYS_ND_BASE<T2,TV_INT>& u,int iterations,int order,int fill_width)
{
    GRID<TV> node_grid(grid.Get_Regular_Grid_At_MAC_Positions());
    Extrapolate_Node(node_grid,phi,inside_mask,ghost,u,iterations,order,fill_width);
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Face(const GRID<TV>& grid,const LEVELSET<TV>& phi,boost::function<bool(const FACE_INDEX<TV::m>& index)> inside_mask,
    int ghost,ARRAY<T2,FACE_INDEX<TV::m> >& u,int iterations,int order,int fill_width)
{
    for(int i=0;i<TV::m;i++){
        GRID<TV> node_grid(grid.Get_Face_Grid(i));
        Extrapolate_Node(node_grid,phi,[&](const TV_INT& index){return inside_mask(FACE_INDEX<TV::m>(i,index));},
            ghost,u.Component(i),iterations,order,fill_width);}
}
//#####################################################################
// Function Smooth_With_Heat_Equation
//#####################################################################
template<class T,class TV,class T2,class TV_INT> void
Smooth_With_Heat_Equation(const GRID<TV>& grid,int ghost,ARRAY<T2,TV_INT>& phi,T frac,int n)
{
    ARRAY<T2,TV_INT> tmp(phi.domain);
    for(int i=0;i<n;i++){
        for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost-1);it.Valid();it.Next()){
            T x=0;
            for(int d=0;d<TV::m;d++){
                TV_INT a(it.index),b(it.index);
                a(d)--;
                b(d)++;
                x+=phi(a)+phi(b);}
            tmp(it.index)=x;}
        phi.Copy(frac/(2*TV::m),tmp,1-frac,phi);

        Dump_Levelset(grid,phi,true,VECTOR<T,3>(1,0,0));
        Flush_Frame<TV>("after heat");}
}
//#####################################################################
// Function Extrapolate_Node_No_Levelset
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Node_No_Levelset(const GRID<TV>& grid,boost::function<bool(const TV_INT& index)> inside_mask,int ghost,ARRAYS_ND_BASE<T2,TV_INT>& u,int iterations,int order,int fill_width,int smooth_steps)
{
    ARRAY<T,TV_INT> phi(grid.Node_Indices(ghost+1),true,1);
    ARRAY<TV_INT> seed_indices;
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next())
        if(inside_mask(it.index))
            phi(it.index)=-1;
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next())
        if(phi(it.index)==-1)
            for(int d=0;d<TV::m;d++){
                TV_INT a(it.index),b(it.index);
                a(d)--;
                b(d)++;
                if(phi(a)==1 || phi(b)==1){
                    phi(it.index)=0;
                    seed_indices.Append(it.index);
                    break;}}

    LEVELSET<TV> levelset(const_cast<GRID<TV>&>(grid),phi,ghost+1);
    Dump_Levelset(grid,phi,true,VECTOR<T,3>(1,0,0));
    // phi.array+=grid.dX.Max()*6;
    // Dump_Levelset(grid,phi,VECTOR<T,3>(0,1,0));
    // phi.array-=grid.dX.Max()*12;
    // Dump_Levelset(grid,phi,VECTOR<T,3>(0,0,1));
    // phi.array+=grid.dX.Max()*6;

    Flush_Frame<TV>("after init");

    Reinitialize(levelset,(ghost+2)*20,(T)0,(ghost+2)*grid.dX.Max()*20,(T)0,(T).9,3,5,0);

    Dump_Levelset(grid,phi,true,VECTOR<T,3>(1,0,0));
    // phi.array+=grid.dX.Max()*6;
    // Dump_Levelset(grid,phi,VECTOR<T,3>(0,1,0));
    // phi.array-=grid.dX.Max()*12;
    // Dump_Levelset(grid,phi,VECTOR<T,3>(0,0,1));
    // phi.array+=grid.dX.Max()*6;

    Flush_Frame<TV>("after reinit");

    Smooth_With_Heat_Equation(grid,ghost,phi,(T).8,smooth_steps);

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost-1);it.Valid();it.Next()){
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,levelset.Normal(it.Location()));}
    Flush_Frame<TV>("normals");

    Reinitialize(levelset,(ghost+2)*20,(T)0,(ghost+2)*grid.dX.Max()*20,(T)0,(T).9,3,5,0);
    Dump_Levelset(grid,phi,true,VECTOR<T,3>(1,0,0));
    // phi.array+=grid.dX.Max()*6;
    // Dump_Levelset(grid,phi,VECTOR<T,3>(0,1,0));
    // phi.array-=grid.dX.Max()*12;
    // Dump_Levelset(grid,phi,VECTOR<T,3>(0,0,1));
    // phi.array+=grid.dX.Max()*6;

    Flush_Frame<TV>("after reinit");

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost-1);it.Valid();it.Next()){
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,levelset.Normal(it.Location()));}
    Flush_Frame<TV>("normals");

    Extrapolate_Node(grid,levelset,inside_mask,ghost,u,iterations,order,fill_width);
}
//#####################################################################
// Function Extrapolate_Cell_No_Levelset
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Cell_No_Levelset(const GRID<TV>& grid,boost::function<bool(const TV_INT& index)> inside_mask,int ghost,ARRAYS_ND_BASE<T2,TV_INT>& u,int iterations,int order,int fill_width,int smooth_steps)
{
    GRID<TV> node_grid(grid.Get_Regular_Grid_At_MAC_Positions());
    Extrapolate_Node_No_Levelset(node_grid,inside_mask,ghost,u,iterations,order,fill_width,smooth_steps);
}
//#####################################################################
// Function Extrapolate_Face_No_Levelset
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Face_No_Levelset(const GRID<TV>& grid,boost::function<bool(const FACE_INDEX<TV::m>& index)> inside_mask,int ghost,ARRAY<T2,FACE_INDEX<TV::m> >& u,int iterations,int order,int fill_width,int smooth_steps)
{
    ARRAY<T,TV_INT> phi(grid.Domain_Indices(ghost+1),true,0);
    ARRAY<TV_INT> seed_indices;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next())
        if(inside_mask(it.Full_Index())){
            phi(it.First_Cell_Index())++;
            phi(it.Second_Cell_Index())++;}
    T scale=grid.dX.Max()/TV::m;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        T& p=phi(it.index);
        if(p!=0 && p!=TV::m*2) seed_indices.Append(it.index);
        p=(p-TV::m)*scale;}

    LEVELSET<TV> levelset(const_cast<GRID<TV>&>(grid),phi,ghost+1);
    FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(levelset,ghost+2);
    fmm.Fast_Marching_Method(phi,grid.dX.Max()*(ghost+2),&seed_indices);
    Reinitialize(levelset,(ghost+2)*20,(T)0,(ghost+2)*grid.dX.Max()*20,(T)10,(T).9,3,5,0);
    Smooth_With_Heat_Equation(grid,ghost,phi,(T).8,smooth_steps);
    Reinitialize(levelset,(ghost+2)*20,(T)0,(ghost+2)*grid.dX.Max()*20,(T)10,(T).9,3,5,0);

    Extrapolate_Face(grid,levelset,inside_mask,ghost,u,iterations,order,fill_width);
}
namespace PhysBAM{
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,1>,float>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,2>,float>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,3>,float>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,1>,double>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,2>,double>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,3>,double>;
}
