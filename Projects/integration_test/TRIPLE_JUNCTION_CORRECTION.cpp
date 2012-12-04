#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "TRIPLE_JUNCTION_CORRECTION.h"
using namespace PhysBAM;
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
    HASHTABLE<TRIPLE<int,int,TV_INT>,T> total_size;
    for(int i=0;i<stencils.m;i++){
        int mask=~0;
        for(int j=0;j<stencils(i).m;j++)
            mask&=pairwise_data(stencils(i)(j)).valid_flags;
        for(int a=0;a<phi.m;a++)
            if(mask&(1<<a))
                for(int b=a+1;b<phi.m;b++)
                    if(mask&(1<<b))
                        for(int c=b+1;c<phi.m;c++)
                            if(mask&(1<<c)){
                                VECTOR<T,TV::m+1> ab(pairwise_phi(a)(b).Subset(stencils(i)));
                                VECTOR<T,TV::m+1> bc(pairwise_phi(b)(c).Subset(stencils(i)));
                                VECTOR<T,TV::m+1> ca(-pairwise_phi(a)(c).Subset(stencils(i)));
                                T A=Bad_Fraction(VECTOR<VECTOR<T,TV::m+1>,3>(ab,bc,ca));
                                for(int j=0;j<stencils(i).m;j++){
                                    total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(a,b,stencils(i)(j)))+=A;
                                    total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(b,c,stencils(i)(j)))+=A;
                                    total_size.Get_Or_Insert(TRIPLE<int,int,TV_INT>(a,c,stencils(i)(j)))-=A;}}}
    for(typename HASHTABLE<TRIPLE<int,int,TV_INT>,T>::ITERATOR it(total_size);it.Valid();it.Next()){
        pairwise_phi(it.Key().x)(it.Key().y)(it.Key().z)-=(T).1*/*Weight_Function*/(it.Data()*grid.dX.Max());
        LOG::cout<<it.Data()<<" "<<grid.dX<<" "<<(T).1*/*Weight_Function*/(it.Data()*grid.dX.Max())<<std::endl;}
}
//#####################################################################
// Function Update_Color_Level_Sets
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Update_Color_Level_Sets()
{
    ARRAY<T> min_phi(phi.m);
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost);it.Valid();it.Next()){
        min_phi.Fill(FLT_MAX);
        int mask=pairwise_data(it.index).valid_flags,color=-1,color_mask=0,not_color_mask=0;
        int num_bits=count_bits(mask);
        if(num_bits==2 && pairwise_data(it.index).trust.Contains(-1)){
            int rmb=rightmost_bit(mask),lmb=mask-rmb,a=integer_log_exact(rmb),b=integer_log_exact(lmb);
            T p=pairwise_phi(a)(b)(it.index);
            phi(a)(it.index)=p;
            phi(b)(it.index)=-p;}
        if(count_bits(mask)<3) continue;
        for(int a=0;a<phi.m;a++)
            if(mask&(1<<a))
                for(int b=a+1;b<phi.m;b++)
                    if(mask&(1<<b))
                        for(int c=b+1;c<phi.m;c++)
                            if(mask&(1<<c)){
                                int new_color=-1;
                                T pab=pairwise_phi(a)(b)(it.index),pac=pairwise_phi(a)(c)(it.index),pbc=pairwise_phi(b)(c)(it.index),new_phi=0;
                                if(pab<=0 && pac<=0){
                                    new_phi=max(pab,pac);
                                    new_color=a;}
                                else if(pab>0 && pbc<=0){
                                    new_phi=max(-pab,pbc);
                                    new_color=b;}
                                else if(pbc>0 && pac>0){
                                    new_phi=max(-pbc,-pac);
                                    new_color=c;}
                                else continue;
                                min_phi(new_color)=max(min_phi(new_color),new_phi);
                                color_mask|=1<<new_color;
                                not_color_mask|=((1<<a)|(1<<b)|(1<<c))&~(1<<new_color);}
        color_mask&=~not_color_mask;
        PHYSBAM_ASSERT(power_of_two(color_mask));
        color=integer_log(color_mask);

        for(int a=0;a<color;a++)
            phi(a)(it.index)=pairwise_phi(a)(color)(it.index);
        for(int a=color+1;a<phi.m;a++)
            phi(a)(it.index)=-pairwise_phi(color)(a)(it.index);
        phi(color)(it.index)=min_phi(color);}
}
//#####################################################################
// Function Initialize_Stencils_Helper
//#####################################################################
template<class TV,class TV_INT> static void
Initialize_Stencils_Helper(const GRID<TV>& grid,ARRAY<VECTOR<TV_INT,3> >& stencils)
{
    TRIANGULATED_AREA<typename TV::SCALAR> ta;
    ta.Initialize_Square_Mesh_And_Particles(grid,true);

    stencils.Resize(ta.mesh.elements.m);
    for(int i=0;i<ta.mesh.elements.m;i++)
        for(int v=0;v<TV::m+1;v++)
            stencils(i)(v)=TV_INT(rint(ta.particles.X(ta.mesh.elements(i)(v))));
}
//#####################################################################
// Function Initialize_Stencils_Helper
//#####################################################################
template<class TV,class TV_INT> static void
Initialize_Stencils_Helper(const GRID<TV>& grid,ARRAY<VECTOR<TV_INT,4> >& stencils)
{
    TETRAHEDRALIZED_VOLUME<typename TV::SCALAR> tv;
    tv.Initialize_Cube_Mesh_And_Particles(grid);

    stencils.Resize(tv.mesh.elements.m);
    for(int i=0;i<tv.mesh.elements.m;i++)
        for(int v=0;v<TV::m+1;v++)
            stencils(i)(v)=TV_INT(rint(tv.particles.X(tv.mesh.elements(i)(v))));
}
//#####################################################################
// Function Initialize_Stencils
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Initialize_Stencils()
{
    GRID<TV> grid_int(grid.Numbers_Of_Cells(),RANGE<TV>(TV(),TV(grid.Numbers_Of_Cells())),false);
    Initialize_Stencils_Helper(grid_int,stencils);
}
//#####################################################################
// Function Initialize_Stencils
//#####################################################################
template<class TV> typename TV::SCALAR TRIPLE_JUNCTION_CORRECTION<TV>::
Bad_Fraction(const VECTOR<VECTOR<T,TV::m+1>,3>& phi)
{
    int pos_mask=0,neg_mask=0;
    for(int c=0;c<3;c++)
        for(int j=0;j<TV::m+1;j++){
            if(phi(c)(j)>0) pos_mask|=1<<c;
            else if(phi(c)(j)<0) neg_mask|=1<<c;}

    if(pos_mask==7){
        if(neg_mask==0) return 1;}
    else if(neg_mask==7){
        if(pos_mask==0) return -1;}
    else return 0;

    for(int c=0;c<3;c++)
        for(int j=0;j<TV::m+1;j++)
            if(T pj=phi(c)(j))
                for(int k=j+1;k<TV::m+1;k++)
                    if(T pk=phi(c)(k))
                        if((pj>0)!=(pk>0)){
                            VECTOR<VECTOR<T,TV::m+1>,3> phi0=phi,phi1=phi;
                            T th=pj/(pj-pk);
                            for(int e=0;e<3;e++){
                                T ej=phi(e)(j),ek=phi(e)(k),el=ej+th*(ek-ej);
                                phi0(e)(j)=ej;
                                phi0(e)(k)=el;
                                phi1(e)(j)=el;
                                phi1(e)(k)=ek;}
                            phi0(c)(k)=phi1(c)(j)=0;
                            return th*Bad_Fraction(phi0)+(1-th)*Bad_Fraction(phi1);}

    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Initialize_Stencils
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Cut_Interface(HASHTABLE<TV_INT,HASH_CELL_DATA>& index_to_cell_data)
{
    for(int i=0;i<stencils.m;i++){
        int mask=~0;
        for(int j=0;j<stencils(i).m;j++)
            mask&=pairwise_data(stencils(i)(j)).valid_flags;
        TV_INT cell(stencils(i)(0));
        for(int j=1;j<stencils(i).m;j++)
            cell=cell.Componentwise_Min(stencils(i)(j));
        if(count_bits(mask)<3) Cut_Stencil_With_Phi(index_to_cell_data,cell,i);
        else Cut_Stencil_With_Pairwise_Phi(index_to_cell_data,cell,i);}
}
//#####################################################################
// Function Cut_Stencil_With_Phi
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Cut_Stencil_With_Phi(HASHTABLE<TV_INT,HASH_CELL_DATA>& index_to_cell_data,const TV_INT& cell,int s)
{
    
}
//#####################################################################
// Function Cut_Stencil_With_Pairwise_Phi
//#####################################################################
template<class TV> void TRIPLE_JUNCTION_CORRECTION<TV>::
Cut_Stencil_With_Pairwise_Phi(HASHTABLE<TV_INT,HASH_CELL_DATA>& index_to_cell_data,const TV_INT& cell,int s)
{
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
template class TRIPLE_JUNCTION_CORRECTION<VECTOR<float,3> >;
template class TRIPLE_JUNCTION_CORRECTION<VECTOR<double,2> >;
template class TRIPLE_JUNCTION_CORRECTION<VECTOR<double,3> >;
}
