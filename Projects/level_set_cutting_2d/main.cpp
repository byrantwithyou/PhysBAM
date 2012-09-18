#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Images/EPS_FILE.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "MARCHING_TETRAHEDRA_CUTTING.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef VECTOR<int,2> TV_INT;
    typedef VECTOR<T,2> TV;

    int seed=time(0),size=4;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-s",&seed,"value","Seed");
    parse_args.Add("-n",&size,"value","Domain size");
    parse_args.Parse();

#if 0
    ARRAY<VECTOR<int,3> > m;
    m.Append(VECTOR<int,3>(0,1,2));
    ARRAY<T> p0;
    p0.Append(0.5);
    p0.Append(0);
    p0.Append(-0.1);
    ARRAY<TV> X(3);
    X(0)=TV(0,0);
    X(1)=TV(1,0);
    X(2)=TV(0,1);
    ARRAY<VECTOR<int,3> > c=m,sp;
    ARRAY<PAIR<VECTOR<int,2>,T> > weights;
    LOG::cout<<m<<"    "<<c<<std::endl;
    MARCHING_TETRAHEDRA_CUTTING<TV>::Query_Case(m,c,sp,p0,weights);
    LOG::cout<<m<<"    "<<c<<"    "<<sp<<"    "<<weights<<std::endl;
    for(int i=0;i<weights.m;i++)
        X.Append(X(weights(i).x.x)*(1-weights(i).y)+X(weights(i).x.y)*weights(i).y);

    EPS_FILE<T> eps("out.eps");
    for(int i=0;i<c.m;i++)
        eps.Draw_Object(X(c(i)(0)),X(c(i)(1)),X(c(i)(2)));
    p0.Resize(X.m);
    for(int i=0;i<p0.m;i++){
        eps.cur_format.fill_color=VECTOR<T,3>(p0(i)<0,p0(i)>0,0);
        eps.cur_format.fill_style=1;
        eps.cur_format.line_style=0;
        eps.Draw_Object(X(i),(T).02);}



    // TETRAHEDRON_MESH tm;
    // for(int i=0;i<cut_mesh.m;i++){
    //     tm.elements.Append(cut_mesh(i).indices);
    //     LOG::cout<<cut_mesh(i).parent<<"  "<<cut_mesh(i).indices<<"  "<<cut_mesh(i).weights<<std::endl;}

    // tm.elements=m;
    // tm.Set_Number_Nodes(tm.elements.Flattened().Max()+1);
    // tm.Initialize_Boundary_Mesh();

    // LOG::cout<<tm.boundary_mesh->elements<<std::endl;

    return 0;

#else

    GRID<TV> grid(TV_INT()+(size+1),RANGE<TV>::Unit_Box());

    TRIANGULATED_AREA<T>& tv=*TRIANGULATED_AREA<T>::Create();
    tv.Initialize_Herring_Bone_Mesh_And_Particles(grid);

    // int k=0;
    // for(int i=0;i<tv.mesh.elements.m;i++)
    //     if(tv.mesh.elements(i).Contains(13))
    //         tv.mesh.elements(k++)=tv.mesh.elements(i);
    // tv.mesh.elements.Resize(k);

//    LOG::cout<<tv.mesh.elements<<std::endl;
//    tv.mesh.elements.Remove_All();
//    tv.mesh.elements.Append(VECTOR<int,4>(30,31,28,19));


    tv.mesh.Initialize_Boundary_Mesh();
    tv.mesh.Initialize_Boundary_Nodes();
    printf("boundary: %i\n",tv.mesh.boundary_mesh->elements.m);

    ARRAY<T> phi0(tv.particles.number),phi1(tv.particles.number);
    RANDOM_NUMBERS<T> random;
    printf("seed %i\n", seed);
    random.Set_Seed(seed);
    random.Fill_Uniform(phi0,-5,5);
    random.Fill_Uniform(phi1,-5,5);
    for(int i=0;i<phi0.m;i++) phi0(i)=::PhysBAM::rint(phi0(i));
    for(int i=0;i<phi1.m;i++) phi1(i)=::PhysBAM::rint(phi1(i));
    phi0.Subset(*tv.mesh.boundary_nodes).Fill(1);
    phi1.Subset(*tv.mesh.boundary_nodes).Fill(1);

    ARRAY<VECTOR<int,3> > m=tv.mesh.elements,sp;
    ARRAY<PAIR<VECTOR<int,2>,T> > weights;
    MARCHING_TETRAHEDRA_CUTTING<TV>::Query_Case(m,tv.mesh.elements,sp,phi0,weights);
    tv.particles.Add_Elements(weights.m);
    for(int i=0;i<weights.m;i++)
        tv.particles.X(i+phi0.m)=tv.particles.X(weights(i).x.x)*(1-weights(i).y)+tv.particles.X(weights(i).x.y)*weights(i).y;

    EPS_FILE<T> eps("out.eps");
    for(int i=0;i<tv.mesh.elements.m;i++)
        eps.Draw_Object(tv.particles.X(tv.mesh.elements(i)(0)),tv.particles.X(tv.mesh.elements(i)(1)),tv.particles.X(tv.mesh.elements(i)(2)));
    phi0.Resize(tv.particles.X.m);
    for(int i=0;i<phi0.m;i++){
        eps.cur_format.fill_color=VECTOR<T,3>(phi0(i)<0,phi0(i)>0,phi0(i)==0);
        eps.cur_format.fill_style=1;
        eps.cur_format.line_style=0;
        eps.Draw_Object(tv.particles.X(i),.01);}

    tv.Update_Number_Nodes();
    tv.mesh.Initialize_Boundary_Nodes();
    printf("boundary: %i\n",tv.mesh.boundary_mesh->elements.m);

    // TETRAHEDRON_MESH tm;
    // for(int i=0;i<cut_mesh.m;i++)
    //     tm.elements.Append(cut_mesh(i).indices);

    // tm.Set_Number_Nodes(tm.elements.Flattened().Max()+1);
    // tm.Initialize_Boundary_Mesh();
    // LOG::cout<<"mesh sizes "<<tm.boundary_mesh->elements.m<<"  "<<tv.mesh.boundary_mesh->elements.m<<std::endl;
    // // LOG::cout<<tm.elements<<std::endl;
    // LOG::cout<<tm.boundary_mesh->elements<<std::endl;
    // LOG::cout<<tv.mesh.boundary_mesh->elements<<std::endl;

    // for(int i=0;i<cut_mesh.m;i++){
    //     VECTOR<TV,4> pts;
    //     for(int j=0;j<4;j++)
    //         for(int k=0;k<4;k++)
    //             pts(j)+=tv.particles.X(cut_mesh(i).parent(k))*cut_mesh(i).weights(j)(k);
    //     if(TETRAHEDRON<T>::Signed_Size(pts)<-1e-15)
    //         LOG::cout<<"NEGATIVE VOLUME  "<<TETRAHEDRON<T>::Signed_Size(pts)<<"   "<<cut_mesh(i).parent<<"   "<<cut_mesh(i).indices<<std::endl;}

    // for(int i=0;i<cut_mesh.m;i++){
    //     VECTOR<T,4> a(phi0.Subset(cut_mesh(i).indices)),b(phi1.Subset(cut_mesh(i).indices));
    //     if(a.Min()<0 && a.Max()>0 && b.Min()<0 && b.Max()>0)
    //         LOG::cout<<"CROSSINGS "<<a<<"  "<<b<<std::endl;}

    // delete &tv;
#endif
    return 0;
}
