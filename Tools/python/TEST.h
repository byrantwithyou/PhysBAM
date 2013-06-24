#ifndef _TEST_H_
#define _TEST_H_
#include <Tools/Log/LOG.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Dynamics/Particles/DEFORMABLE_BODY_PARTICLES.h>

using namespace PhysBAM;

void Initialize_Incident_Old(TRIANGLE_MESH& mesh)
{
    LOG::Time("Build");
    mesh.Initialize_Incident_Elements();
    LOG::Time("Query Sequential");
    for(int e=0;e<mesh.elements.m;e++){
        const VECTOR<int,3>& nodes=mesh.elements(e);
        for(int i=0;i<3;i++) 
            if(!(*mesh.incident_elements)(nodes[i]).Contains(e)) PHYSBAM_FATAL_ERROR();}
}

void Initialize_Incident_New(TRIANGLE_MESH& mesh)
{
    LOG::Time("Build");
    int min_index=2000000000;
    int max_index=0;

    for(int e=0;e<mesh.elements.m;e++){
        const VECTOR<int,3>& nodes=mesh.elements(e);
        for(int i=0;i<3;i++){
            min_index=min(min_index,nodes[i]);
            max_index=max(max_index,nodes[i]);}}
    LOG::cout<<"min="<<min_index<<" max="<<max_index<<std::endl;
    ARRAY<ARRAY<int> > reverse(max_index-min_index+1);
    for(int e=0;e<mesh.elements.m;e++){
        const VECTOR<int,3>& nodes=mesh.elements(e);
        for(int i=0;i<3;i++){
            reverse(nodes[i]-min_index+1).Append(e);}}
    LOG::Time("Query Sequential");
    int base_index=min_index-1;
    for(int e=0;e<mesh.elements.m;e++){
        const VECTOR<int,3>& nodes=mesh.elements(e);
        for(int i=0;i<3;i++){
            int lookup=nodes[i]-base_index;
            if(!reverse(lookup).Contains(e)) PHYSBAM_FATAL_ERROR();}}

}

void Initialize_Incident_Hash(TRIANGLE_MESH& mesh)
{
    HASHTABLE<int,ARRAY<int> > hash(mesh.elements.m);

    LOG::Time("Build");
    for(int e=0;e<mesh.elements.m;e++){
        const VECTOR<int,3>& nodes=mesh.elements(e);
        for(int i=0;i<3;i++) 
            hash.Get_Or_Insert(nodes[i]).Append(e);}
    LOG::Time("Query Sequential");
    for(int e=0;e<mesh.elements.m;e++){
        const VECTOR<int,3>& nodes=mesh.elements(e);
        for(int i=0;i<3;i++) 
            if(!hash.Get(nodes[i]).Contains(e)) PHYSBAM_FATAL_ERROR();}
}

template<class TV>
void Scale_Restlength(LINEAR_SPRINGS<TV>& linear_springs,const typename TV::SCALAR scale)
{
    linear_springs.Invalidate_CFL();
    for(int i=0;i<linear_springs.segment_mesh.elements.m;i++){
        linear_springs.restlength(i)=linear_springs.visual_restlength(i)=linear_springs.restlength(i)*scale;}
}

template<class TV>
void Project(const ARRAY<int>& indices,RIGID_BODY<TV>& body,DEFORMABLE_BODY_PARTICLES<TV>& particle)
{
    for(int i=0;i<indices.m;i++){
        TV& X=particle.X(indices(i));
        if(!body.Oriented_Bounding_Box().Lazy_Inside(X)) continue;
        typename TV::SCALAR phi;
        TV normal=body.Implicit_Geometry_Extended_Normal(X,phi);
        if(TV::Dot_Product(normal,TV(0,-1,0))<.3) continue;
        X+=phi*normal;
    }
}


#endif
