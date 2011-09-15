//#####################################################################
// Copyright 2002, 2003, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CAPE_COLLISIONS
//##################################################################### 
#ifndef __CAPE_COLLISIONS__
#define __CAPE_COLLISIONS__

namespace PhysBAM{

template<class T>
class CAPE_COLLISIONS
{
public:
    TRIANGLE_MESH &triangle_mesh_1,&triangle_mesh_2;
    PARTICLES<T,VECTOR_3D<T> > &particles_1,&particles_2;
    TRIANGULATED_SURFACE<T> &surface_1,&surface_2;
    TRIANGLE_HIERARCHY<T> *hierarchy_1,*hierarchy_2;
    ARRAY<T> distances;
    ARRAY<int> close_triangles;
    ARRAY<VECTOR_3D<T> > close_points,close_weights;
    T thickness,tolerance;
    T bound;
    T coefficient_of_friction;
    ARRAY<bool> *disabled;
    
    ARRAY<bool> *enforce_collision_velocity;
    ARRAY<VECTOR_3D<T> > *collision_normal;
    ARRAY<T> *collision_velocity;
    
    T ground_thickness;
    
    /*
    T flypaper_thickness;
    ARRAY<bool> flypaper;
    ARRAY<VECTOR_3D<T> > flypaper_X;
    int flies;
    */
    
    int debug;
    
    CAPE_COLLISIONS(TRIANGULATED_SURFACE<T>& surface_1_input,TRIANGULATED_SURFACE<T>& surface_2_input)
        :triangle_mesh_1(surface_1_input.triangle_mesh),triangle_mesh_2(surface_2_input.triangle_mesh),
        particles_1(surface_1_input.particles),particles_2(surface_2_input.particles),
        surface_1(surface_1_input),surface_2(surface_2_input),hierarchy_1(0),hierarchy_2(0),
        distances(1,particles_1.number),close_triangles(1,particles_1.number),close_points(1,particles_1.number),close_weights(1,particles_1.number),
        thickness((T)1e-3),tolerance((T)1e-7),bound(0),coefficient_of_friction(0),disabled(0),enforce_collision_velocity(0),collision_normal(0),collision_velocity(0),debug(0),
        ground_thickness(0)
        //flypaper_thickness(0),flypaper(1,particles_1.number),flypaper_X(1,particles_1.number)
    {
        if(!surface_1.triangle_hierarchy)surface_1.Initialize_Triangle_Hierarchy();
        if(!surface_2.triangle_hierarchy)surface_2.Initialize_Triangle_Hierarchy();
        hierarchy_1=surface_1.triangle_hierarchy;hierarchy_2=surface_2.triangle_hierarchy;
        if(!triangle_mesh_2.incident_triangles)triangle_mesh_2.Initialize_Incident_Triangles();
        if(!triangle_mesh_2.adjacent_triangles)triangle_mesh_2.Initialize_Adjacent_Triangles();
    }
    
    void Enable_Constraints(ARRAY<bool>& enforce_collision_velocity_input,ARRAY<VECTOR_3D<T> >& collision_normal_input,ARRAY<T>& collision_velocity_input)
    {enforce_collision_velocity=&enforce_collision_velocity_input;collision_normal=&collision_normal_input;collision_velocity=&collision_velocity_input;}
    
    void Disable_Currently_Colliding_Particles()
    {if(!disabled)disabled=new ARRAY<bool>();disabled->Resize_List(particles_1.number);
    Find_Closest_Points();
    ARRAY<int> queue;
    for(int p=1;p<=particles_1.number;p++)if(close_triangles(p) && Close_Point_Phi(p)<0){
        queue.Append(p);(*disabled)(p)=true;}
    bool neighbor_nodes_defined=triangle_mesh_1.neighbor_nodes!=0;if(!neighbor_nodes_defined)triangle_mesh_1.Initialize_Neighbor_Nodes();
    for(int i=1;i<=queue.m;i++){
        int p=queue(i);
        for(int a=1;a<=(*triangle_mesh_1.neighbor_nodes)(p).m;a++){
            int q=(*triangle_mesh_1.neighbor_nodes)(p)(a);
            if(!(*disabled)(q) && !close_triangles(q)){
                queue.Append(q);(*disabled)(q)=true;}}}
    if(!neighbor_nodes_defined){delete triangle_mesh_1.neighbor_nodes;triangle_mesh_1.neighbor_nodes=0;}}
    
    void Process_Collisions(T dt)
    {std::cout<<"cape collisions:\n";
    Estimate_Bound(dt);
    std::cout<<"    bound "<<bound<<"\n";
    std::cout<<"    finding closest points...\n";
    Find_Closest_Points();
    std::cout<<"    processing collisions...\n";
    // ARRAY<bool>::copy(false,flypaper);
    int close_count=0,interactions=0;//flies=0;
    for(int p=1;p<=particles_1.number;p++)if(close_triangles(p) && !(*enforce_collision_velocity)(p) && !(disabled && (*disabled)(p))){
        close_count++;interactions+=Process_Collision(dt,p);}
    std::cout<<"    close points = "<<close_count<<"\n";
    std::cout<<"    interactions = "<<interactions<<"\n";}
//    std::cout<<"    flies = "<<flies<<"\n";}
        
    int Process_Collision(T dt,int p)
    {if(particles_1.X(p).y<ground_thickness)return 0;
    VECTOR_3D<T> normal;
    T phi=Close_Point_Phi_And_Normal(p,normal);
    if(phi>0)return 0;T depth=-phi;
    VECTOR_3D<T> &X=particles_1.X(p),&V=particles_1.V(p);
    X+=depth*normal;
    T VN=VECTOR_3D<T>::Dot_Product(V,normal);VECTOR_3D<T> VT=V-VN*normal;
    VECTOR_3D<T> V_rigid_body=Surface_2_Velocity(close_triangles(p),close_weights(p));
    T VN_rigid_body=VECTOR_3D<T>::Dot_Product(V_rigid_body,normal);
    T VN_final=max(VN,VN_rigid_body);
    if(coefficient_of_friction){
        T normal_force=VN_final-VN+depth/dt;
        VECTOR_3D<T> VT_rigid_body=V_rigid_body-VN_rigid_body*normal,VT_relative=VT-VT_rigid_body;T VT_relative_magnitude=VT_relative.Magnitude();
        T friction_magnitude=1;if(coefficient_of_friction*normal_force < VT_relative_magnitude) 
            friction_magnitude=coefficient_of_friction*normal_force/VT_relative_magnitude;
        VT=VT_rigid_body+VT_relative*(1-friction_magnitude);}
    V=VN_final*normal+VT;
    assert(enforce_collision_velocity && collision_normal && collision_velocity);
    if(enforce_collision_velocity){   // constraints
        (*enforce_collision_velocity)(p)=true;
        (*collision_normal)(p)=normal;
        (*collision_velocity)(p)=VECTOR_3D<T>::Dot_Product(V,normal);}
    return 1;}
    
    VECTOR_3D<T> Surface_2_Velocity(int t,VECTOR_3D<T>& weights)
    {int i,j,k;triangle_mesh_2.triangles.Get(t,i,j,k);
    return TRIANGLE_3D<T>::Point_From_Barycentric_Coordinates(weights,particles_2.V(i),particles_2.V(j),particles_2.V(k));}
    
    void Estimate_Bound(T dt)
    {T velocity=particles_1.Maximum_Speed(0)+particles_2.Maximum_Speed(0);
    bound=2*dt*velocity+thickness;}
    
    void Find_Closest_Point_Base(const int p,const int t)
    {VECTOR_3D<T> x=particles_1.X(p);
    TRIANGLE_3D<T>& triangle=(*surface_2.triangle_list)(t);
    VECTOR_3D<T> w,close=triangle.Closest_Point(x,w);
    T distance=(close-x).Magnitude();
    if(distance<distances(p)){
        distances(p)=distance;close_triangles(p)=t;
        close_points(p)=close;close_weights(p)=w;}}
    
    void Find_Closest_Points_In_Triangles(const int t1,const int t2)
    {int i,j,k;triangle_mesh_1.triangles.Get(t1,i,j,k);
    Find_Closest_Point_Base(i,t2);Find_Closest_Point_Base(j,t2);Find_Closest_Point_Base(k,t2);}
    
    void Find_Closest_Points_In_Boxes(const int box1,const int box2)
    {assert(1<=box1 && box1<=hierarchy_1->root && 1<=box2 && box2<=hierarchy_2->root);
    if(!hierarchy_1->box_hierarchy(box1).Intersection(hierarchy_2->box_hierarchy(box2),bound))return;
    T radius_1=hierarchy_1->box_radius(box1);
    T radius_2=hierarchy_2->box_radius(box2);
    /*
    T center_distance=(hierarchy_1->box_hierarchy(box1).Center()-hierarchy_2->box_hierarchy(box2).Center()).Magnitude_Squared();
    if(radius_1+radius_2+bound<center_distance)return;
    */
    if(box1<=hierarchy_1->leaves){
        if(box2<=hierarchy_2->leaves)Find_Closest_Points_In_Triangles(box1,box2);
        else{for(int c=1;c<=2;c++)Find_Closest_Points_In_Boxes(box1,hierarchy_2->children(c,box2-hierarchy_2->leaves));}}
    else if(box2<=hierarchy_2->leaves || radius_1<radius_2)
        for(int c=1;c<=2;c++)Find_Closest_Points_In_Boxes(hierarchy_1->children(c,box1-hierarchy_1->leaves),box2);
    else for(int c=1;c<=2;c++)Find_Closest_Points_In_Boxes(box1,hierarchy_2->children(c,box2-hierarchy_2->leaves));}
    
    void Find_Closest_Points()
    {hierarchy_1->Update_Boxes();hierarchy_2->Update_Boxes();
    hierarchy_1->Update_Box_Radii();hierarchy_2->Update_Box_Radii();
    surface_2.Update_Triangle_List();
    ARRAY<T>::copy(bound,distances);ARRAY<int>::copy(0,close_triangles);
    Find_Closest_Points_In_Boxes(hierarchy_1->root,hierarchy_2->root);}
    
    bool Close_Point_Inside(int p)
    {return surface_2.Inside(particles_1.X(p));}

/*
    bool Close_Point_Inside(int p)
    {VECTOR_3D<T> location=particles_1.X(p);
    TRIANGLE_3D<T>& triangle=(*surface_2.triangle_list)(close_triangles(p));
    int region_id,region=triangle.Region(close_points(p),region_id,thickness);
    if(debug){std::cout<<"region "<<region<<"\n";}
    if(region == 1){ // vertex
        VECTOR_3D<T> direction=location-close_points(p);
        VECTOR_3D<T> point=close_points(p)+(T)1e-3*direction;
        if(debug){
            std::cout<<"cp "<<close_points(p)<<", distance "<<surface_2.Calculate_Signed_Distance(close_points(p),1e-7)<<"\n";
            std::cout<<location<<", distance "<<surface_2.Calculate_Signed_Distance(location,1e-7)<<"\n";
            std::cout<<point<<", distance "<<surface_2.Calculate_Signed_Distance(point,1e-7)<<"\n";}
        if(surface_2.Signed_Solid_Angle_Of_Triangle_Web(point,triangle_mesh_2.triangles(region_id,close_triangles(p))) > 0) return true;}
    else if(region == 2) { // edge
        int neighbor=(*triangle_mesh_2.adjacent_triangles)(close_triangles(p))(region_id);
        const TRIANGLE_3D<T>& triangle2=(*surface_2.triangle_list)(neighbor);
        int convex=0;
        if(region_id == 1){if(VECTOR_3D<T>::Dot_Product(triangle2.normal,triangle.x3-triangle2.x1) >= 0) convex=1;}
        else if(region_id == 2){if(VECTOR_3D<T>::Dot_Product(triangle2.normal,triangle.x1-triangle2.x1) >= 0) convex=1;}
        else if(region_id == 3){if(VECTOR_3D<T>::Dot_Product(triangle2.normal,triangle.x2-triangle2.x1) >= 0) convex=1;}
        if(convex){if(triangle.Inside(location) && triangle2.Inside(location)) return true;} // inside both - can use location or point
        else{if(triangle.Inside(location) || triangle2.Inside(location)) return true;}} // inside either - can use location or point
    else{if(triangle.Inside(location)) return true;} // region=3 - face - can use location or point
    return false;}*/
    
    T Close_Point_Phi(int p)
    {VECTOR_3D<T> delta=particles_1.X(p)-close_points(p);
    T phi=delta.Magnitude();if(Close_Point_Inside(p))phi=-phi;
    phi-=thickness;return phi;}

    T Close_Point_Phi_And_Normal(int p,VECTOR_3D<T>& normal)
    {VECTOR_3D<T> delta=particles_1.X(p)-close_points(p);
    T phi=delta.Magnitude();if(Close_Point_Inside(p))phi=-phi;
    normal=delta/phi;phi-=thickness;return phi;}

//#####################################################################
};
}
#endif

