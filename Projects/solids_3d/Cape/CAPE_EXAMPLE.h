//#####################################################################
// Copyright 2002, 2003, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CAPE_EXAMPLE
//##################################################################### 
#ifndef __CAPE_EXAMPLE__
#define __CAPE_EXAMPLE__

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include "../CLOTH_EXAMPLE.h"
#include "CAPE_COLLISIONS.h"
namespace PhysBAM{

template<class T>
class CAPE_EXAMPLE:public CLOTH_EXAMPLE<T>
{
public:
    int number_side_panels;
    T aspect_ratio,side_length;
    T ball_1_radius;
    VECTOR_3D<T> ball_1_start,ball_1_velocity;
    //int cloth_particles,cloth_triangles;
    TRIANGLE_MESH* full_triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >* full_particles;
    TETRAHEDRON_MESH buddha_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> > buddha_particles;
    TETRAHEDRALIZED_VOLUME<T> buddha_volume;
    int buddha_frame;
    char* buddha_directory;
    char* public_data_directory;
    char *cape_file;
    ARRAY<VECTOR_3D<T> > barycentric_coordinates;
    int collision_triangles;
    T cloth_height;
    int example;
    
    T buddha_fractional_time;
    
    int resolution;
    T beginning_of_time;
    T time_scale;
    
    T length_scale;
    T cape_thickness;
    T cape_tolerance;
    T cape_friction;
    CAPE_COLLISIONS<T>* cape_collisions;
    ARRAY<int> attachments;
    
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> > visual_particles;
    TRIANGLE_MESH visual_mesh;
    TRIANGULATED_SURFACE<T> visual_surface;
    ARRAY<int,VECTOR<int,1> > visual_mapping;
    
    TRIANGULATED_SURFACE<T>* collision_surface;
    int max_collision_loops;

    CAPE_EXAMPLE()
        :number_side_panels(40),aspect_ratio((T)1),side_length(1),ball_1_radius(.25),buddha_volume(buddha_mesh,buddha_particles),
        buddha_fractional_time(-1),visual_surface(visual_mesh,visual_particles),
        cape_thickness((T)1e-2),cape_tolerance((T)1e-7),cape_friction((T).03),cape_collisions(0)
    {
        resolution=150;
        beginning_of_time=1;
        time_scale=1;
        output_directory="Cape/output_2";
        cape_file="Cape/capes/capeflare_long_wide.tri";
        
        final_time=2;
        frame_rate=120*4;
        cfl_number=(T)1;
        cg_tolerance=(T)1e-3;
        cg_iterations=1000;
        max_strain_per_time_step=(T).1;
        
        use_masses_and_springs=false;use_altitude_springs=false;use_fvm=false;
        
        use_diagonalized_fvm=true;
        use_linear_elasticity=false;
        youngs_modulus=200;
        poissons_ratio=(T).1;
        Rayleigh_coefficient=(T).02;
        bending_stiffness=(T).002;
        bending_damping=(T).002;
        
        example=1;
        
        if(example==1){
            ball_1_start=VECTOR_3D<T>((T).05,1,-1);
            ball_1_velocity=VECTOR_3D<T>(0,0,15);
            buddha_directory="f:/data/buddha_weak";}
        else if(example==2){
            ball_1_start=VECTOR_3D<T>((T).05,1,-1);
            ball_1_velocity=VECTOR_3D<T>(0,0,10);
            buddha_directory="f:/data/buddha_strong";}
            
        cloth_height=(T).987684;
        
        public_data_directory="../../Public_Data";
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=true;
        //allow_intersections=true;
        //allow_intersections_tolerance=1e-8;
        gravity=10;
        
        
//#####################################################################

        restart_step_number=658;  // sane again
        restart_step_number=663;  // gave the ground priority some thickness
        restart_step_number=664;  // doubled it
        
//#####################################################################
    }

    ~CAPE_EXAMPLE()
    {}
    
    void Read_Fractional_Buddha(int frame,T fraction,ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& X,ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V)
    {std::fstream input;char filename[200];
    ARRAY<VECTOR_3D<T> ,VECTOR<int,1> > Q;
    sprintf(filename,"%s/position.%d",buddha_directory,frame);input.open(filename,std::ios::in|std::ios::binary);Q.Read_Objects_Float(input);input.close();
    sprintf(filename,"%s/position.%d",buddha_directory,frame+1);input.open(filename,std::ios::in|std::ios::binary);buddha_particles.X.Read_Objects_Float(input);input.close();
    for(int p=0;p<buddha_particles.array_collection->Size();p++)X(p)=(1-fraction)*Q(p)+fraction*buddha_particles.X(p);
    sprintf(filename,"%s/velocity.%d",buddha_directory,frame);input.open(filename,std::ios::in|std::ios::binary);Q.Read_Objects_Float(input);input.close();
    sprintf(filename,"%s/velocity.%d",buddha_directory,frame+1);input.open(filename,std::ios::in|std::ios::binary);buddha_particles.V.Read_Objects_Float(input);input.close();
    for(int p=0;p<buddha_particles.array_collection->Size();p++)V(p)=(1-fraction)*Q(p)+fraction*buddha_particles.V(p);}
    
    void Read_Fractional_Buddha(T time)
    {std::cout<<"reading fractional "<<time<<"\n";
    T real_time=max((T)(0),time-beginning_of_time);
    if(real_time==buddha_fractional_time)return;buddha_fractional_time=real_time;
    int frame=(int)(real_time*120);T frame_time=(T)frame/120,fraction=(real_time-frame_time)*120;
    Read_Fractional_Buddha(frame,fraction,buddha_particles.X,buddha_particles.V);
    if(time<beginning_of_time)ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::copy(VECTOR_3D<T>(0),buddha_particles.V);
    std::cout<<"done\n";}
    
    bool Corner_Search(ARRAY<ARRAY<int> >& graph,DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles,ARRAY<VECTOR_3D<T> >& edge,ARRAY<bool>& marks,int node,int goal)
    {marks(node)=true;
    std::cout<<"node "<<node<<", degree "<<graph(node).m<<"\n";
    for(int a=0;a<graph(node).m;a++){
        int p=graph(node)(a);
        if(p==goal){std::cout<<"found "<<p<<"\n";edge.Append(particles.X(p));edge.Append(particles.X(node));return true;}
        else if(!marks(p) && Corner_Search(graph,particles,edge,marks,p,goal)){edge.Append(particles.X(node));return true;}}
    std::cout<<"node "<<node<<" failed\n";
    return false;}
    
    void Retriangulate(TRIANGULATED_SURFACE<T>& triangulated_surface)
    {TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& X=particles.X;
    triangle_mesh.Initialize_Neighbor_Nodes();
    triangle_mesh.Initialize_Boundary_Mesh();
    triangle_mesh.boundary_mesh->Initialize_Neighbor_Nodes();
    ARRAY<ARRAY<int> >& neighbor_nodes=*triangle_mesh.neighbor_nodes;
    ARRAY<ARRAY<int> >& boundary_neighbor_nodes=*triangle_mesh.boundary_mesh->neighbor_nodes;
    ARRAY<int> corner;
    corner.Append(6);
    corner.Append(12);
    corner.Append(1);
    corner.Append(5);
    /*for(int p=0;p<particles.array_collection->Size();p++){
        std::cout<<p<<" - "<<particles.X(p)<<"\n";/*
        if(!boundary_neighbor_nodes(p).m)continue;
        if(boundary_neighbor_nodes(p).m==2)corner.Append(p);}
        int a=boundary_neighbor_nodes(p)(1),b=boundary_neighbor_nodes(p)(2);
        if(fabs(VECTOR_3D<T>::Dot_Product(X(a)-X(p),X(b)-X(p)))<(X(a)-X(p)).Magnitude_Squared())corner.Append(p);*///}
    std::cout<<"Found "<<corner.m<<" corners.\n";
    std::cout<<corner(1)<<" "<<corner(2)<<" "<<corner(3)<<" "<<corner(4)<<"\n";
    for(int a=0;a<4;a++)std::cout<<neighbor_nodes(corner(a)).m<<" ";std::cout<<"\n";
    if(corner.m!=4){assert(false);exit(1);}
    ARRAY<VECTOR_3D<T> > cx(1,4);for(int a=0;a<4;a++)cx(a)=particles.X(corner(a));
    T d2=(cx(1)-cx(2)).Magnitude(),d3=(cx(1)-cx(3)).Magnitude(),d4=(cx(1)-cx(4)).Magnitude();
    T min_d=min(d2,d3,d4);
    if(d3==min_d){exchange(corner(2),corner(3));exchange(cx(2),cx(3));}
    else if(d4==min_d){exchange(corner(2),corner(4));exchange(cx(2),cx(4));}
    if((cx(1)-cx(3)).Magnitude()>(cx(1)-cx(4)).Magnitude()){exchange(corner(3),corner(4));exchange(cx(3),cx(4));}
    T dist[4][4];for(int a=0;a<4;a++)for(int b=0;b<4;b++)dist[a][b]=(cx(a+1)-cx(b+1)).Magnitude();
    ARRAY<bool> marks(1,particles.array_collection->Size());
    ARRAY<VECTOR_3D<T> > edge1,edge2;
    for(int a=0;a<4;a++)marks(corner(a))=true;
    bool result1=Corner_Search(boundary_neighbor_nodes,particles,edge1,marks,corner(1),corner(2));assert(result1);
    bool result2=Corner_Search(boundary_neighbor_nodes,particles,edge2,marks,corner(3),corner(4));assert(result2);
    edge1.Compact();edge2.Compact();
    std::cout<<"edge1 "<<edge1.m<<", edge2 "<<edge2.m<<"\n";
    triangle_mesh.Clean_Memory();
    particles.Clean_Memory();particles.Store_Position();
    GRID<TV> edge1_grid(edge1.m,0,1),edge2_grid(edge2.m,0,1);
    LINEAR_INTERPOLATION<T,VECTOR_3D<T> > interpolation;
    int m=resolution,n=2*m;
    triangle_mesh.Initialize_Herring_Bone_Mesh(m,n);
    particles.Increase_Array_Size(m*n);
    for(int j=0;j<n;j++)for(int i=0;i<m;i++){
        VECTOR_3D<T> x1=interpolation.Clamped_To_Array(edge1_grid,edge1.array,(T)(i-1)/m);
        VECTOR_3D<T> x2=interpolation.Clamped_To_Array(edge2_grid,edge2.array,(T)(i-1)/m);
        particles.X(particles.array_collection->Add_Element())=x1+T(j-1)/n*(x2-x1);}
    std::fstream output("Cape/tmp.tri",std::ios::binary|std::ios::out);triangulated_surface.Write_Float(output);output.close();}
    
    void Extract_Component(TRIANGULATED_SURFACE<T>& triangulated_surface,int seed,TRIANGULATED_SURFACE<T>& component)
    {TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    if(!triangle_mesh.neighbor_nodes)triangle_mesh.Initialize_Neighbor_Nodes();
    int count=0;ARRAY<int> mark(1,particles.array_collection->Size());ARRAY<int> stack;
    stack.Append(seed);mark(seed)=++count;
    while(stack.m){
        int p;stack.Pop(p);
        for(int a=0;a<(*triangle_mesh.neighbor_nodes)(p).m;a++){
            int q=(*triangle_mesh.neighbor_nodes)(p)(a);
            if(!mark(q)){stack.Append(q);mark(q)=++count;}}}
    int pmin=particles.array_collection->Size(),pmax=1;
    for(int p=0;p<particles.array_collection->Size();p++)if(mark(p)){pmin=min(pmin,mark(p));pmax=max(pmax,mark(p));}
    if(pmax-pmin+1!=count){std::cout<<"Lazy\n";exit(1);}
    component.Clean_Memory();component.particles.Clean_Memory();component.triangle_mesh.Clean_Memory();
    component.particles.Increase_Array_Size(count);
    component.particles.Store_Position();
    for(int p=pmin;p<=pmax;p++){component.particles.X(component.particles.array_collection->Add_Element())=particles.X(p);}
    component.triangle_mesh.triangles.Append_Elements(triangle_mesh.triangles);
    component.triangle_mesh.number_nodes=count;
    for(int t=0;t<triangle_mesh.triangles.m;t++)for(int a=0;a<3;a++){
        int p=component.triangle_mesh.triangles(a,t);
        component.triangle_mesh.triangles(a,t)=p>=pmin&&p<=pmax?p-pmin+1:0;}
    component.triangle_mesh.Delete_Triangles_With_Missing_Nodes();}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data(TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    full_triangle_mesh=&triangle_mesh;full_particles=&particles;
    
 //   if(!restart_step_number){
        std::ifstream input(cape_file,std::ios::binary);assert(input.is_open());triangulated_surface.Read_Float(input);input.close();
        Retriangulate(triangulated_surface);
        for(int p=0;p<particles.array_collection->Size();p++)particles.X(p).y+=cloth_height;
/*    else{
        std::ifstream input;char filename[200];
        sprintf(filename,"%s/triangle_mesh",output_directory);input.open(filename,std::ios::binary);assert(input);visual_mesh.Read(input);input.close();
        sprintf(filename,"%s/particle_class_state",output_directory);input.open(filename,std::ios::binary);assert(input);visual_particles.Read_State(input);input.close();
        sprintf(filename,"%s/position.0",output_directory);input.open(filename,std::ios::binary);assert(input);visual_particles.X.Read_Objects_Float(input);input.close();
        Extract_Component(visual_surface,1,triangulated_surface);
        visual_particles.Clean_Memory();visual_mesh.Clean_Memory();visual_surface.Clean_Memory();} */
        
    particles.Delete_Velocity_And_Acceleration();particles.Delete_Mass();
    particles.Update_Position_And_Velocity();particles.Store_Mass();
    buddha_particles.Update_Position_And_Velocity();
    triangulated_surface.Set_Density(1);triangulated_surface.Set_Mass_Of_Particles(true);
    
    char filename[200];sprintf(filename,"%s/tetrahedron_mesh",buddha_directory);
    {std::ifstream input;input.open(filename,std::ios::binary);assert(input);buddha_mesh.Read(input);input.close();}
    buddha_particles.Increase_Array_Size(buddha_mesh.number_nodes);
    for(int p=0;p<buddha_mesh.number_nodes;p++)buddha_particles.array_collection->Add_Element();
    buddha_mesh.number_nodes=buddha_particles.array_collection->Size();
    buddha_volume.Initialize_Triangulated_Surface();

    triangle_mesh.triangles.Compact();
    std::cout<<"allocating cape_collisions\n";
    cape_collisions=new CAPE_COLLISIONS<T>(triangulated_surface,*buddha_volume.triangulated_surface);
    length_scale=triangulated_surface.Minimum_Edge_Length();
    cape_collisions->thickness=(T).8*length_scale;
    cape_collisions->tolerance=cape_tolerance;
    cape_collisions->coefficient_of_friction=cape_friction;
//    cape_collisions->flypaper_thickness=cape_collisions->thickness*3;
    cape_collisions->ground_thickness=2*cape_collisions->thickness;
    
    Read_Fractional_Buddha(0);
    barycentric_coordinates.Resize(1,particles.array_collection->Size());
    attachments.Resize(1,particles.array_collection->Size());
    std::cout<<"Turning off collisions for points along the seam...";
    cape_collisions->bound=8*length_scale;
    cape_collisions->thickness*=2;
    cape_collisions->Disable_Currently_Colliding_Particles();
    cape_collisions->thickness/=2;
    std::cout<<"done\n";
    BOX_3D<T> attachment_box(100,-100,100,-100,100,-100);ARRAY<int> attachment_tets;
    std::cout<<"Attaching";
    for(int p=0;p<particles.array_collection->Size();p++)if((*cape_collisions->disabled)(p))
        attachment_box.Enlarge_To_Include_Point(particles.X(p));
    attachment_box.Scale_About_Center((T)1.1);
    for(int t=0;t<buddha_mesh.tetrahedrons.m;t++){
        int i,j,k,l;buddha_mesh.tetrahedrons.Get(t,i,j,k,l);
        BOX_3D<T> box(buddha_particles.X(i));box.Enlarge_To_Include_Point(buddha_particles.X(j));
        box.Enlarge_To_Include_Point(buddha_particles.X(k));box.Enlarge_To_Include_Point(buddha_particles.X(l));
        if(box.Intersection(attachment_box))attachment_tets.Append(t);}
    for(int p=0;p<particles.array_collection->Size();p++)if((*cape_collisions->disabled)(p)){
        VECTOR_3D<T> X=particles.X(p),w;
        attachments(p)=-1;
        std::cout<<".";
        for(int a=0;a<attachment_tets.m;a++){
            int t=attachment_tets(a);
            int i,j,k,l;buddha_mesh.tetrahedrons.Get(t,i,j,k,l);
            w=TETRAHEDRON<T>::Barycentric_Coordinates(X,buddha_particles.X(i),buddha_particles.X(j),buddha_particles.X(k),buddha_particles.X(l));
            if(w.x>=0 && w.y>=0 && w.z>=0 && w.x+w.y+w.z<=1){attachments(p)=t;barycentric_coordinates(p)=w;break;}}
        if(attachments(p)==-1){
            if(buddha_volume.triangulated_surface->Inside(particles.X(p),(T)1e-7)){
                std::cout<<"\nbroken\n";exit(1);}}}
    std::cout<<"done\n";
    
    bool neighbor_nodes_defined=triangle_mesh.neighbor_nodes!=0;if(!neighbor_nodes_defined)triangle_mesh.Initialize_Neighbor_Nodes();
    ARRAY<int> disabled1(1,particles.array_collection->Size());
    ARRAY<int> disabled2(1,particles.array_collection->Size());
    ARRAY<int>::copy(attachments,disabled1);
    for(int p=0;p<particles.array_collection->Size();p++)if(disabled1(p)>0)
        for(int a=0;a<(*triangle_mesh.neighbor_nodes)(p).m;a++){
            disabled2((*triangle_mesh.neighbor_nodes)(p)(a))=1;}
    ARRAY<int>::copy(0,disabled1);
    for(int p=0;p<particles.array_collection->Size();p++)if(disabled2(p)>0)
        for(int a=0;a<(*triangle_mesh.neighbor_nodes)(p).m;a++){
            disabled1((*triangle_mesh.neighbor_nodes)(p)(a))=1;}
    if(!neighbor_nodes_defined){delete triangle_mesh.neighbor_nodes;triangle_mesh.neighbor_nodes=0;}
    
    collision_surface->triangle_mesh.triangles.Append_Elements(triangle_mesh.triangles);
    for(int t=0;t<collision_surface->triangle_mesh.triangles.m;t++){
        int i,j,k;triangle_mesh.triangles.Get(t,i,j,k);
        if(disabled1(i)>0 || disabled1(j)>0 || disabled1(k)>0)collision_surface->triangle_mesh.triangles(1,t)=0;}
    collision_surface->triangle_mesh.Delete_Triangles_With_Missing_Nodes();
    std::cout<<"collision triangles: "<<collision_surface->triangle_mesh.triangles.m<<"\n";
    collision_surface->triangle_mesh.number_nodes=particles.array_collection->Size();
    
    for(int t=0;t<triangle_mesh.triangles.m;t++){
        int i,j,k;triangle_mesh.triangles.Get(t,i,j,k);
        if(!attachments(i) && !attachments(j) && !attachments(k))exchange(triangle_mesh.triangles(1,t),triangle_mesh.triangles(2,t));}

    visual_mesh.triangles.Append_Elements(triangle_mesh.triangles);
    visual_mesh.triangles.Append_Elements(buddha_volume.triangulated_surface->triangle_mesh.triangles);
    for(int t=0;t<buddha_volume.triangulated_surface->triangle_mesh.triangles.m;t++)for(int a=0;a<3;a++)
        visual_mesh.triangles(a,t+triangle_mesh.triangles.m)+=particles.array_collection->Size();
    visual_mesh.number_nodes=particles.array_collection->Size()+buddha_particles.array_collection->Size();
    visual_particles.Initialize_Particles(particles);
    visual_particles.Increase_Array_Size(buddha_particles.array_collection->Size());
    for(int p=0;p<buddha_particles.array_collection->Size();p++)visual_particles.array_collection->Add_Element();
    visual_surface.Discard_Valence_Zero_Particles_And_Renumber(visual_mapping);
    
    std::cout<<"triangles = "<<triangle_mesh.triangles.m<<"\n";
}
//#####################################################################
// Function Read_Data_Files
//#####################################################################
void Read_Data_Files(TRIANGULATED_SURFACE<T>& triangulated_surface,T& time,const int frame)
{
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    std::ifstream input;char filename[256];
    visual_particles.Read_Deformable_Dynamic_Variables_Float(output_directory,frame);
    std::cout<<"reading "<<particles.array_collection->Size()<<"\n";
    ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::copy_up_to(visual_particles.X,particles.X,particles.array_collection->Size());
    ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::copy_up_to(visual_particles.V,particles.V,particles.array_collection->Size());
    sprintf(filename,"%s/time.%d",output_directory,frame);
    if(FILE_UTILITIES::File_Exists(filename)){
        std::ifstream input(filename,std::ios::binary);Read_Binary_Float(input,time);input.close();
        initial_time=time-restart_step_number/frame_rate;}
}
//#####################################################################
// Function Write_Data_Files
//#####################################################################
void Write_Data_Files(const TRIANGULATED_SURFACE<T>& triangulated_surface,const ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const T time,const int frame)
{      
    bool visual=true;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    std::fstream output;char filename[256];
    if(frame == 0){
        sprintf(filename,"%s/triangle_mesh",output_directory);output.open(filename,std::ios::out|std::ios::binary);
        assert(output.is_open());
        if(visual)visual_mesh.Write(output);else triangulated_surface.triangle_mesh.Write(output);
        output.close();
        if(visual)visual_particles.Write_Deformable_Static_Variables_Float(output_directory);
        else particles.Write_Deformable_Static_Variables_Float(output_directory);}
    if(visual){
        ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::copy(VECTOR_3D<T>(0),visual_particles.X);
        ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::copy_up_to(particles.X,visual_particles.X,particles.array_collection->Size());
        ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >::copy_up_to(particles.V,visual_particles.V,particles.array_collection->Size());
        Read_Fractional_Buddha(time);
        for(int p=0;p<buddha_particles.array_collection->Size();p++)visual_particles.X(visual_mapping(p+particles.array_collection->Size()))=buddha_particles.X(p);
        visual_particles.Write_Deformable_Dynamic_Variables_Float(output_directory,frame);}
    else particles.Write_Deformable_Dynamic_Variables_Float(output_directory,frame);

    if(frame == 0){triangulated_surface_list.Write_Float(output_directory);implicit_surface_list.Write_Float(output_directory);}
    RIGID_BODY<TV>::Write_Float(output_directory,rigid_bodies,frame);
    
    sprintf(filename,"%s/LAST_VALID_FRAME.txt",output_directory);output.open(filename,std::ios::out);output << frame;output.close();
    if(time){sprintf(filename,"%s/time.%d",output_directory,frame);output.open(filename,std::ios::out|std::ios::binary);Write_Binary(output,(float)time);output.close();}
}
//#####################################################################
// Function Initialize_Diagonalized_Finite_Volume_Model
//#####################################################################
void Initialize_Diagonalized_Finite_Volume_Model(DEFORMABLE_OBJECT<T,VECTOR_3D<T> >& deformable_object,TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    
    strain=new STRAIN_MEASURE_S3D<T>(triangulated_surface);
    diagonalized_constitutive_model=new DIAGONALIZED_LINEAR_FVM_2D<T>(youngs_modulus,poissons_ratio,Rayleigh_coefficient);
    
    diagonalized_fvm=new DIAGONALIZED_FINITE_VOLUME_S3D<T>(*strain,*diagonalized_constitutive_model);
    solid_body_collection.solids_forces.Append(diagonalized_fvm);
    diagonalized_fvm->Limit_Time_Step_By_Strain_Rate(true,max_strain_per_time_step);
    diagonalized_fvm->Use_Rest_State_For_Strain_Rate();

    // dummy edge springs for collisions
    delete triangulated_surface.triangle_mesh.segment_mesh;triangulated_surface.triangle_mesh.Initialize_Segment_Mesh();
    ls=new LINEAR_SPRINGS<T,VECTOR_3D<T> >(*triangulated_surface.triangle_mesh.segment_mesh,particles);
    ls->Set_Restlength_From_Particles();ls->Set_Stiffness((T)150000);
}
//#####################################################################
// Function Initialize_Bending_Elements
//#####################################################################
void Initialize_Bending_Elements(DEFORMABLE_OBJECT<T,VECTOR_3D<T> >& deformable_object,TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    bend=new TRIANGLE_BENDING_ELEMENTS<T>(triangulated_surface.particles);
    solid_body_collection.solids_forces.Append(bend);
    bend->Set_Quadruples_From_Triangle_Mesh(triangulated_surface.triangle_mesh);
    bend->Set_Sine_Half_Rest_Angle(0);
    bend->Set_Stiffness(bending_stiffness);
    bend->Set_Damping(bending_damping);
    //if(use_bending_plasticity)bend->Enable_Plasticity(bending_plastic_yield,bending_plastic_hardening);
    //if(bending_cutoff_ratio)bend->Set_Normal_Cutoff_From_Triangulated_Surface(triangulated_surface,bending_cutoff_ratio);
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
void Initialize_Rigid_Bodies(RIGID_BODY_LIST_3D<T>& solids_parameters.rigid_body_parameters.list)
{
    int index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T)1);

    index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/sphere",ball_1_radius);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->position=ball_1_start;
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->velocity=ball_1_velocity;
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T)0);
}
//#####################################################################
// Function Update_Rigid_Body_Positions
//#####################################################################
void Update_Rigid_Body_Positions(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const T time)
{
    T real_time=max((T)0,time-beginning_of_time);
    rigid_bodies(2)->position=ball_1_start+real_time*ball_1_velocity;
}
//#####################################################################
// Function Update_Rigid_Body_Velocities
//#####################################################################
void Update_Rigid_Body_Velocities(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const T time)
{}
//#####################################################################
// Function Cape_Set_Positions
//#####################################################################
void Cape_Set_Positions(const T time)
{
    //Read_Fractional_Buddha(time);
    for(int p=0;p<full_particles->number;p++){
        int t=attachments(p);if(t<=0)continue;
        int i,j,k,l;buddha_mesh.tetrahedrons.Get(t,i,j,k,l);
        VECTOR_3D<T> w=barycentric_coordinates(p);
        full_particles->X(p)=w.x*buddha_particles.X(i)+w.y*buddha_particles.X(j)+w.z*buddha_particles.X(k)+(1-w.x-w.y-w.z)*buddha_particles.X(l);}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
// for external forces and velocities
void Set_External_Velocities(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    Read_Fractional_Buddha(time);
    for(int p=0;p<full_particles->number;p++){
        int t=attachments(p);if(t<=0)continue;
        int i,j,k,l;buddha_mesh.tetrahedrons.Get(t,i,j,k,l);
        VECTOR_3D<T> w=barycentric_coordinates(p);
        V(p)=w.x*buddha_particles.V(i)+w.y*buddha_particles.V(j)+w.z*buddha_particles.V(k)+(1-w.x-w.y-w.z)*buddha_particles.V(l);}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
// for external forces and velocities
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int p=0;p<full_particles->number;p++)if(attachments(p)>0)V(p)=VECTOR_3D<T>(0);
}
//#####################################################################
// Function Initialize_Triangulated_Surface_Collisions
//#####################################################################
void Initialize_Triangulated_Surface_Collisions(TRIANGULATED_SURFACE<T>& triangulated_surface,TRIANGLE_COLLISIONS<T>*& triangle_collisions)
{
    // set up repulsion springs for collisions
    if(!collisions_repulsion_spring_constant_over_mass_times_length){
        T average_restlength=ARRAY<T,VECTOR<int,1> >::sum(ls->restlength)/ls->restlength.m;
        T mass_node=ARRAY<T,VECTOR<int,1> >::sum(triangulated_surface.particles.mass)/triangulated_surface.particles.array_collection->Size();
        collisions_repulsion_spring_constant_over_mass_times_length=2*ls->constant_youngs_modulus/(mass_node*average_restlength);
        collisions_repulsion_spring_constant_over_mass_times_length/=100;  // FIDDLING
        }
    
    // triangle collisions and friction - note that the initial bounding boxes are initialized
    triangle_collisions=new TRIANGLE_COLLISIONS<T>(triangulated_surface);
    triangle_collisions->Allow_Intersections(allow_intersections);
    triangle_collisions->Set_Allow_Intersections_Tolerance(allow_intersections_tolerance);
    triangle_collisions->Set_Small_Number(collisions_small_number);
    triangle_collisions->Set_Repulsion_Thickness(collisions_repulsion_thickness);
    //triangle_collisions->Set_Repulsion_Thickness(.2*length_scale);
    triangle_collisions->Clamp_Repulsion_Thickness_With_Triangulated_Surface((T).1);//collisions_repulsion_clamp_fraction); // new
    triangle_collisions->Set_Collision_Thickness(collisions_collision_thickness);
    triangle_collisions->Set_Up_Repulsion_Spring(collisions_repulsion_spring_constant_over_mass_times_length);
    if(collisions_output_repulsion_results) triangle_collisions->Output_Repulsion_Results();
    if(collisions_repulsion_thickness) triangle_collisions->Output_Collision_Results();
    if(collisions_output_number_checked) triangle_collisions->Output_Number_Checked();
    triangle_collisions->Set_Friction_Coefficient(self_collision_friction_coefficient);
    triangle_collisions->Set_Attempts_For_Nonrigid_Collisions(collisions_nonrigid_collision_attempts);
}
//#####################################################################
};
}
#endif

