//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRON_GEAR
//#####################################################################
#ifndef __HEXAHEDRON_GEAR__
#define __HEXAHEDRON_GEAR__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Forces_And_Torques/BODY_FORCES_3D.h>
namespace PhysBAM{

template<class T,class RW>
class HEXAHEDRON_GEAR:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::frame_rate;using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;

    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int gear[2];
    T cylinder_time;
    VECTOR_3D<T> cylinder_start,cylinder_velocity;
    QUATERNION<T> roller_orientation;
    T roller_speed,roller_friction;
    bool gear_triangle_collisions;
    T gear_repulsion_stiffness;
    int gear_repulsion_sampling;
    
    HEXAHEDRON_GEAR()
        :BASE(0,fluids_parameters.NONE),initial_height((T)2.234),initial_orientation(T(pi/2),VECTOR_3D<T>(0,0,1)),initial_velocity(),initial_angular_velocity(),
        gear_triangle_collisions(false),gear_repulsion_stiffness(0),gear_repulsion_sampling(5)
    {
        solids_parameters.collisions_repulsion_thickness=(T)2e-3;
        solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
        solids_parameters.collision_repulsion_spring_multiplier=100;
        solids_parameters.self_collision_friction_coefficient=1;

        frame_rate=48;last_frame=700;
        solids_parameters.cfl=(T)2;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.perform_self_collision=true;
        solids_parameters.collide_with_interior=true;
        solids_parameters.max_collision_loops=1;

        int example=1;

        if(example==1){
            restart=true;restart_frame=107;
            solids_parameters.collide_with_interior=false;
            solids_parameters.allow_intersections=true;
            output_directory="Hexahedron_Gear/output";}
        else if(example==2){
            restart=true;restart_frame=90;
            solids_parameters.collide_with_interior=false;
            solids_parameters.allow_intersections=true;
            solids_parameters.enforce_tangential_collision_velocity=true;
            output_directory="Hexahedron_Gear/output2";}
        else if(example==3){
            restart=true;restart_frame=107;
            solids_parameters.allow_intersections=true;
            output_directory="Hexahedron_Gear/output3";}

        roller_speed=(T)2;
        roller_friction=1;
        roller_orientation=QUATERNION<T>(0,VECTOR_3D<T>(1,0,0));

        cylinder_time=(T).2;
        cylinder_start=VECTOR_3D<T>(0,2,0);
        cylinder_velocity=VECTOR_3D<T>(0,-5,0);

        //initial_height-=.2;

        gear_repulsion_stiffness=2e4;
    }

    ~HEXAHEDRON_GEAR()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Hexahedralized_Volume();
    HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).hexahedralized_volume;
    HEXAHEDRON_MESH& hexahedron_mesh=hexahedralized_volume.hexahedron_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=hexahedralized_volume.particles;

    int resolution=20;
    GRID<TV> grid(3*resolution+1,3*resolution+1,3*resolution+1,(T)-1,(T)1,(T)-1,(T)1,(T)-1,(T)1);
    hexahedralized_volume.Initialize_Cube_Mesh_And_Particles(grid);
    T radius=(T)(one_third+1e-5);
    for(int h=1;h<=hexahedron_mesh.hexahedrons.m;h++)for(int k=1;k<=8;k++){
        VECTOR_3D<T> r=abs(particles.X(hexahedron_mesh.hexahedrons(k,h)));
        if((r.x>radius)+(r.y>radius)+(r.z>radius)>1)hexahedron_mesh.hexahedrons(k,h)=0;}
    hexahedron_mesh.Delete_Hexahedrons_With_Missing_Nodes();
    hexahedralized_volume.Discard_Valence_Zero_Particles_And_Renumber();
    std::cout<<"particles = "<<particles.array_collection->Size()<<", hexahedrons = "<<hexahedron_mesh.hexahedrons.m<<std::endl;

    hexahedralized_volume.Initialize_Triangulated_Surface();hexahedralized_volume.triangulated_surface->triangle_mesh.Initialize_Segment_Mesh();
    hexahedron_mesh.Initialize_Incident_Hexahedrons();
    hexahedron_mesh.Initialize_Node_On_Boundary();
    int count=0;for(int p=1;p<=particles.array_collection->Size();p++)count+=(*hexahedron_mesh.incident_hexahedrons)(p).m;
    std::cout<<count<<std::endl;
    count=0;for(int p=1;p<=particles.array_collection->Size();p++)if((*hexahedron_mesh.node_on_boundary)(p)) count++;
    std::cout<<count<<std::endl;

    hexahedralized_volume.Set_Density(1000);hexahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    hexahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(hexahedralized_volume.bounding_box->Center());T bottom=hexahedralized_volume.bounding_box->ymin;
    for(int i=1;i<=particles.array_collection->Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}

    RIGID_BODY<TV>* rigid_body;
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(solids_parameters.rigid_body_parameters.list.rigid_bodies.m);
    rigid_body->frame.t=VECTOR_3D<T>(0,-1,0);
    rigid_body->coefficient_of_friction=(T)1;

    if(gear_triangle_collisions) solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/gear_remeshed",.375,true,false,false,false);
    else solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/gear",.375,true,true,false,false);
    gear[0]=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(gear[0]);
    rigid_body->frame.t=VECTOR_3D<T>(-(T).4,1.5,-.75);
    rigid_body->frame.r=roller_orientation;
    rigid_body->angular_velocity=-roller_speed*VECTOR_3D<T>(0,0,1);
    rigid_body->coefficient_of_friction=roller_friction;

    if(gear_triangle_collisions) solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/gear_remeshed",.375,true,false,false,false);
    else solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/gear",.375,true,true,false,false);
    gear[1]=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(gear[1]);
    rigid_body->frame.t=VECTOR_3D<T>((T).4,1.5,-.75);
    rigid_body->frame.r=roller_orientation;
    rigid_body->angular_velocity=roller_speed*VECTOR_3D<T>(0,0,1);
    rigid_body->coefficient_of_friction=roller_friction;

    if(!gear_triangle_collisions) solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    else{
        solids_parameters.collision_body_list.Add_Body(solids_parameters.rigid_body_parameters.list(1)); // only do rigid body collisions with ground
        for(int g=0;g<2;g++){    
            RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list.rigid_bodies(gear[g]);
            TRIANGULATED_SURFACE<T>* gear_surface=TRIANGULATED_SURFACE<T>::Create();
            gear_surface->triangle_mesh.Initialize_Triangle_Mesh(rigid_body.triangulated_surface->triangle_mesh);
            gear_surface->particles.Initialize_Particles(rigid_body.triangulated_surface->particles);
            gear_surface->particles.Store_Velocity();gear_surface->particles.Store_Mass();
            ARRAY<T>::copy(1e20,gear_surface->particles.mass.array);
            solids_parameters.extra_collision_surfaces.Append(gear_surface);}
        Update_Collision_Body_Positions_And_Velocities(restart?restart_frame/frame_rate:0);}
}
//#####################################################################
// Function Initialize_Deformable_Object
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume=*deformable_object.hexahedralized_volume;
    solid_body_collection.Add_Force(Create_Body_Forces<T>(hexahedralized_volume));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(hexahedralized_volume,new DIAGONALIZED_SPLINE_MODEL_3D<T>((T)50000,(T).45,(T)0,(T)10,(T).03),true,(T).15));
}    
//#####################################################################
// Function Update_Collision_Body_Positions_And_Velocities
//#####################################################################
void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE
{
    solids_parameters.rigid_body_parameters.list.rigid_bodies(gear[0])->frame.r=QUATERNION<T>(-roller_speed*time,VECTOR_3D<T>(0,0,1))*roller_orientation;
    solids_parameters.rigid_body_parameters.list.rigid_bodies(gear[1])->frame.r=QUATERNION<T>(roller_speed*time,VECTOR_3D<T>(0,0,1))*roller_orientation;
    if(gear_triangle_collisions)
        for(int g=0;g<2;g++){
            PARTICLES<T,VECTOR_3D<T> >& particles=solids_parameters.extra_collision_surfaces(g+1)->particles;
            RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list.rigid_bodies(gear[g]);
            for(int p=1;p<=particles.array_collection->Size();p++){
                particles.X(p)=rigid_body.World_Space_Point(rigid_body.triangulated_surface->particles.X(p));
                particles.V(p)=rigid_body.Pointwise_Object_Velocity(rigid_body.triangulated_surface->particles.X(p));}}
}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
void Add_External_Forces(ARRAY<VECTOR_3D<T> >& F,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(!gear_repulsion_stiffness) return;
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    SEGMENT_MESH& segment_mesh=*deformable_object.hexahedralized_volume->triangulated_surface->triangle_mesh.segment_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=deformable_object.particles;
    T stiffness=gear_repulsion_stiffness/(gear_repulsion_sampling-1);
    if(time<42/frame_rate) stiffness=0;
    else if(time<60/frame_rate) stiffness*=(time*frame_rate-42)/(60-42);
    for(int g=0;g<2;g++){RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list.rigid_bodies(gear[g]);
        for(int s=1;s<=segment_mesh.segments.m;s++){
            int i,j;segment_mesh.segments.Get(s,i,j);VECTOR_3D<T> Xi=particles.X(i),Xj=particles.X(j);
            for(int a=1;a<gear_repulsion_sampling;a++){
                T t=(T)a/gear_repulsion_sampling,phi;VECTOR_3D<T> X=(1-t)*Xi+t*Xj;
                if(rigid_body.Implicit_Geometry_Lazy_Inside_And_Value(X,phi)){
                    VECTOR_3D<T> normal=X-rigid_body.frame.t;normal.z=0;normal.Normalize();
                    VECTOR_3D<T> force=-stiffness*phi*normal;
                    F(i)+=(1-t)*force;F(j)+=t*force;}}}}
}
//#####################################################################
// Function Remesh_Gears
//#####################################################################
void Remesh_Gears()
{
    TRIANGLE_MESH bad_mesh;PARTICLES<T,VECTOR_3D<T> > bad_particles;TRIANGULATED_SURFACE<T> bad_surface(bad_mesh,bad_particles);
    TRIANGLE_MESH good_mesh;PARTICLES<T,VECTOR_3D<T> > good_particles;TRIANGULATED_SURFACE<T> good_surface(good_mesh,good_particles);
    FILE_UTILITIES::Read_From_File<RW>(data_directory+"/Rigid_Bodies/gear.tri",bad_surface);
    bad_surface.Update_Bounding_Box();BOX_3D<T> box=*bad_surface.bounding_box;
    std::cout<<"box = "<<box<<std::endl;
    ARRAY<VECTOR_3D<T> > curve;
    for(int p=1;p<=bad_particles.array_collection->Size();p++)if(bad_particles.X(p).z == 0 && bad_particles.X(p).Magnitude()>.1)curve.Append(bad_particles.X(p));
    int n=curve.m;std::cout<<"n = "<<n<<std::endl;
    curve.Append(curve(1));
    std::cout<<"curve.m = "<<curve.m<<std::endl;

    for(;;){
        T deviation=1000;int di=0;
        for(int i=1;i<=n;i++){
            int ip=i==1?n:i-1,in=i+1;SEGMENT_3D<T> segment(curve(ip),curve(in));
            T d=segment.Distance_From_Point_To_Segment(curve(i));
            if(deviation>d){deviation=d;di=i;}}
        if(deviation>.001) break;
        //std::cout<<"deviation "<<di<<" = "<<deviation<<std::endl;
        curve.Remove_Index(di);n--;curve(n+1)=curve(1);}

    std::cout<<"n = "<<n<<std::endl;

    T min_distance=1000,max_radius=0;
    for(int i=1;i<=n;i++){
        min_distance=min(min_distance,(curve(i)-curve(i+1)).Magnitude());
        max_radius=max(max_radius,curve(i).Magnitude());}
    std::cout<<"min distance = "<<min_distance<<", max radius = "<<max_radius<<std::endl;
    
    T aspect=10;
    int m_body=(int)(4/min_distance/aspect);
    int m_end=(int)(max_radius/min_distance/aspect);
    std::cout<<"m_body = "<<m_body<<", m_end = "<<m_end<<std::endl;
    good_particles.Increase_Array_Size(4*n*m_body+2*n*m_end);

    std::cout<<"correct number of particles = "<<n*(m_body+1)+2*n*(m_end-2)+2<<std::endl;
    std::cout<<"end particles = "<<2*n*(m_end-1)+2<<std::endl;
    std::cout<<"body particles = "<<n*(m_body+1)<<std::endl;

    good_mesh.Initialize_Circle_Mesh(m_end,n);
    int circle_triangles=good_mesh.triangles.m,circle_particles=m_end*n;
    for(int t=1;t<=circle_triangles;t++){
        int i,j,k;good_mesh.triangles.Get(t,i,j,k);
        good_mesh.triangles.Append(i+circle_particles,k+circle_particles,j+circle_particles);}
    for(int j=0;j<m_end;j++)for(int i=1;i<=n;i++)good_particles.X(good_particles.array_collection->Add_Element())=(T)j/(m_end-1)*curve(i);
    for(int j=0;j<m_end;j++)for(int i=1;i<=n;i++)good_particles.X(good_particles.array_collection->Add_Element())=(T)j/(m_end-1)*curve(i)+VECTOR_3D<T>(0,0,4);

    for(int i=1;i<=n;i++)for(int j=0;j<m_body;j++){
        int a=good_particles.array_collection->Add_Element(),b=good_particles.array_collection->Add_Element(),c=good_particles.array_collection->Add_Element(),d=good_particles.array_collection->Add_Element();
        good_particles.X(a)=curve(i)+VECTOR_3D<T>(0,0,(T)4*j/m_body);
        good_particles.X(b)=curve(i)+VECTOR_3D<T>(0,0,(T)4*(j+1)/m_body);
        good_particles.X(c)=curve(i+1)+VECTOR_3D<T>(0,0,(T)4*j/m_body);
        good_particles.X(d)=curve(i+1)+VECTOR_3D<T>(0,0,(T)4*(j+1)/m_body);
        good_mesh.triangles.Append(a,c,b);good_mesh.triangles.Append(b,d,c);}

    good_mesh.number_nodes=good_particles.array_collection->Size();
    good_surface.Close_Surface(true,1e-3,false);
    good_surface.Remove_Degenerate_Triangles();
    good_mesh.Make_Orientations_Consistent();
    good_surface.Print_Statistics(LOG::cout);

    FILE_UTILITIES::Write_To_File<RW>("Hexahedron_Gear/gear_fixed.tri",good_surface);
    exit(0);
}
//#####################################################################
};
}
#endif
