//#####################################################################
// Copyright 2007, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANDOM_PLACEMENT
//#####################################################################
#ifndef __RANDOM_PLACEMENT__
#define __RANDOM_PLACEMENT__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
namespace PhysBAM{
template<class T>
void Generate_Orientation(RANDOM_NUMBERS<T>& random_numbers,T max_orientation_angle,ROTATION<VECTOR<T,3> >& orientation)
{
    typedef VECTOR<T,3> TV;
    T s=random_numbers.Get_Uniform_Number((T)0,(T)max_orientation_angle);
    TV dir=random_numbers.template Get_Direction<TV>();
    orientation=ROTATION<TV>(s,dir);
}

template<class T>
void Generate_Orientation(RANDOM_NUMBERS<T>& random_numbers,T max_orientation_angle,ROTATION<VECTOR<T,2> >& orientation)
{
    T s=random_numbers.Get_Uniform_Number((T)0,(T)max_orientation_angle);
    orientation=ROTATION<VECTOR<T,2> >::From_Angle(s);
}

template<class T>
void Generate_Angular_Velocity(RANDOM_NUMBERS<T>& random_numbers,VECTOR<T,3>& spin)
{
    spin=random_numbers.template Get_Direction<VECTOR<T,3> >();
}

template<class T>
void Generate_Angular_Velocity(RANDOM_NUMBERS<T>& random_numbers,VECTOR<T,1>& spin)
{
    spin=VECTOR<T,1>((T)1);
}

template<class TV>
class RANDOM_PLACEMENT
{
public:
    typedef typename TV::SCALAR T;

    RANDOM_PLACEMENT()
        :fixed_scale(1),max_orientation_angle((T)pi),min_speed(0),max_speed(0),min_angular_speed(0),max_angular_speed(0)
    {}

    virtual ~RANDOM_PLACEMENT()
    {}

    virtual typename TV::SCALAR Random_Scale(RANDOM_NUMBERS<T>& random_numbers)
    {if(fixed_scale>0) return fixed_scale;else return (T)0.25*random_numbers.Get_Uniform_Number((T)1,(T)8);}

    void Random_Placement(RANDOM_NUMBERS<T>& random_numbers,RIGID_BODY<TV>& body)
    {Assign_Position(random_numbers,body);Assign_Orientation(random_numbers,body);Assign_Velocity(random_numbers,body);Assign_Angular_Velocity(random_numbers,body);
    body.Update_Angular_Momentum();}

    void Set_Fixed_Scale(T scale) 
    {fixed_scale=scale;}

    void Set_Random_Scale() 
    {fixed_scale=0;}

    void Set_Max_Orientation_Angle(T max_orientation_angle_in) // in radians
    {max_orientation_angle=max_orientation_angle_in;}

    void Set_Speed_Range(T min_speed_in,T max_speed_in) 
    {min_speed=min_speed_in;max_speed=max_speed_in;}

    void Set_Angular_Speed_Range(T min_angular_speed_in,T max_angular_speed_in)
    {min_angular_speed=min_angular_speed_in;max_angular_speed=max_angular_speed_in;}

protected:
    virtual void Assign_Position(RANDOM_NUMBERS<T>& random_numbers,RIGID_BODY<TV>& body)=0;

    virtual void Assign_Orientation(RANDOM_NUMBERS<T>& random_numbers,RIGID_BODY<TV>& body)
    {Generate_Orientation<T>(random_numbers,max_orientation_angle,body.Rotation());}

    virtual void Assign_Velocity(RANDOM_NUMBERS<T>& random_numbers,RIGID_BODY<TV>& body)
    {body.Twist().linear=random_numbers.template Get_Direction<TV>();
    body.Twist().linear*=random_numbers.Get_Uniform_Number((T)min_speed,(T)max_speed);}

    virtual void Assign_Angular_Velocity(RANDOM_NUMBERS<T>& random_numbers,RIGID_BODY<TV>& body)
    {Generate_Angular_Velocity(random_numbers,body.Twist().angular);
    body.Twist().angular*=random_numbers.Get_Uniform_Number((T)min_angular_speed,(T)max_angular_speed);
    body.Update_Angular_Momentum();}

    T fixed_scale;
    T max_orientation_angle; // in radians
    T min_speed,max_speed;
    T min_angular_speed,max_angular_speed;
};

template<class TV> 
class SPHERICAL_RANDOM_PLACEMENT:public RANDOM_PLACEMENT<TV>
{
    typedef typename TV::SCALAR T;
    using RANDOM_PLACEMENT<T>::min_speed;using RANDOM_PLACEMENT<T>::max_speed;
public:
    SPHERICAL_RANDOM_PLACEMENT(T radius,TV center=TV(0,0,0))
       :boundary_only(false),radius(radius),center(center),velocity_towards_center(false)
    {}

    virtual ~SPHERICAL_RANDOM_PLACEMENT()
    {}

    void Set_Boundary_Only(bool boundary_only_in=true)
    {boundary_only=boundary_only_in;}

    void Set_Velocity_Towards_Center()
    {velocity_towards_center=true;}

protected:
    void Assign_Position(RANDOM_NUMBERS<T>& random_numbers,RIGID_BODY<TV>& body) PHYSBAM_OVERRIDE
    {T r=(boundary_only)?radius:random_numbers.Get_Uniform_Number((T)0,(T)radius);
    TV d=random_numbers.template Get_Direction<TV>();
    body.frame.t=center+r*d;}

    void Assign_Velocity(RANDOM_NUMBERS<T>& random_numbers,RIGID_BODY<TV>& body) PHYSBAM_OVERRIDE
    {if(velocity_towards_center) body.velocity=(center-body.frame.t).Normalized()*random_numbers.Get_Uniform_Number((T)min_speed,(T)max_speed);
    else RANDOM_PLACEMENT<T>::Assign_Velocity(random_numbers,body);}

private:
    bool boundary_only;
    T radius;
    TV center;
    bool velocity_towards_center;
};

template<class TV>
class CYLINDRICAL_RANDOM_PLACEMENT:public RANDOM_PLACEMENT<TV>
{
    typedef typename TV::SCALAR T;
public:
    CYLINDRICAL_RANDOM_PLACEMENT(T radius,T height,TV base=TV(0,0,0))
       :boundary_only(false),radius(radius),height(height),base(base)
    {}

    virtual ~CYLINDRICAL_RANDOM_PLACEMENT()
    {}

    void Set_Boundary_Only(bool boundary_only_in=true)
    {boundary_only=boundary_only_in;}

protected:
    void Assign_Position(RANDOM_NUMBERS<T>& random_numbers,RIGID_BODY<TV>& body) PHYSBAM_OVERRIDE
    {T r=(boundary_only)?radius:random_numbers.Get_Uniform_Number((T)0,(T)radius);
    TV d=random_numbers.template Get_Direction<TV>();
    d.y=0;if(d.Magnitude()>0) d.Normalize();
    T h=random_numbers.Get_Uniform_Number((T)0,(T)height);
    body.X()=base+r*d+TV(0,h,0);}

private:
    bool boundary_only;
    T radius,height;
    TV base;
};

template<class TV>
class RECTANGULAR_RANDOM_PLACEMENT:public RANDOM_PLACEMENT<TV>
{
public:
    typedef typename TV::SCALAR T;

    RECTANGULAR_RANDOM_PLACEMENT(const BOX<TV>& box)
        :box(box)
    {}

    virtual ~RECTANGULAR_RANDOM_PLACEMENT()
    {}

protected:
    void Assign_Position(RANDOM_NUMBERS<T>& random_numbers,RIGID_BODY<TV>& body) PHYSBAM_OVERRIDE
    {body.X().x=random_numbers.Get_Uniform_Number((T)box.min_corner.x,(T)box.max_corner.x);
    body.X().y=random_numbers.Get_Uniform_Number((T)box.min_corner.y,(T)box.max_corner.y);}

private:
    BOX<TV> box;
};

//#####################################################################
// Function Random_Scene_Generator
//#####################################################################
template<class TV>
void Random_Scene_Generator(const std::string& filename,const int number_objects,const int random_seed,RANDOM_PLACEMENT<TV>& random_placement,
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGIDS_STANDARD_TESTS<TV>& tests)
{
    ARRAY<std::string> filenames(number_objects);ARRAYS_COMPUTATIONS::Fill(filenames,filename);
    Random_Scene_Generator(filenames,random_seed,random_placement,rigid_body_collection,tests);
}

template<class TV>
void Random_Scene_Generator(const ARRAY<std::string>& filenames,const int random_seed,RANDOM_PLACEMENT<TV>& random_placement,
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGIDS_STANDARD_TESTS<TV>& tests)
{
    typedef typename TV::SCALAR T;
    
    RANDOM_NUMBERS<T> random_numbers;
    random_numbers.Set_Seed(random_seed);
    
    int i;
    int start_body=rigid_body_collection.rigid_body_particle.array_collection->Size()+1;
    for(i=1;i<=filenames.m;i++){// create all the objects to get their bounding boxes
        T scale=random_placement.Random_Scale(random_numbers);
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body(filenames(i),scale,1);
        rigid_body.Set_Name(STRING_UTILITIES::string_sprintf("%s %d",filenames(i).c_str(),i));
        random_placement.Random_Placement(random_numbers,rigid_body);}
    
    RIGID_BODY_INTERSECTIONS<TV> intersections(rigid_body_collection);
    LOG::cout<<"Resolving random scene collisions...";
    
    for(int i=start_body;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){
        rigid_body_collection.Rigid_Body(i).Update_Bounding_Box();
        std::cout<<i<<" "<<std::flush;
        bool found_collision=true;int counter=0;
        while(found_collision){
            found_collision=false;
            for(int j(1);j<i;j++) if(rigid_body_collection.Is_Active(j)){
                while(intersections.Bounding_Boxes_Intersect(i,j)){
                    random_placement.Random_Placement(random_numbers,rigid_body_collection.Rigid_Body(i));
                    rigid_body_collection.Rigid_Body(i).Update_Bounding_Box();
                    found_collision=true;}
                if(found_collision) break;}
            counter++;
            if(counter%10000==0) std::cerr<<"Warning: problems generating random scene"<<std::endl;}}
    LOG::cout<<std::endl;
}
}
#endif
