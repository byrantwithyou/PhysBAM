#ifndef __RANDOM_PLACEMENT_2D__
#define __RANDOM_PLACEMENT_2D__

#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>

namespace PhysBAM
{

template<class T>
class RANDOM_PLACEMENT_2D
{
public:
    RANDOM_PLACEMENT_2D() : fixed_scale(1), max_orientation_angle(pi), 
                         min_speed(0), max_speed(0), max_angular_speed(0)
     { }

    virtual T Random_Scale(RANDOM_NUMBERS &random_numbers)
    {
        if (fixed_scale > 0) return fixed_scale;
        else return 0.25 * random_numbers.Get_Uniform_Number((int)1,(int)8);
    }

    void Random_Placement(RANDOM_NUMBERS &random_numbers, RIGID_BODY<TV> &body)
    { 
        Assign_Position(random_numbers, body);
        Assign_Orientation(random_numbers, body);
        Assign_Velocity(random_numbers, body);
        Assign_Angular_Velocity(random_numbers, body);
        body.Update_Angular_Momentum();
    }

    void Set_Fixed_Scale(T scale) 
        { fixed_scale = scale; }
    void Set_Random_Scale() 
        { fixed_scale = 0; }
    void Set_Max_Orientation_Angle(T max_orientation_angle_in)    // in radians
        { max_orientation_angle = max_orientation_angle_in; }
    void Set_Speed_Range(T min_speed_in, T max_speed_in) 
        { min_speed = min_speed_in; max_speed = max_speed_in; }
    void Set_Max_Angular_Speed(T max_angular_speed_in)
        { max_angular_speed = max_angular_speed_in; }
    
protected:
    virtual void Assign_Position(RANDOM_NUMBERS &random_numbers, RIGID_BODY<TV> &body) = 0;

    virtual void Assign_Orientation(RANDOM_NUMBERS &random_numbers, RIGID_BODY<TV> &body)
    {
        body.orientation = random_numbers.Get_Uniform_Number((T)-max_orientation_angle, (T)max_orientation_angle);
    }

    virtual void Assign_Velocity(RANDOM_NUMBERS &random_numbers, RIGID_BODY<TV> &body)
    {
        body.velocity = random_numbers.Get_Direction_2D();
        body.velocity *= random_numbers.Get_Uniform_Number((T)min_speed, (T)max_speed);
    }

    virtual void Assign_Angular_Velocity(RANDOM_NUMBERS &random_numbers, RIGID_BODY<TV> &body)
    {
        body.angular_velocity = random_numbers.Get_Uniform_Number((T)-max_angular_speed, (T)max_angular_speed);
        body.Update_Angular_Momentum();
    }

    T fixed_scale;
    T max_orientation_angle;    // in radians
    T min_speed, max_speed;
    T max_angular_speed;
};

template<class T> 
class RECTANGULAR_RANDOM_PLACEMENT : public RANDOM_PLACEMENT_2D<T>
{
public:
    RECTANGULAR_RANDOM_PLACEMENT(const BOX_2D<T>& box)
        : box(box)
    {}

protected:
    void Assign_Position(RANDOM_NUMBERS &random_numbers, RIGID_BODY<TV> &body) PHYSBAM_OVERRIDE
    {
        body.position.x = random_numbers.Get_Uniform_Number((T)box.xmin,(T)box.xmax);
        body.position.y = random_numbers.Get_Uniform_Number((T)box.ymin,(T)box.ymax);
    }

private:
    BOX_2D<T> box;
};

}

#endif
