#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <iostream>
#include <string>
#include "../../Public_Library/Rigid_Bodies/RIGID_BODY_LIST_2D.h"
#include "../../Public_Library/Rigid_Bodies/RIGID_BODY_LIST_3D.h"

using namespace PhysBAM;
using namespace FILE_UTILITIES;

PARSE_ARGS parse_args;
int body_index = 0; // 0 is for all
int num_frames = 0;
int start_frame = 0;
bool get_max = false;
double frame_rate = 30;
bool plot_position = false, plot_velocity = false, plot_angular_velocity = false;
std::string input_directory="";

bool frame_valid(const std::string &input_directory, int frame)
{
    int last_frame;
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/last_frame",last_frame);
    return (frame <= last_frame);
}

template<class T> void print_components(const VECTOR<T,3> &vector, const VECTOR<T,3> &tangent, const VECTOR<T,3> &normal, const VECTOR<T,3> &binormal)
{
    std::cout << "t = " << VECTOR<T,3>::Dot_Product(vector, tangent);
    std::cout << ", n = " << VECTOR<T,3>::Dot_Product(vector, normal);
    std::cout << ", b = " << VECTOR<T,3>::Dot_Product(vector, binormal);
}

template<class T> void print_components(const VECTOR<T,2> &vector, const VECTOR<T,2> &tangent, const VECTOR<T,2> &normal)
{
    std::cout << "t = " << VECTOR<T,2>::Dot_Product(vector, tangent);
    std::cout << ", n = " << VECTOR<T,2>::Dot_Product(vector, normal);
}

template<class T> void Do_It_2D()
{
    bool use_frame = false;
    T incline_angle = 0;
    VECTOR<T,2> tangent, normal;

    if (parse_args.Is_Value_Set("-f"))
    {
#if 0
        use_frame = true;
        VECTOR<T,3> rotation_vector(parse_args.Get_Vector_3D_Value("-f"));
        QUATERNION<T> orientation = QUATERNION<T>::From_Rotation_Vector(rotation_vector);
        orientation.Get_Rotated_Frame(tangent, normal, binormal);
#endif
    }
    if (parse_args.Is_Value_Set("-incline"))
    {
#if 0
        use_frame = true;
        incline_angle = parse_args.Get_Double_Value("-incline");
        QUATERNION<T> orientation(incline_angle * pi / 180.0, VECTOR<T,3>(0,0,1));
        orientation.Get_Rotated_Frame(tangent, normal, binormal);
#endif
    }

    RIGID_BODY_LIST_2D<T> rigid_body_list;

    for (int frame = start_frame; ; frame++)
    {
        if (num_frames > 0 && frame > num_frames) break;
        if (!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/rigid_bodies.%d",input_directory.c_str(),frame))) break;
        if (!frame_valid(input_directory, frame)) break;

        rigid_body_list.template Read<T>(input_directory, frame);

        T time=(T)frame/frame_rate;
        std::string time_file=STRING_UTILITIES::string_sprintf("%s/time.%d",input_directory.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(time_file)) FILE_UTILITIES::Read_From_File<T>(time_file,time);

        if (frame == start_frame && body_index != 0 && !(plot_velocity || plot_position || plot_angular_velocity))
            std::cout << "SIMULATION OUTPUT FOR BODY " << body_index << " (" << rigid_body_list(body_index)->name << ")" << std::endl << std::endl;

        T max_speed = -1, max_angular_speed = -1;
        int max_speed_index = 0, max_angular_speed_index = 0;
        if (!plot_position && !plot_velocity && !plot_angular_velocity) std::cout << "time = " << time << " (frame = " << frame << ")" << std::endl;
        for (int i = 1; i <= rigid_body_list.Number_Of_Active_Elements(); i++) if(rigid_body_list.Is_Active(i)){
            if (body_index == 0 || i == body_index){
                if (plot_position) std::cout << time << " " << -VECTOR<T,2>::Dot_Product(rigid_body_list(i)->frame.t, tangent) << std::endl;
                else if (plot_velocity) std::cout << time << " " << -VECTOR<T,2>::Dot_Product(rigid_body_list(i)->velocity, tangent) << std::endl;
                else if (plot_angular_velocity) std::cout << time << " " << rigid_body_list(i)->angular_velocity << std::endl;
                else if (get_max)
                {
                    if (rigid_body_list(i)->velocity.Magnitude() > max_speed)
                    {
                        max_speed = rigid_body_list(i)->velocity.Magnitude();
                        max_speed_index = i;
                    }

                    rigid_body_list(i)->Update_Angular_Velocity();
                    if (fabs(rigid_body_list(i)->angular_velocity) > max_angular_speed)
                    {
                        max_angular_speed = fabs(rigid_body_list(i)->angular_velocity);
                        max_angular_speed_index = i;
                    }
                }
                else
                {
                    if (body_index == 0) std::cout << "body " << i << " (" << rigid_body_list(i)->name << ")" << std::endl;
                    std::cout << "position = " << rigid_body_list(i)->frame.t << std::endl;
                    if (use_frame)
                    {
                        std::cout << "position components: ";
                        print_components(rigid_body_list(i)->frame.t, tangent, normal);
                        std::cout << std::endl;
                    }
                    std::cout << "orientation = " << rigid_body_list(i)->frame.r << std::endl;
                    std::cout << "velocity = " << rigid_body_list(i)->velocity << std::endl;
                    if (use_frame)
                    {
                        std::cout << "velocity components: ";
                        print_components(rigid_body_list(i)->velocity, tangent, normal);
                        std::cout << std::endl;
                    }
                    rigid_body_list(i)->Update_Angular_Velocity();
                    std::cout << "angular_velocity = " << rigid_body_list(i)->angular_velocity << std::endl;
                    std::cout << std::endl;
                }
            }
        }

        if (get_max)
        {
            std::cout << "body " << max_speed_index << " (" << rigid_body_list(max_speed_index)->name << ") has max speed: " << max_speed << std::endl;
            std::cout << "body " << max_angular_speed_index << " (" << rigid_body_list(max_angular_speed_index)->name << ") has max angular speed: " << max_angular_speed << std::endl;
        }
    }

}
template<class T> void Do_It_3D()
{
    bool use_frame = false;
    T incline_angle = 0;
    VECTOR<T,3> tangent, normal, binormal;

    if (parse_args.Is_Value_Set("-f"))
    {
        use_frame = true;
        VECTOR<T,3> rotation_vector(parse_args.Get_Vector_3D_Value("-f"));
        QUATERNION<T> orientation = QUATERNION<T>::From_Rotation_Vector(rotation_vector);
        orientation.Get_Rotated_Frame(tangent, normal, binormal);
    }
    if (parse_args.Is_Value_Set("-incline"))
    {
        use_frame = true;
        incline_angle = parse_args.Get_Double_Value("-incline");
        QUATERNION<T> orientation(incline_angle * pi / 180.0, VECTOR<T,3>(0,0,1));
        orientation.Get_Rotated_Frame(tangent, normal, binormal);
    }

    RIGID_BODY_LIST_3D<T> rigid_body_list;

    for (int frame = start_frame; ; frame++)
    {
        if (num_frames > 0 && frame > num_frames) break;
        if (!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/rigid_bodies.%d",input_directory.c_str(),frame))) break;
        if (!frame_valid(input_directory, frame)) break;

        rigid_body_list.template Read<T>(input_directory, frame);

        T time=(T)frame/frame_rate;
        std::string time_file=STRING_UTILITIES::string_sprintf("%s/time.%d",input_directory.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(time_file)) FILE_UTILITIES::Read_From_File<T>(time_file,time);

        if (frame == start_frame && body_index != 0 && !(plot_velocity || plot_position || plot_angular_velocity))
        {
            std::cout << "SIMULATION OUTPUT FOR BODY " << body_index << " (" << rigid_body_list(body_index)->name << ")" << std::endl << std::endl;
        }

        T max_speed = -1, max_angular_speed = -1;
        int max_speed_index = 0, max_angular_speed_index = 0;
        if (!plot_position && !plot_velocity && !plot_angular_velocity) std::cout << "time = " << time << " (frame = " << frame << ")" << std::endl;
        for (int i = 1; i <= rigid_body_list.Number_Of_Active_Elements(); i++) if(rigid_body_list.Is_Active(i)){
            if (body_index == 0 || i == body_index){
                if (plot_position) std::cout << time << " " << -VECTOR<T,3>::Dot_Product(rigid_body_list(i)->frame.t, tangent) << std::endl;
                else if (plot_velocity) std::cout << time << " " << -VECTOR<T,3>::Dot_Product(rigid_body_list(i)->velocity, tangent) << std::endl;
                else if (get_max)
                {
                    if (rigid_body_list(i)->velocity.Magnitude() > max_speed)
                    {
                        max_speed = rigid_body_list(i)->velocity.Magnitude();
                        max_speed_index = i;
                    }

                    rigid_body_list(i)->Update_Angular_Velocity();
                    if (rigid_body_list(i)->angular_velocity.Magnitude() > max_angular_speed)
                    {
                        max_angular_speed = rigid_body_list(i)->angular_velocity.Magnitude();
                        max_angular_speed_index = i;
                    }
                }
                else
                {
                    if (body_index == 0) std::cout << "body " << i << " (" << rigid_body_list(i)->name << ")" << std::endl;
                    std::cout << "position = " << rigid_body_list(i)->frame.t << std::endl;
                    if (use_frame)
                    {
                        std::cout << "position components: ";
                        print_components(rigid_body_list(i)->frame.t, tangent, normal, binormal);
                        std::cout << std::endl;
                    }
                    std::cout << "orientation = " << rigid_body_list(i)->frame.r << std::endl;
                    std::cout << "velocity = " << rigid_body_list(i)->velocity << std::endl;
                    if (use_frame)
                    {
                        std::cout << "velocity components: ";
                        print_components(rigid_body_list(i)->velocity, tangent, normal, binormal);
                        std::cout << std::endl;
                    }
                    rigid_body_list(i)->Update_Angular_Velocity();
                    std::cout << "angular_velocity = " << rigid_body_list(i)->angular_velocity << std::endl;
                    std::cout << std::endl;
                }
            }
        }

        if (get_max)
        {
            std::cout << "body " << max_speed_index << " (" << rigid_body_list(max_speed_index)->name << ") has max speed: " << max_speed << std::endl;
            std::cout << "body " << max_angular_speed_index << " (" << rigid_body_list(max_angular_speed_index)->name << ") has max angular speed: " << max_angular_speed << std::endl;
        }
    }
}

int main(int argc, char *argv[])
{
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-3d");
    parse_args.Add_Option_Argument("-2d");
    parse_args.Add_Integer_Argument("-b", body_index, "body index");
    parse_args.Add_Vector_3D_Argument("-f", VECTOR<double,3>(0,0,0), "frame rotation");
    parse_args.Add_Double_Argument("-incline", 0, "inline angle (degrees)");
    parse_args.Add_Integer_Argument("-n", num_frames, "num frames");
    parse_args.Add_Integer_Argument("-start_frame", start_frame, "start frames");
    parse_args.Add_Option_Argument("-pv", "plot velocity");
    parse_args.Add_Option_Argument("-pp", "plot position");
    parse_args.Add_Option_Argument("-pw", "plot angular_velocity");
    parse_args.Add_Option_Argument("-m", "get max");
    parse_args.Add_Double_Argument("-r", frame_rate, "frame rate");
    parse_args.Set_Extra_Arguments(-1, "<input_directory>");

    parse_args.Parse(argc, argv);
    
    body_index = parse_args.Get_Integer_Value("-b");
    get_max = parse_args.Get_Option_Value("-m");
    num_frames = parse_args.Get_Integer_Value("-n");
    start_frame = parse_args.Get_Integer_Value("-start_frame");
    plot_position = parse_args.Get_Option_Value("-pp");
    plot_velocity = parse_args.Get_Option_Value("-pv");
    plot_angular_velocity = parse_args.Get_Option_Value("-pw");
    frame_rate = parse_args.Get_Double_Value("-r");

    if(parse_args.Num_Extra_Args()) input_directory=parse_args.Extra_Arg(1);
    else input_directory=".";

    bool type_double = false;   // float by default
    if (PARSE_ARGS::Find_And_Remove("-float", argc, argv)) type_double = false;
    if (PARSE_ARGS::Find_And_Remove("-double", argc, argv)) type_double = true;

    bool do_3d = true;
    if (PARSE_ARGS::Find_And_Remove("-3d", argc, argv)) do_3d = true;
    if (PARSE_ARGS::Find_And_Remove("-2d", argc, argv)) do_3d = false;

    if(!type_double){if(do_3d)Do_It_3D<float>();else Do_It_2D<float>();}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else{if(do_3d)Do_It_3D<double>();else Do_It_2D<double>();}
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
