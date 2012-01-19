#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_ROTATION.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void
Get_Mass_Properties(const std::string& filename,PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,3> TV;

    LOG::SCOPE scope("mass properties","mass properties");
    LOG::cout<<"filename = "<<filename<<std::endl;
    LOG::cout<<"using "<<(IS_SAME<T,float>::value?"floats":"doubles")<<std::endl;

    TRIANGLE_MESH mesh;
    PARTICLES<TV> particles;
    TRIANGULATED_SURFACE<T> surface(mesh,particles);
    FILE_UTILITIES::Read_From_File<RW>(filename,surface);
    surface.Update_Triangle_List();

    RIGID_BODY_COLLECTION<TV> rigid_body_collection(0,0);
    RIGID_BODY<TV> rigid_body(rigid_body_collection);
    rigid_body.thin_shell=parse_args.Get_Option_Value("-thin_shell");
    rigid_body.surface_roughness=(T)1e-8;
    SYMMETRIC_MATRIX<T,3> world_space_inertia_tensor;

    MASS_PROPERTIES<TV> mass_properties(surface,parse_args.Get_Option_Value("-thin_shell"));
    if(parse_args.Is_Value_Set("-density")) mass_properties.Set_Density((T)parse_args.Get_Double_Value("-density"));
    else mass_properties.Set_Mass((T)parse_args.Get_Double_Value("-mass"));
    rigid_body.Mass()=mass_properties.Mass();
    FRAME<TV> frame(rigid_body.Frame());
    mass_properties.Transform_To_Object_Frame(frame,rigid_body.Inertia_Tensor());

    LOG::cout<<"mass = "<<rigid_body.Mass()<<std::endl;
    LOG::cout<<(rigid_body.thin_shell?"surface area = ":"volume = ")<<mass_properties.Volume()<<std::endl;
    LOG::cout<<"center of mass = "<<rigid_body.Frame().t<<std::endl;
    LOG::cout<<"orientation = "<<rigid_body.Frame().r<<std::endl;
    //LOG::cout<<"world space inertia_tensor = "<<mass_properties.Inertia_Tensor()<<std::endl;
    LOG::cout<<"object space inertia_tensor = "<<rigid_body.Inertia_Tensor()<<std::endl;

    if(parse_args.Get_Option_Value("-generate")){
        {LOG::SCOPE scope("transforming surface","transforming surface");
        for(int p=0;p<surface.particles.array_collection->Size();p++) surface.particles.X(p)=rigid_body.Frame().Inverse_Times(surface.particles.X(p));
        FILE_UTILITIES::Write_To_File<RW>(filename,surface);}

        if(parse_args.Is_Value_Set("-secondary_surface")){
            LOG::SCOPE scope("transforming secondary surface","transforming secondary surface");
            std::string secondary_filename=parse_args.Get_String_Value("-secondary_surface");
            LOG::cout<<"filename = "<<secondary_filename<<std::endl;
            FILE_UTILITIES::Read_From_File<RW>(secondary_filename,surface);
            for(int p=0;p<surface.particles.array_collection->Size();p++) surface.particles.X(p)=rigid_body.Frame().Inverse_Times(surface.particles.X(p));
            FILE_UTILITIES::Write_To_File<RW>(secondary_filename,surface);}

        LOG::SCOPE scope("generating rgd file","generating rgd file");
        FILE_UTILITIES::Write_To_File<RW>(FILE_UTILITIES::Get_Basename(filename)+".rgd",rigid_body.Mass(),rigid_body.Frame());}
}

int main(int argc,char *argv[])
{
    bool type_double=false,write_double=false;

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-write_float");
    parse_args.Add_Option_Argument("-write_double");
    parse_args.Add_Option_Argument("-thin_shell");
    parse_args.Add_Option_Argument("-generate");
    parse_args.Add_Double_Argument("-mass",1,"mass");
    parse_args.Add_Double_Argument("-density",1,"density");
    parse_args.Add_String_Argument("-secondary_surface","");
    parse_args.Set_Extra_Arguments(1,"<filename>");
    parse_args.Parse(argc,argv);

    std::string filename=parse_args.Extra_Arg(1);

    if(parse_args.Get_Option_Value("-float")) type_double=false;
    else if(parse_args.Get_Option_Value("-double")) type_double=true;

    if(parse_args.Get_Option_Value("-write_float")) write_double=false;
    else if(parse_args.Get_Option_Value("-write_double")) write_double=true;

    if(!type_double && !write_double) Get_Mass_Properties<float,float>(filename,parse_args);
#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else PHYSBAM_FATAL_ERROR("double support not enabled.");
#else
    else if(!type_double && write_double) Get_Mass_Properties<float,double>(filename,parse_args);
    else if(type_double && !write_double) Get_Mass_Properties<double,float>(filename,parse_args);
    else Get_Mass_Properties<double,double>(filename,parse_args);
#endif
    LOG::cout<<std::flush;
    return 0;
}
