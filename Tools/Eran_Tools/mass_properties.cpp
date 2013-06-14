#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void
Get_Mass_Properties(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,3> TV;

    bool opt_generate=false,opt_density=false,thin_shell=false;
    T mass=1,density=1;
    std::string secondary_filename,filename;
    parse_args.Add("-thin_shell",&thin_shell,"Thin shell");
    parse_args.Add("-generate",&opt_generate,"Generate");
    parse_args.Add("-mass",&mass,"mass","mass");
    parse_args.Add("-density",&density,&opt_density,"density","density");
    parse_args.Add("-secondary_surface",&secondary_filename,"file","secondary filename");
    parse_args.Extra(&filename,"filename","filename");
    parse_args.Parse();

    LOG::SCOPE scope("mass properties","mass properties");
    LOG::cout<<"filename = "<<filename<<std::endl;
    LOG::cout<<"using "<<(IS_SAME<T,float>::value?"floats":"doubles")<<std::endl;

    TRIANGLE_MESH mesh;
    DEFORMABLE_PARTICLES<TV> particles;
    TRIANGULATED_SURFACE<T> surface(mesh,particles);
    FILE_UTILITIES::Read_From_File<RW>(filename,surface);
    surface.Update_Triangle_List();

    RIGID_BODY_COLLECTION<TV> rigid_body_collection(0);
    RIGID_BODY<TV> rigid_body(rigid_body_collection);
    rigid_body.surface_roughness=(T)1e-8;
    SYMMETRIC_MATRIX<T,3> world_space_inertia_tensor;

    rigid_body.thin_shell=thin_shell;
    MASS_PROPERTIES<TV> mass_properties(surface,thin_shell);
    if(opt_density) mass_properties.Set_Density(density);
    else mass_properties.Set_Mass(mass);
    rigid_body.Mass()=mass_properties.Mass();
    FRAME<TV> frame(rigid_body.Frame());
    mass_properties.Transform_To_Object_Frame(frame,rigid_body.Inertia_Tensor());

    LOG::cout<<"mass = "<<rigid_body.Mass()<<std::endl;
    LOG::cout<<(rigid_body.thin_shell?"surface area = ":"volume = ")<<mass_properties.Volume()<<std::endl;
    LOG::cout<<"center of mass = "<<rigid_body.Frame().t<<std::endl;
    LOG::cout<<"orientation = "<<rigid_body.Frame().r<<std::endl;
    //LOG::cout<<"world space inertia_tensor = "<<mass_properties.Inertia_Tensor()<<std::endl;
    LOG::cout<<"object space inertia_tensor = "<<rigid_body.Inertia_Tensor()<<std::endl;

    if(opt_generate){
        {LOG::SCOPE scope("transforming surface","transforming surface");
        for(int p=0;p<surface.particles.Size();p++) surface.particles.X(p)=rigid_body.Frame().Inverse_Times(surface.particles.X(p));
        FILE_UTILITIES::Write_To_File<RW>(filename,surface);}

        if(secondary_filename.size()){
            LOG::SCOPE scope("transforming secondary surface","transforming secondary surface");
            LOG::cout<<"filename = "<<secondary_filename<<std::endl;
            FILE_UTILITIES::Read_From_File<RW>(secondary_filename,surface);
            for(int p=0;p<surface.particles.Size();p++) surface.particles.X(p)=rigid_body.Frame().Inverse_Times(surface.particles.X(p));
            FILE_UTILITIES::Write_To_File<RW>(secondary_filename,surface);}

        LOG::SCOPE scope("generating rgd file","generating rgd file");
        FILE_UTILITIES::Write_To_File<RW>(FILE_UTILITIES::Get_Basename(filename)+".rgd",rigid_body.Mass(),rigid_body.Frame());}
}

int main(int argc,char *argv[])
{
    bool type_double=false,write_double=false;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add_Not("-write_float",&type_double,"Write floats");
    parse_args.Add("-write_double",&type_double,"Write doubles");
    parse_args.Parse(true);

    if(!type_double && !write_double) Get_Mass_Properties<float,float>(parse_args);
#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else PHYSBAM_FATAL_ERROR("double support not enabled.");
#else
    else if(!type_double && write_double) Get_Mass_Properties<float,double>(parse_args);
    else if(type_double && !write_double) Get_Mass_Properties<double,float>(parse_args);
    else Get_Mass_Properties<double,double>(parse_args);
#endif
    LOG::cout<<std::flush;
    return 0;
}
