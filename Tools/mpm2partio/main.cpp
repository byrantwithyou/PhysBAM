#include <Partio.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;
using namespace std;

template <typename T>
void writePartio(const std::string& output_filename,const MPM_PARTICLES<VECTOR<T,3> >& particles,bool dump_valid,bool dump_extra)
{
    Partio::ParticlesDataMutable* parts=Partio::create();

    Partio::ParticleAttribute posH,vH,FxH,FyH,FzH,mH,volH,validH;
    posH=parts->addAttribute("position",Partio::VECTOR,3);
    vH=parts->addAttribute("v",Partio::VECTOR,3);
    if(dump_extra){
        FxH=parts->addAttribute("Fx",Partio::VECTOR,3);
        FyH=parts->addAttribute("Fy",Partio::VECTOR,3);
        FzH=parts->addAttribute("Fz",Partio::VECTOR,3);
        volH=parts->addAttribute("vol",Partio::VECTOR,1);}
    mH=parts->addAttribute("m",Partio::VECTOR,1);
    if(dump_valid) validH=parts->addAttribute("valid",Partio::INT,1);

    for(int i=0;i<particles.number;++i){
        if(!dump_valid && !particles.valid(i)) continue;
        int idx=parts->addParticle();
        float *p=parts->dataWrite<float>(posH,idx);
        float *v=parts->dataWrite<float>(vH,idx);
        if(dump_extra){
            float *Fx=parts->dataWrite<float>(FxH,idx);
            float *Fy=parts->dataWrite<float>(FyH,idx);
            float *Fz=parts->dataWrite<float>(FzH,idx);
            float *vol=parts->dataWrite<float>(volH,idx);

            MATRIX<T,3> F=particles.F(i);
            for(int k=0;k<3;++k) Fx[k]=F(k,0);
            for(int k=0;k<3;++k) Fy[k]=F(k,1);
            for(int k=0;k<3;++k) Fz[k]=F(k,2);
            vol[0]=particles.volume(i);
        } 

        float *m=parts->dataWrite<float>(mH,idx);

        for(int k=0;k<3;++k) p[k]=particles.X(i)[k];
        for(int k=0;k<3;++k) v[k]=particles.V(i)[k];

        m[0]=particles.mass(i);

        if(dump_valid) *parts->dataWrite<int>(validH,idx)=particles.valid(i);}

    Partio::write(output_filename.c_str(),*parts);
    parts->release();

    cout << STRING_UTILITIES::string_sprintf("Exported %d particles to %s\n",particles.number,output_filename.c_str());
}

template <typename T>
void Convert(const std::string& input_filename,const std::string& output_filename,bool dump_valid,bool dump_extra)
{
    MPM_PARTICLES<VECTOR<T,3> > particles;
    FILE_UTILITIES::Read_From_File<T>(input_filename, particles);

    writePartio<T>(output_filename,particles,dump_valid,dump_extra);
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=true;
    bool dump_valid=false,dump_extra=false;
    std::string input_filename,output_filename;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-dump_valid_flag",&dump_valid,"Dump all particles along with a valid flag");
    parse_args.Add("-dump_extra",&dump_extra,"Dump extra particles parameters: F and volume");
    parse_args.Extra(&input_filename,"debug particles file","debug particles to convert");
    parse_args.Extra(&output_filename,"partio file","output partio file name");
    parse_args.Parse();

    if(dump_valid) cout<<"Dumping valid flags along particles\n";
    if(dump_extra) cout<<"Dumping F and volume\n";

    if(type_double){
        cout << "Assuming doubles in input files\n";
        Convert<double>(input_filename,output_filename,dump_valid,dump_extra);}
    else{
        cout << "Assuming floats in input files\n";
        Convert<float>(input_filename,output_filename,dump_valid,dump_extra);}
}
