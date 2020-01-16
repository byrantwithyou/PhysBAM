//#####################################################################
// Copyright 2015, Andre Pradhana and Gergely Klar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Partio.h>
#include <regex>
#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace PhysBAM;

HASHTABLE<std::string,std::string> name_map;

std::string get_name(const std::string& name)
{
    std::string s;
    if(name_map.Get(name,s)) return s;
    return name;
}

template<class U,int d,class TV,class IT>
bool Try_Write_Vector(Partio::ParticlesDataMutable* parts,const MPM_PARTICLES<TV>& particles,ATTRIBUTE_INDEX id,const IT& it)
{
    ARRAY_COLLECTION_ELEMENT<VECTOR<U,d> >* c=dynamic_cast<ARRAY_COLLECTION_ELEMENT<VECTOR<U,d> >*>(particles.arrays(id));
    if(!c) return false;
    ARRAY_VIEW<VECTOR<U,d> > a(*c->array);
    Partio::ParticleAttribute pa=parts->addAttribute(get_name(c->name).c_str(),Partio::VECTOR,d);
    IT k(it);
    for(int i=0;i<particles.number;i++)
    {
        VECTOR<U,d> w=a(i);
        float *v=parts->dataWrite<float>(pa,k.index);
        for(int j=0;j<d;j++) v[j]=w(j);
        ++k;
    }
    LOG::printf("Wrote attribute %P as vector of size %P named %P.\n",c->name,d,get_name(c->name));
    return true;
}

template<class U,class TV,class IT>
bool Try_Write_Float(Partio::ParticlesDataMutable* parts,const MPM_PARTICLES<TV>& particles,ATTRIBUTE_INDEX id,const IT& it)
{
    ARRAY_COLLECTION_ELEMENT<U>* c=dynamic_cast<ARRAY_COLLECTION_ELEMENT<U>*>(particles.arrays(id));
    if(!c) return false;
    ARRAY_VIEW<U> a(*c->array);
    Partio::ParticleAttribute pa=parts->addAttribute(get_name(c->name).c_str(),Partio::FLOAT,1);
    IT k(it);
    for(int i=0;i<particles.number;i++)
    {
        *parts->dataWrite<float>(pa,k.index)=a(i);
        ++k;
    }
    LOG::printf("Wrote attribute %P as float named %P.\n",c->name,get_name(c->name));
    return true;
}

template<class U,class TV,class IT>
bool Try_Write_Int(Partio::ParticlesDataMutable* parts,const MPM_PARTICLES<TV>& particles,ATTRIBUTE_INDEX id,const IT& it)
{
    ARRAY_COLLECTION_ELEMENT<U>* c=dynamic_cast<ARRAY_COLLECTION_ELEMENT<U>*>(particles.arrays(id));
    if(!c) return false;
    ARRAY_VIEW<U> a(*c->array);
    Partio::ParticleAttribute pa=parts->addAttribute(get_name(c->name).c_str(),Partio::INT,1);
    IT k(it);
    for(int i=0;i<particles.number;i++)
    {
        *parts->dataWrite<int>(pa,k.index)=a(i);
        ++k;
    }
    LOG::printf("Wrote attribute %P as int named %P.\n",c->name,get_name(c->name));
    return true;
}

template<class U,int d,class TV,class IT>
bool Try_Write_Matrix(Partio::ParticlesDataMutable* parts,const MPM_PARTICLES<TV>& particles,ATTRIBUTE_INDEX id,const IT& it)
{
    ARRAY_COLLECTION_ELEMENT<MATRIX<U,d> >* c=dynamic_cast<ARRAY_COLLECTION_ELEMENT<MATRIX<U,d> >*>(particles.arrays(id));
    if(!c) return false;
    ARRAY_VIEW<MATRIX<U,d> > a(*c->array);
    Partio::ParticleAttribute pa[d];
    VECTOR<std::string,d> names;
    for(int i=0;i<d;i++)
    {
        names(i)=get_name(c->name+"_"+"xyzabcd"[i]);
        pa[i]=parts->addAttribute(names(i).c_str(),Partio::VECTOR,d);
    }
    IT k(it);
    for(int i=0;i<particles.number;i++)
    {
        MATRIX<U,d> w=a(i);
        for(int l=0;l<d;l++)
        {
            float *v=parts->dataWrite<float>(pa[l],k.index);
            for(int j=0;j<d;j++) v[j]=w(j,l);
        }
        ++k;
    }
    LOG::printf("Wrote attribute %P as column vectors of size %P named %P.\n",c->name,d,names);
    return true;
}

template<class T,int N>
void writePartio(const std::string& output_filename,const MPM_PARTICLES<VECTOR<T,N> >& particles,bool dump_valid)
{
    Partio::ParticlesDataMutable* parts=Partio::create();
    auto it=parts->addParticles(particles.number);
    for(ATTRIBUTE_INDEX id(0);id<particles.arrays.m;id++)
    {
        if(Try_Write_Float<T>(parts,particles,id,it)) continue;
        if(Try_Write_Int<int>(parts,particles,id,it)) continue;
        if(Try_Write_Int<bool>(parts,particles,id,it)) continue;
        if(Try_Write_Vector<T,1>(parts,particles,id,it)) continue;
        if(Try_Write_Vector<T,2>(parts,particles,id,it)) continue;
        if(Try_Write_Vector<T,3>(parts,particles,id,it)) continue;
        if(Try_Write_Matrix<T,1>(parts,particles,id,it)) continue;
        if(Try_Write_Matrix<T,2>(parts,particles,id,it)) continue;
        if(Try_Write_Matrix<T,3>(parts,particles,id,it)) continue;
        LOG::printf("Attribute %P has unsupported type.\n",particles.arrays(id)->name);
    }
    Partio::write(output_filename.c_str(),*parts);
    parts->release();

    LOG::printf("Exported %d particles to %s\n",particles.number,output_filename);
}

template<class T,int N>
void Convert(const std::string& input,const std::string& output_filename_pattern,int start_at,bool dump_valid,bool attentive)
{
    bool has_format_str=std::regex_match(output_filename_pattern,std::regex(".*%[0-9]*d.*"));
    if(Directory_Exists(input)){
        if(!has_format_str){
            LOG::printf("Format string missing from output file name! E.g. %%04d\n");
            exit(-1);}

        int last_frame;
        Read_From_Text_File(input+"/common/last_frame",last_frame);
        PHYSBAM_ASSERT(start_at>=0 && start_at<=last_frame);
#pragma omp parallel for
        for(int i=start_at;i<=last_frame;++i){
            MPM_PARTICLES<VECTOR<T,N> > particles;
            Read_From_File(LOG::sprintf("%s/%d/mpm_particles.gz",input,i),particles);
            writePartio<T,N>(LOG::sprintf(output_filename_pattern.c_str(),i),particles,dump_valid);
            if(attentive)
                Read_From_Text_File(input+"/common/last_frame",last_frame);}}
    else{
        if(has_format_str){
            LOG::printf("Format string found in output file name! Did you want to convert a whole sim?\n");
            exit(-1);}
        MPM_PARTICLES<VECTOR<T,N> > particles;
        Read_From_File(input,particles);
        writePartio<T,N>(output_filename_pattern,particles,dump_valid);}
}

int main(int argc,char *argv[])
{
    bool type_double=true;
    bool dump_valid=false;
    bool use_3d=false;
    bool attentive=true;
    int start_at=0,threads=1;
    std::string input,output_filename_pattern;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-dump_valid_flag",&dump_valid,"Dump all particles along with a valid flag; don't prune invalid particles");
    parse_args.Add("-3d",&use_3d,"Convert 3D examples");
    parse_args.Add("-start",&start_at,"start_at","Start conversion from this frame");
    parse_args.Add_Not("-no_attentive",&attentive,"Don't look for new frames during conversion");
    parse_args.Add("-threads",&threads,"threads","Number of threads");
    parse_args.Extra(&input,"input","Simulation output directory or single file");
    parse_args.Extra(&output_filename_pattern,"output","Output partio file name. Use format string (e.g. %04d), if converting whole directory");
    parse_args.Parse();

    name_map.Set("X","position");
    
#ifdef USE_OPENMP
    omp_set_num_threads(threads);
#pragma omp parallel
#pragma omp single
    {
        if(omp_get_num_threads()!=threads) PHYSBAM_FATAL_ERROR();
        LOG::cout<<"Running on "<<threads<<" threads"<<std::endl;
    }
#else
    PHYSBAM_ASSERT(threads==1);
#endif

    if(output_filename_pattern.length()<5||output_filename_pattern.substr(output_filename_pattern.length()-5,5)!=".bgeo"){
        LOG::printf("Missing \".bgeo\" from output file name. Appending.\n");
        output_filename_pattern+=".bgeo";}

    if(type_double)
        if(use_3d) Convert<double,3>(input,output_filename_pattern,start_at,dump_valid,attentive);
        else Convert<double,2>(input,output_filename_pattern,start_at,dump_valid,attentive);
    else
        if(use_3d) Convert<float,3>(input,output_filename_pattern,start_at,dump_valid,attentive);
        else Convert<float,2>(input,output_filename_pattern,start_at,dump_valid,attentive);

    return 0;
}
