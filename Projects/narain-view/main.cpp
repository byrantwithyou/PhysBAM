//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_DILATE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <fstream>
#include <btBulletCollisionCommon.h>
#include <btBulletDynamicsCommon.h>
#include <btBulletWorldImporter.h>

using namespace PhysBAM;

typedef float T;

template<int d>
void Bullet_To_PhysBAM(const btVector3& v,VECTOR<float,d>& o)
{
    for(int i=0;i<d;i++) o(i)=v.m_floats[i];
}

void Bullet_To_PhysBAM(const btQuaternion& q,ROTATION<VECTOR<float,3> >& r)
{
    r.q=QUATERNION<T>(q[0],q[1],q[2],q[3]);
}

void Bullet_To_PhysBAM(const btQuaternion& q,ROTATION<VECTOR<float,2> >& r)
{
    r=ROTATION<VECTOR<float,2> >::From_Angle(q.getAngle()*q.getAxis()[2]);
}

enum modType {SPLIT, MERGE, MOVE};
struct SplitMerge
{
    modType m;
    int pID[2];
};

template<class TV>
LEVELSET_IMPLICIT_OBJECT<TV>* Read_Implicit_Object(const std::string& filename)
{
    static HASHTABLE<std::string,LEVELSET_IMPLICIT_OBJECT<TV>*> implicit_object_hash;
    if(LEVELSET_IMPLICIT_OBJECT<TV>** io=implicit_object_hash.Get_Pointer(filename))
        return *io;
    typedef VECTOR<int,TV::m> TV_INT;
    TV_INT size;
    T dx;
    std::fstream file(filename.c_str(), std::ios::in|std::ios::binary);
    file.read((char*)&size, sizeof(size));
    file.read((char*)&dx, sizeof(dx));
    ARRAY<T> data(size.Size());
    Read_Binary_Array<T>(file,data.Get_Array_Pointer(),data.m);
    LEVELSET_IMPLICIT_OBJECT<TV> *lio=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    lio->levelset.grid.Initialize(size,RANGE<TV>(TV(),TV(size)*dx),false);
    lio->levelset.phi.Resize(RANGE<TV_INT>(TV_INT(),size));
    TV_INT stride(TV_INT()+1);
    for(int i=1;i<TV_INT::m;i++)
        stride(i)=stride(i-1)*size(i);
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),size));it.Valid();it.Next())
        lio->levelset.phi(it.index)=data(it.index.Dot(stride));
    LOG::printf("DUMP:%P %P %P\n%P\n",size,dx,filename,lio->levelset.phi);
    implicit_object_hash.Set(filename,lio);
    return lio;
}

template<class T>
void Read_Rotation(ROTATION<VECTOR<T,2> >& rot,std::fstream& in)
{
    T q=0;
    in.read((char*)&q, 4);
    rot=ROTATION<VECTOR<T,2> >::From_Angle(q);
}

template<class T>
void Read_Rotation(ROTATION<VECTOR<T,3> >& rot,std::fstream& in)
{
    in.read((char*)&rot, sizeof(rot));
}

template<class TV>
void Narain_To_PhysBAM(PARSE_ARGS& parse_args)
{
    typedef VECTOR<int,TV::m> TV_INT;
    static const int d=TV::m;
    std::string input_directory;
    std::string input_name;
    std::string models_directory;
    std::string output_directory;
    parse_args.Add("-i",&input_directory,"dir","input simulation directory");
    parse_args.Add("-n",&input_name,"name","input simulation name");
    parse_args.Add("-m",&models_directory,"dir","models directory");
    parse_args.Add("-o",&output_directory,"dir","output simulation directory");
    parse_args.Parse();

    GRID<TV> grid;
    SOLID_BODY_COLLECTION<TV> solid_body_collection;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particles;
    DEBUG_PARTICLES<TV> debug_particles;
    ARRAY_VIEW<MATRIX<T,3> > particle_A;
    deformable_body_collection.particles.Add_Array(ATTRIBUTE_ID(100),&particle_A);
    ARRAY<T,FACE_INDEX<d> > face_velocities;

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/first_frame",0);

    for(int frame=0;/*TODO*/;frame++)
    {
        std::string simfile = LOG::sprintf("%s/%s-%04i.sand",input_directory,input_name,frame);
        std::fstream file(simfile.c_str(), std::ios::in|std::ios::binary);
        if(!file)
        {
            break;
        }
        int32_t frame_=0,zb=0,np=0,nb=0,nsm=0;
        file.read((char*)&frame_, 4);
        PHYSBAM_ASSERT(frame==frame_);
        T dt=0,simtime=0,dx=0,rho_max=0,r=0,mu=0;
        TV gravity;
        TV_INT size;
        char buff[100];
        file.read((char*)&dt, 4);
        file.read((char*)&simtime, 4);
        file.read((char*)&gravity, d*4);
        file.read((char*)&size, d*4);
        file.read((char*)&dx, 4);
        file.read((char*)&rho_max, 4);
        file.read((char*)&r, 4);
        file.read((char*)&mu, 4);
        file.read((char*)&zb, 4);
        grid.Initialize(size,RANGE<TV>(TV(),dx*TV(size)),true);
        face_velocities.Resize(grid);
        file.read((char*)&np, 4);
        deformable_body_collection.particles.Resize(np);
        for(int i = 0; i < np; i++)
        {
            file.read((char*)&deformable_body_collection.particles.X(i), d*4);
            file.read((char*)&deformable_body_collection.particles.V(i), d*4);
            file.read((char*)&deformable_body_collection.particles.mass(i), 4);
            file.read((char*)&particle_A(i), d*d*4);
        }

        ARRAY<FRAME<TV> > bullet_frames;
        {
            // load Bullet state
            btDefaultCollisionConfiguration cc;
            btCollisionDispatcher dispatcher(&cc);
            btDbvtBroadphase broadphase;
            btSequentialImpulseConstraintSolver solver;
            btDiscreteDynamicsWorld bt(&dispatcher, &broadphase, &solver, &cc);
            btBulletWorldImporter importer(&bt);
            std::string btfile = LOG::sprintf("%s/%s-%04i.bullet",input_directory,input_name,frame);
            PHYSBAM_ASSERT(importer.loadFile(btfile.c_str()));
            int nc = bt.getCollisionObjectArray().size();
            for(int c = 0; c < nc; c++)
            {
                FRAME<TV> frame;
                btCollisionObject* brb=bt.getCollisionObjectArray()[c];
                Bullet_To_PhysBAM(brb->getWorldTransform().getOrigin(),frame.t);
                Bullet_To_PhysBAM(brb->getWorldTransform().getRotation(),frame.r);
                bullet_frames.Append(frame);
            }
        }

        file.read((char*)&nb, 4);
        for(int i = 0, bf = 0; i < nb; i++)
        {
            file.read(buff, 8);
            LOG::printf("dynamic %i\n",buff[0]);
            FRAME<TV> frame;
            if(buff[0]==0)
            {
                file.read((char*)&frame.t, d*4);
                Read_Rotation(frame.r,file);
            }
            else
            {
                frame=bullet_frames(bf++);
            }
            int32_t len;
            char src[256]="";
            T scale,offset;
            file.read((char*)&len, 4);
            file.read(src, len);

            file.read((char*)&scale, 4);
            file.read((char*)&offset, 4);
            // TODO:  scale or offset first?

            if(!rigid_body_collection.Exists(i))
            {
                RIGID_BODY<TV>* rb=new RIGID_BODY<TV>(rigid_body_collection,false,i);
                rb->thin_shell=false;
            
                rb->Update_Angular_Momentum();
                LEVELSET_IMPLICIT_OBJECT<TV>* lio=Read_Implicit_Object<TV>(models_directory+"/"+src+".dist");
                typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
                T_SURFACE* surface=T_SURFACE::Create();
                MARCHING_CUBES<TV>::Create_Surface(*surface,lio->levelset.grid,lio->levelset.phi,-offset*scale);
                surface->particles.X*=scale;
                LOG::printf("%P %P\n",surface->particles.X,surface->mesh.elements);
                IMPLICIT_OBJECT<TV>* io=lio;
                if(scale!=1) io=new IMPLICIT_OBJECT_TRANSFORMED<TV,T>(io,false,TV(),scale);
                if(offset) io=new IMPLICIT_OBJECT_DILATE<TV>(io,-offset);
                rb->Add_Structure(*io);

                rb->Add_Structure(*surface);

                rigid_body_collection.Add_Rigid_Body_And_Geometry(rb);
            }
            rigid_body_particles.mass(i)=1;
            rigid_body_particles.inertia_tensor(i)=DIAGONAL_MATRIX<T,TV::SPIN::m>()+1;
            rigid_body_particles.frame(i)=frame;
        }
        file.read((char*)&nsm, 4);
        ARRAY<SplitMerge> smmap(nsm);
        for(int i = 0; i < nsm; i++)
            file.read((char*)&smmap(i),sizeof(SplitMerge));

        if(d==2) file.read(buff, 4*7);

//        std::string f=LOG::sprintf("%d",frame);
        FILE_UTILITIES::Create_Directory(output_directory+LOG::sprintf("/%d",frame));
        FILE_UTILITIES::Write_To_File(STREAM_TYPE((T)0),output_directory+"/common/grid",grid);
        if(!system(LOG::sprintf("rm -f %s/%d/mpm_particles.gz ;  ln -s ./deformable_object_particles.gz %s/%d/mpm_particles.gz",output_directory.c_str(),frame,output_directory.c_str(),frame).c_str())){}
        solid_body_collection.Write(STREAM_TYPE((T)0),output_directory,frame,0,frame==0,true,true,false,false);
        FILE_UTILITIES::Write_To_File(STREAM_TYPE((T)0),LOG::sprintf("%s/%d/mac_velocities",output_directory.c_str(),frame),face_velocities);
        for(int i=0;i<deformable_body_collection.particles.X.m;i++){
            Add_Debug_Particle(deformable_body_collection.particles.X(i),VECTOR<T,3>(1,0,1));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,deformable_body_collection.particles.V(i));}
        debug_particles.Write_Debug_Particles(STREAM_TYPE((T)0),output_directory,frame);
        FILE_UTILITIES::Write_To_File(STREAM_TYPE((T)0),LOG::sprintf("%s/%d/frame_title",output_directory.c_str(),frame),"");
        FILE_UTILITIES::Write_To_File(STREAM_TYPE((T)0),output_directory+"/common/grid",grid);
        FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame-1);
    }
}

int main(int argc, char* argv[])
{
    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"Use 3D");
    parse_args.Parse(true);
    
    if(use_3d)
        Narain_To_PhysBAM<VECTOR<T,3> >(parse_args);
    else
        Narain_To_PhysBAM<VECTOR<T,2> >(parse_args);

    return 0;
}



/*

  common/last_frame
  common/grid.gz
  common/first_frame
  common/deformable_object_structures.gz (??)

  ###/debug_particles.gz
  ###/frame_title
  ###/face_velocities.gz
  ###/deformable_object_particles.gz
  ###/mpm_particles.gz -> deformable_object_particles.gz


*/



