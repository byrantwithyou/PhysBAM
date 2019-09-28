//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Arrays/PROJECTED_ARRAY.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/RALEIGH_DAMPING_FORCE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>

#define lapack_complex_double std::complex<double>
#include <lapacke.h>

using namespace PhysBAM;

void Compute_Eigenvalues(MATRIX_MXN<std::complex<double> > M,
    ARRAY<std::complex<double> >& eig,MATRIX_MXN<std::complex<double> >& W)
{
    typedef double T;
    PHYSBAM_ASSERT(M.m==M.n);
    char jobvl = 'N';
    char jobvr = 'V';
    int n = M.m, size_work=10*n;
    ARRAY<std::complex<T> > work(size_work);
    ARRAY<T> rwork(2*n);
    int info = 0;
    M.Transpose();

    LAPACK_zgeev(&jobvl, &jobvr, &n, &M(0,0), &n, &eig(0), 0, &n, &W(0,0), &n,
        &work(0), &size_work, &rwork(0), &info);

    if(info) printf("ZGEEV FAILED: %i\n",info);
    PHYSBAM_ASSERT(info==0);
}
void Compute_Eigenvalues(MATRIX_MXN<std::complex<float> > M,
    ARRAY<std::complex<float> >& eig,MATRIX_MXN<std::complex<float> >& W)
{
    PHYSBAM_FATAL_ERROR();
}

template<class TV>
struct EXPLICIT_EXAMPLE
{
    typedef typename TV::SCALAR T;
    std::string out_dir;
    STREAM_TYPE stream_type;
    int last_frame,frame,step,resolution;
    T time,default_dt,frame_dt,dt;
    T alpha;
    T damping,poissons_ratio,stiffness;
    bool test_matrix,dump_matrix,compute_max_sound_speed,dump_sound_speed;

    EXPLICIT_EXAMPLE():stream_type((T())),last_frame(100),resolution(8),default_dt(1e-3),frame_dt(1.0/24),alpha(1),
        damping(0),poissons_ratio(0.45),stiffness(1e4),
        test_matrix(false),dump_matrix(false),compute_max_sound_speed(false),dump_sound_speed(false)
    {}
};

template<class T,int d>
T Max_Sound_Speed(const FINITE_VOLUME<VECTOR<T,d>,d>& fvm,SOLID_BODY_COLLECTION<VECTOR<T,d> >& sbc,bool dump)
{
    typedef VECTOR<T,d> TV;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_OBJECT;
    auto st=sbc.deformable_body_collection.template Find_Structure<T_OBJECT*>(0);
    DEFORMABLE_PARTICLES<TV>& particles=sbc.deformable_body_collection.particles;
    INTERPOLATED_COLOR_MAP<T> icm;
    icm.Initialize_Colors((T)1e-8,1e3,false,false,false);
    T max_c=0;
    //SYMMETRIC_MATRIX<T,TV::m> H; // H(i,k) = x_iikk
    //typename TV::SPIN B,C; // B = x_ikik; C = x_ikki; order: (1D: none; 2D: 01; 3D: 12 20 01)
    for(int t:fvm.force_elements){
        const DIAGONAL_MATRIX<T,d>& F=fvm.Fe_hat(t);
        const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dPdF=(*fvm.dPi_dFe)(t);
        T rho=fvm.use_uniform_density?fvm.density:(*fvm.density_list)(t),J=F.Determinant();
        SYMMETRIC_MATRIX<T,d> A=SYMMETRIC_MATRIX<T,d>::Conjugate(F,dPdF.H);
        T k=A.Eigenvalues().x.Max_Abs();
        for(int j=0;j<TV::SPIN::m;j++){
            int i=(j+1)%d,l=(j+2)%d;
            SYMMETRIC_MATRIX<T,2> B;
            B(0,0)=dPdF.B(j)*sqr(F(i));
            B(1,1)=dPdF.B(j)*sqr(F(l));
            B(0,1)=dPdF.C(j)*F(i)*F(l);
            k=max(k,B.Eigenvalues().x.Max_Abs());}
        T c=sqrt(k/(J*rho));
        if(dump){
            VECTOR<int,d+1>& elem=st->mesh.elements(t);
            TV X=particles.X.Subset(elem).Average();
            VECTOR<T,3> color=icm(c);
            Add_Debug_Particle(X,color);
            Debug_Particle_Set_Attribute<TV>("V",c+TV());}
        max_c=max(max_c,c);
    }
    return max_c;
}

template<class TV>
void Step(SOLID_BODY_COLLECTION<TV>& sbc,EXPLICIT_EXAMPLE<TV>& example)
{
    typedef typename TV::SCALAR T;
    sbc.Update_Position_Based_State(example.time,false,false);
    for(int i=0;example.compute_max_sound_speed && i<sbc.deformable_body_collection.deformables_forces.m;i++){
        FINITE_VOLUME<TV,TV::m>* fvm=dynamic_cast<FINITE_VOLUME<TV,TV::m>*>(sbc.deformable_body_collection.deformables_forces(i));
        if(fvm) LOG::printf("Max c: %f\n",Max_Sound_Speed(*fvm,sbc,example.dump_sound_speed));}
    DEFORMABLE_PARTICLES<TV>& particles=sbc.deformable_body_collection.particles;
    ARRAY<TV> col(particles.number);
    auto idx=[](int p,int dim){return p*TV::m+dim;};
    int base=particles.number*TV::m;
    MATRIX_MXN<T> M;

    for(int i=0;i<sbc.deformable_body_collection.deformables_forces.m;i++)
        if(LAGGED_FORCE<TV>* lf=dynamic_cast<LAGGED_FORCE<TV>*>(sbc.deformable_body_collection.deformables_forces(i)))
            lf->Lagged_Update_Position_Based_State(example.time);

    if(example.dump_matrix || example.test_matrix){
        MATRIX_MXN<T> A(particles.number*TV::m,particles.number*TV::m);
        M.Resize(particles.number*TV::m*2,particles.number*TV::m*2);
        ARRAY<TV> V(particles.number);
        for(int i=0;i<particles.number;i++){
            for(int j=0;j<TV::m;j++){
                V*=0;
                V(i)(j)=1;
                col*=0;
                sbc.deformable_body_collection.Add_Implicit_Velocity_Independent_Forces(V,col,example.time);
                for(int k=0;k<particles.number;k++){
                    for(int l=0;l<TV::m;l++){
                        T a=col(k)(l)/particles.mass(k);
                        M(idx(k,l),idx(i,j))=(idx(k,l)==idx(i,j)?1:0)+example.dt*example.dt*example.alpha*a;
                        M(base+idx(k,l),idx(i,j))=example.dt*a;
                        if(idx(k,l)==idx(i,j)){
                            M(idx(k,l),base+idx(i,j))=example.dt;
                            M(base+idx(k,l),base+idx(i,j))=1;}
                        A(idx(k,l),idx(i,j))=col(k)(l);}}}}
        OCTAVE_OUTPUT<T>(LOG::sprintf("%s/A%d.txt",example.out_dir,example.step).c_str()).Write("A",A);
        OCTAVE_OUTPUT<T>(LOG::sprintf("%s/M%d.txt",example.out_dir,example.step).c_str()).Write("M",M);
        ARRAY<std::complex<T> > eig(M.m);
        MATRIX_MXN<std::complex<T> > cM(M.m),W(M.m),P(M.m);
        for(int i=0;i<M.x.m;i++) cM.x(i)=(std::complex<T>)M.x(i);
        Compute_Eigenvalues(cM,eig,W);
        eig.Sort([](const auto& a,const auto& b){return std::abs(a)>std::abs(b);});
        P=cM-eig(0);
        ARRAY<std::complex<T> > v(M.m),v1(M.m);
        int n=0;
        for(int i=0;i<W.m;i++){
            for(int j=0;j<W.m;j++) v(j)=W(i,j);
            v1=P*v;
            T r=0;
            for(int j=0;j<v1.m;j++) r=max(r,abs(v1(j)));
            if(r<1e-12){
                //LOG::printf("limiting eigenvector: %P\n",v);
                ++n;}}
        LOG::printf("largest eigenvalue: %P magnitude: %f multiplicity: %d\n",eig(0),abs(eig(0)),n);}

    ARRAY<T> U,U1;
    if(example.test_matrix){
        U.Resize(particles.number*TV::m*2);
        U1.Resize(particles.number*TV::m*2);
        for(int i=0;i<particles.number;i++){
            for(int j=0;j<TV::m;j++){
                U(idx(i,j))=particles.X(i)(j);
                U(base+idx(i,j))=particles.V(i)(j);}}
        U1=M*U;}

    col*=0;
    sbc.Add_Velocity_Independent_Forces(col,ARRAY_VIEW<TWIST<TV> >(),example.time);
    for(int i=0;i<particles.number;i++){
        particles.X(i)+=example.dt*example.dt*example.alpha/particles.mass(i)*col(i)+example.dt*particles.V(i);
        particles.V(i)+=example.dt/particles.mass(i)*col(i);}

    if(example.test_matrix){
        T x_error=-1,v_error=-1;
        for(int i=0;i<particles.number;i++){
            for(int j=0;j<TV::m;j++){
                x_error=max(x_error,abs(U1(idx(i,j))-particles.X(i)(j)));
                v_error=max(v_error,abs(U1(base+idx(i,j))-particles.V(i)(j)));
            }
        }
        LOG::printf("L-inf error x: %P v: %P\n",x_error,v_error);}
}

template<class T,int d,class T_OBJECT>
void Add_Constitutive_Model(SOLID_BODY_COLLECTION<VECTOR<T,d> >& sbc,EXPLICIT_EXAMPLE<VECTOR<T,d> >& example,T_OBJECT& st)
{
    typedef VECTOR<T,d> TV;
    sbc.Add_Force(Create_Finite_Volume(st,new COROTATED_FIXED<T,d>(example.stiffness,example.poissons_ratio,example.damping)));
    DEFORMABLE_PARTICLES<TV>& particles=sbc.deformable_body_collection.particles;
    if(example.damping){
        DEFORMABLES_FORCES<TV>* force=Create_Finite_Volume(st,new COROTATED_FIXED<T,d>(
            example.stiffness,example.poissons_ratio,example.damping));
        force->Update_Position_Based_State(0,true,true);
        sbc.Add_Force(new RALEIGH_DAMPING_FORCE<TV>(particles,force,example.damping,1,example.dt));}
}

template<class T_OBJECT,class TV>
typename TV::SCALAR Min_Edge(const T_OBJECT& st,ARRAY_VIEW<TV> X)
{
    typedef typename TV::SCALAR T;
    SEGMENT_MESH& segment_mesh=st.mesh.Get_Segment_Mesh();
    T len=1e6;
    for(int i=0;i<segment_mesh.elements.m;i++)
        len=min(len,(X(segment_mesh.elements(i)(0))-X(segment_mesh.elements(i)(1))).Magnitude());
    return len;
}

template<class T>
void Setup(SOLID_BODY_COLLECTION<VECTOR<T,3> >& sbc,EXPLICIT_EXAMPLE<VECTOR<T,3> >& example,
    SOLIDS_STANDARD_TESTS<VECTOR<T,3> >& tests,int test_case)
{
    typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<T,3> TV;
    T density=1;
    RIGID_BODY_STATE<TV> initial_state;
    std::string data_directory="../../Public_Data";
    switch(test_case){
    case 0:{
            TETRAHEDRALIZED_VOLUME<T>& st=tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
                initial_state,true,true,density,1);
            Add_Constitutive_Model(sbc,example,st);
            sbc.deformable_body_collection.particles.X.template Project<T,&TV::y>()*=0.9;
            LOG::printf("Min edge: %f\n",Min_Edge(st,sbc.deformable_body_collection.particles.X));}
        break;
    case 1:{
            GRID<TV> mattress_grid(TV_INT()+example.resolution,RANGE<TV>(TV(),TV()+1),true);
            TETRAHEDRALIZED_VOLUME<T>& st=tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            Add_Constitutive_Model(sbc,example,st);
            sbc.deformable_body_collection.particles.X.template Project<T,&TV::y>()*=0.9;
            LOG::printf("Min edge: %f\n",Min_Edge(st,sbc.deformable_body_collection.particles.X));}
        break;
    default:PHYSBAM_FATAL_ERROR();}
}
template<class T>
void Setup(SOLID_BODY_COLLECTION<VECTOR<T,2> >& sbc,EXPLICIT_EXAMPLE<VECTOR<T,2> >& example,
    SOLIDS_STANDARD_TESTS<VECTOR<T,2> >& tests,int test_case)
{
    typedef VECTOR<int,2> TV_INT;
    typedef VECTOR<T,2> TV;
    T density=1;
    RIGID_BODY_STATE<TV> initial_state;
    std::string data_directory="../../Public_Data";
    switch(test_case){
    case 1:{
            GRID<TV> mattress_grid(TV_INT()+example.resolution,RANGE<TV>(TV(),TV()+1),true);
            TRIANGULATED_AREA<T>& st=tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            Add_Constitutive_Model(sbc,example,st);
            sbc.deformable_body_collection.particles.X.template Project<T,&TV::y>()*=0.9;
            LOG::printf("Min edge: %f\n",Min_Edge(st,sbc.deformable_body_collection.particles.X));}
        break;
    case 2:{
            GRID<TV> mattress_grid(TV_INT(8,1)*example.resolution,RANGE<TV>(TV(0,0),TV(8,1)),true);
            TRIANGULATED_AREA<T>& st=tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            Add_Constitutive_Model(sbc,example,st);
            sbc.deformable_body_collection.particles.X.template Project<T,&TV::x>()*=1.1;
            LOG::printf("Min edge: %f\n",Min_Edge(st,sbc.deformable_body_collection.particles.X));}
        break;
    default:PHYSBAM_FATAL_ERROR();}
}

template<class TV>
void Run(PARSE_ARGS& parse_args,STREAM_TYPE stream_type,const std::string& output_dir)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    EXPLICIT_EXAMPLE<TV> example;
    example.out_dir=output_dir;
    int test_number;
    bool write_substep=false;
    parse_args.Add("-test_matrix",&example.test_matrix,"Test matrix");
    parse_args.Add("-dump_matrix",&example.dump_matrix,"Dump transit matrix");
    parse_args.Add("-substep",&write_substep,"Write substep");
    parse_args.Add("-dump_c",&example.dump_sound_speed,"Dump sound speed");
    parse_args.Add("-compute_max_c",&example.compute_max_sound_speed,"Compute sound speed");
    parse_args.Add("-frame_dt",&example.frame_dt,"rate","Number of frames per second");
    parse_args.Add("-dt",&example.default_dt,"time","dt");
    parse_args.Add("-last_frame",&example.last_frame,"frame","number of frames to simulate");
    parse_args.Add("-resolution",&example.resolution,"resolution","resolution");
    parse_args.Add("-a",&example.alpha,"alpha","alpha");
    parse_args.Add("-damping",&example.damping,"damping","damping");
    parse_args.Add("-stiffness",&example.stiffness,"stiffness","stiffness");
    parse_args.Add("-poissons_ratio",&example.poissons_ratio,"poissons_ratio","poissons_ratio");
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Parse();

    std::string data_directory="../../Public_Data";
    GRID<TV> grid;
    VIEWER_OUTPUT<TV> vo(stream_type,grid,output_dir);
    SOLID_BODY_COLLECTION<TV> sbc;
    SOLIDS_STANDARD_TESTS<TV> tests(stream_type,data_directory,sbc);
    Setup(sbc,example,tests,test_number);
    Flush_Frame<TV>("init");
    sbc.Write(stream_type,output_dir,vo.frame-1,0,true,true,true,true,true);

    sbc.Update_Simulated_Particles();
    example.time=0;
    example.step=0;
    for(example.frame=0;example.frame<example.last_frame;example.frame++){
        T frame_end=(example.frame+1)*example.frame_dt;
        while(example.time<frame_end){
            example.dt=example.default_dt;
            T end_time=example.time+example.dt;
            if(end_time>=frame_end){
                example.dt=frame_end-example.time;
                end_time=frame_end;}
            Step(sbc,example);
            if(write_substep){
                Flush_Frame<TV>("Substep");
                sbc.Write(stream_type,output_dir,vo.frame-1,0,true,true,true,true,true);}
            example.time=end_time;
            example.step++;}
        if(!write_substep) Flush_Frame<TV>("Frame");
        sbc.Write(stream_type,output_dir,vo.frame-1,0,true,true,true,true,true);}
}

int main(int argc, char* argv[])
{
    bool type_double=true,use_3d=false;
    std::string output_dir="output";
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-o",&output_dir,"dir","Output directory");
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-3d",&use_3d,"Run 3D examples");
    parse_args.Parse(true);

    Create_Directory(output_dir+"/common");
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::printf("%s\n",parse_args.Print_Arguments());
    if(type_double){
        STREAM_TYPE stream_type((double()));
        if(use_3d) Run<VECTOR<double,3> >(parse_args,stream_type,output_dir);
        else Run<VECTOR<double,2> >(parse_args,stream_type,output_dir);}
    else{
        STREAM_TYPE stream_type((float()));
        if(use_3d) Run<VECTOR<float,3> >(parse_args,stream_type,output_dir);
        else Run<VECTOR<float,2> >(parse_args,stream_type,output_dir);}
    LOG::Finish_Logging();
    return 0;
}

