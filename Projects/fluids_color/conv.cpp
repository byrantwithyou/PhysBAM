#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Level_Sets/LEVELSET.h>

using namespace PhysBAM;

template<class T>
T Approx_Exact(const ARRAY<int>& n,const ARRAY<T>& x,int order)
{
    ARRAY<int> w(x.m);
    MATRIX_MXN<T> A(x.m,order+1);
    for(int i=0;i<x.m;i++){
        T w=(T)1/n(i), p=1;
        for(int j=0;j<order+1;j++){
            A(i,j)=p;
            p*=w;}}

    ARRAY<T> coeffs=(A.Transpose_Times(A)).Inverse_Times(A.Transpose_Times(x));
    return coeffs(0);
}

template<class TV>
struct SIM_DATA
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    T time;
    int res;
    GRID<TV> grid;
    VECTOR<GRID<TV>,TV::m> face_grids;
    ARRAY<int,FACE_INDEX<TV::dimension> > face_color;
    ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> > > face_velocities;
    ARRAY<T,TV_INT> pressure,phi;
    LEVELSET<TV> levelset;

    SIM_DATA(STREAM_TYPE& stream_type,const std::string& dir,int frame)
        :levelset(grid,phi)
    {
        ARRAY<int,FACE_INDEX<TV::m> > prev_face_color;
        ARRAY<ARRAY<T,FACE_INDEX<TV::m> > > prev_face_velocities;
        FILE_UTILITIES::Read_From_File(stream_type,dir+"/common/grid",grid);
        FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/levelset_%d",dir.c_str(),frame,0),levelset);
        FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/pressure",dir.c_str(),frame),pressure);
        FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/restart_data",dir.c_str(),frame),
            time,face_color,prev_face_color,face_velocities,prev_face_velocities);
        res=grid.counts.x;
        for(int d=0;d<TV::m;d++) face_grids(d)=grid.Get_Face_Grid(d);
    }
};

template<class TV>
bool compare_data(const SIM_DATA<TV>* a,const SIM_DATA<TV>* b) {return a->res<b->res;}

template<class TV>
void Analyze(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    STREAM_TYPE stream_type((T()));
    TV X0(TV()+FLT_MAX);
    int frame=1;
    ARRAY<std::string> sim_dirs;
    parse_args.Add("-frame",&frame,"frame","Frame to test");
    parse_args.Add("-x",&X0.x,"value","point to test");
    parse_args.Add("-y",&X0.y,"value","point to test");
    if(TV::m==3) parse_args.Add("-z",&X0(TV::m-1),"value","point to test");
    parse_args.Extra_Optional(&sim_dirs,"sim dirs","simulation directories");
    parse_args.Parse();

    ARRAY<SIM_DATA<TV>*> sim_data;

    for(int i=0;i<sim_dirs.m;i++)
    {
        try
        {
            sim_data.Append(new SIM_DATA<TV>(stream_type,sim_dirs(i),frame));
        }
        catch(...)
        {
            LOG::cout<<"Failed to load: "<<sim_dirs(i)<<std::endl;
        }
    }

    sim_data.Sort(compare_data<TV>);
    ARRAY<int> res;
    for(int i=0;i<sim_data.m;i++)
        res.Append(sim_data(i)->res);

    for(int i=0;i<sim_data.m;i++){
        T p_ave=0;
        int cnt=0;
        for(RANGE_ITERATOR<TV::m> it(sim_data(i)->grid.Domain_Indices());it.Valid();it.Next()){
            p_ave+=sim_data(i)->pressure(it.index);
            cnt++;}
        p_ave/=cnt;
        sim_data(i)->pressure.array-=p_ave;}

    ARRAY<int> samples(sim_data.m);
    ARRAY<T> u_linf(sim_data.m);
    ARRAY<T> u_2(sim_data.m);
    CUBIC_MN_INTERPOLATION_UNIFORM<TV,T> interp;
    const GRID<TV>& base_grid(sim_data(0)->grid);
    const ARRAY<T,TV_INT>& base_phi=sim_data(0)->phi;
    for(CELL_ITERATOR<TV> it(base_grid,-1);it.Valid();it.Next()){
        int col=base_phi(it.index)>0;
        for(int d=0;d<TV::m;d++){
            ARRAY<T> values;
            for(int i=0;i<sim_data.m;i++){
                values.Append(interp.Periodic(sim_data(i)->face_grids(d),sim_data(i)->face_velocities(col).Component(d),it.Location()));}
            T exact=Approx_Exact(res,values,2);
            for(int i=0;i<sim_data.m;i++){
                T v=abs(values(i)-exact);
                samples(i)++;
                u_linf(i)=std::max(u_linf(i),v);
                u_2(i)+=sqr(v);}}}

    T Sx=0,Sy=0,Sz=0,Sxx=0,Sxy=0,Sxz=0;
    for(int i=0;i<sim_data.m;i++){
        T u2=sqrt(u_2(i)/samples(i)),uinf=u_linf(i),x=log10((T)res(i)),y=log10(u2),z=log10(uinf);
        Sx+=x;
        Sy+=y;
        Sz+=z;
        Sxx+=x*x;
        Sxy+=x*y;
        Sxz+=x*z;
        printf("%g %g\n", u2, uinf);}

    printf("order %g %g\n", -(Sy*Sx-sim_data.m*Sxy)/(-Sxx*sim_data.m+sqr(Sx)), -(Sz*Sx-sim_data.m*Sxz)/(-Sxx*sim_data.m+sqr(Sx)));

    LOG::cout<<"DIRECT ORDER TEST"<<std::endl;

    for(int i=0;i<res.m;i++){
        int a=res.Find(res(i)*2),b=res.Find(res(i)*4);
        if(a<0 || b<0) continue;
        T L2a=0,Linfa=0,L2b=0,Linfb=0;
        int cnt=0;
        for(CELL_ITERATOR<TV> it(sim_data(i)->grid,-1);it.Valid();it.Next()){
            int col=sim_data(i)->levelset.Phi(it.Location())>0;
            if(col!=(sim_data(a)->levelset.Phi(it.Location())>0)) continue;
            if(col!=(sim_data(b)->levelset.Phi(it.Location())>0)) continue;
            for(int d=0;d<TV::m;d++){
                T e0=interp.Periodic(sim_data(i)->face_grids(d),sim_data(i)->face_velocities(col).Component(d),it.Location());
                T e1=interp.Periodic(sim_data(a)->face_grids(d),sim_data(a)->face_velocities(col).Component(d),it.Location());
                T e2=interp.Periodic(sim_data(b)->face_grids(d),sim_data(b)->face_velocities(col).Component(d),it.Location());
                T f0=abs(e1-e0),f1=abs(e2-e1);
                cnt++;
                L2a+=sqr(f0);
                L2b+=sqr(f1);
                Linfa=max(Linfa,f0);
                Linfb=max(Linfb,f1);}}
        if(cnt)
        {
            L2a=sqrt(L2a/cnt);
            L2b=sqrt(L2b/cnt);
            printf("ORDER %3d %.3f %.3f\n", res(i), log(Linfa/Linfb)/log(2), log(L2a/L2b)/log(2));
        }
    }

    LOG::cout<<"DIRECT ORDER PRESSURE TEST"<<std::endl;

    for(int i=0;i<res.m;i++){
        int a=res.Find(res(i)*2),b=res.Find(res(i)*4);
        if(a<0 || b<0) continue;
        T L2a=0,Linfa=0,L2b=0,Linfb=0;
        int cnt=0;
        for(CELL_ITERATOR<TV> it(sim_data(i)->grid,-1);it.Valid();it.Next()){
            int col=sim_data(i)->levelset.Phi(it.Location())>0;
            if(col!=(sim_data(a)->levelset.Phi(it.Location())>0)) continue;
            if(col!=(sim_data(b)->levelset.Phi(it.Location())>0)) continue;
            T e0=interp.Periodic(sim_data(i)->grid,sim_data(i)->pressure,it.Location());
            T e1=interp.Periodic(sim_data(a)->grid,sim_data(a)->pressure,it.Location());
            T e2=interp.Periodic(sim_data(b)->grid,sim_data(b)->pressure,it.Location());
            T f0=abs(e1-e0),f1=abs(e2-e1);
            cnt++;
            L2a+=sqr(f0);
            L2b+=sqr(f1);
            Linfa=max(Linfa,f0);
            Linfb=max(Linfb,f1);}
        if(cnt)
        {
            L2a=sqrt(L2a/cnt);
            L2b=sqrt(L2b/cnt);
            printf("ORDER %3d %.3f %.3f\n", res(i), log(Linfa/Linfb)/log(2), log(L2a/L2b)/log(2));
        }
    }

    if(X0.Max()<FLT_MAX){
        LOG::cout<<"SAMPLES"<<std::endl;

        LOG::cout<<"sample at "<<X0<<std::endl;
        int col=0;
        for(int d=0;d<TV::m;d++){
            ARRAY<T> values;
            for(int i=0;i<sim_data.m;i++) values.Append(interp.Periodic(sim_data(i)->face_grids(d),sim_data(i)->face_velocities(col).Component(d),X0));
            LOG::cout<<"u("<<d<<")  "<<values<<"    "<<Approx_Exact(res,values,2)<<std::endl;}

        ARRAY<T> values,values2;
        for(int i=0;i<sim_data.m;i++) values.Append(interp.Periodic(sim_data(i)->grid,sim_data(i)->levelset.phi,X0));
        LOG::cout<<"level set  "<<values<<"    "<<Approx_Exact(res,values,2)<<std::endl;

        for(int i=0;i<sim_data.m;i++) values2.Append(interp.Periodic(sim_data(i)->grid,sim_data(i)->pressure,X0));
        LOG::cout<<"pressure  "<<values2<<"    "<<Approx_Exact(res,values2,2)<<std::endl;}
}

int main(int argc,char *argv[])
{
    bool use_3d=false,use_double=false;
    PARSE_ARGS parse_args(argc,argv);
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;
    parse_args.Add("-3d",&use_3d,"analyze 3D data");
    parse_args.Add("-double",&use_double,"analyze 3D data");
    parse_args.Parse(true);

    if(use_3d){
        if(use_double) Analyze<VECTOR<double,3> >(parse_args);
        else Analyze<VECTOR<float,3> >(parse_args);}
    else{
        if(use_double) Analyze<VECTOR<double,2> >(parse_args);
        else Analyze<VECTOR<float,2> >(parse_args);}

    LOG::Finish_Logging();
    return 0;
}
