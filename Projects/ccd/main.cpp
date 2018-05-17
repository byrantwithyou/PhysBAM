#include <algorithm>
#include <ctime>
#include <set>
#include <functional>
#include <vector>
#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include "DOUBLE_ERROR.h"

using namespace PhysBAM;

typedef double EXACT;
typedef VECTOR<EXACT,3> TV;
typedef VECTOR<DOUBLE_ERROR,3> ERROR_POINT;
typedef VECTOR<EXACT,3> EXACT_POINT;
typedef TRIANGLE_3D<EXACT> EXACT_TRIANGLE;

struct POINT_ERROR{
    ERROR_POINT parent;     // Center of the range 'point' was sampled from
    ERROR_POINT point;      // Sampled point
    DOUBLE_ERROR distance;  // Distance from point to triangle

    POINT_ERROR(const ERROR_POINT& a,const ERROR_POINT& b,DOUBLE_ERROR c):parent(a),point(b),distance(c){}
};

unsigned int SEED=time(NULL);
int SAMPLES=10000; // Number of poins sampled within a given range
int CUT_TO_TOP=100; // Explore the top 'CUT_TO_TOP' POINT_ERRORs
int OPTION=1; // 1-Minimize Distance 2-Maximize Error
EXACT lower_bound=-1;
EXACT upper_bound=1;

RANDOM_NUMBERS<double> rng(SEED);
POINT_ERROR minimum_distance(ERROR_POINT(),ERROR_POINT(),DOUBLE_ERROR(DBL_MAX));
POINT_ERROR maximum_error(ERROR_POINT(),ERROR_POINT(),DOUBLE_ERROR(DBL_MIN));
POINT_ERROR maximum_error_failure(ERROR_POINT(),ERROR_POINT(),DOUBLE_ERROR(DBL_MIN));

//############################################################################
// Function - Generate_Triangle
// Generates an exact random triangle with vertices in the range [lb,ub] and a
// distance of at most one from the origin
//#############################################################################
void Generate_Triangle(EXACT_TRIANGLE& t)
{
    for(EXACT_POINT& vertex:t.X){
        TV direction;
        rng.Fill_Uniform(direction,-1,1);
        direction.Normalize();
        vertex=direction*rng.Get_Uniform_Number(0,upper_bound-lower_bound);}
}

//####################################
// Function - Compute_VV
//####################################        
DOUBLE_ERROR Compute_VV(const ERROR_POINT& p,const EXACT_POINT& v)
{
    ERROR_POINT p1(v);
    return sqrt((p-p1).Magnitude_Squared());
}

//####################################
// Function - Compute_VE
//####################################        
DOUBLE_ERROR Compute_VE(const ERROR_POINT& p,const EXACT_POINT& v1,const EXACT_POINT& v2)
{
    ERROR_POINT p1(v1),p2(v2);
    return abs((p2-p1).Cross(p-p1).Magnitude()/(p2-p1).Magnitude());
}

//####################################
// Function - Compute_VF
//####################################        
DOUBLE_ERROR Compute_VF(const ERROR_POINT& p,const EXACT_TRIANGLE& t)
{
    ERROR_POINT p0(t.X[0]),p1(t.X[1]),p2(t.X[2]);
    ERROR_POINT temp=(p1-p0).Cross(p2-p0);
    return abs(temp.Dot(p-p0))/temp.Magnitude();
}

//####################################
// Function - Compute_Collision
//####################################
DOUBLE_ERROR Compute_Collision(const ERROR_POINT& p,const EXACT_TRIANGLE& t)
{
    std::function<DOUBLE_ERROR(const DOUBLE_ERROR&,const DOUBLE_ERROR&)> obj_fun=(OPTION==1) ? 
        [](const DOUBLE_ERROR& e1,const DOUBLE_ERROR& e2)->DOUBLE_ERROR{return std::min(e1,e2);}: 
        [](const DOUBLE_ERROR& e1,const DOUBLE_ERROR& e2)->DOUBLE_ERROR{return std::max(e1.error(),e2.error());};

    DOUBLE_ERROR obj = DOUBLE_ERROR(DBL_MAX);
    for(EXACT_POINT vertex:t.X)
        obj=obj_fun(obj,Compute_VV(p,vertex));
    for(int i=0;i<3;i++)
        obj=obj_fun(obj,Compute_VE(p,t.X[i%3],t.X[(i+1)%3]));
    obj=obj_fun(obj,Compute_VF(p,t));
    return obj;
}

//##############################################################################
// Function - Generate_Range
// Given a center point and side-length 'delta',return an EXACT range to
// sample points from.
//#############################################################################
void Generate_Range(const ERROR_POINT p,EXACT delta,RANGE<EXACT_POINT>& r)
{
    EXACT_POINT center(p[0].dbl,p[1].dbl,p[2].dbl);
    r=RANGE<EXACT_POINT>(center-delta/2,center+delta/2); 
}

//##############################################################################
// Function - Sample_Range
// 
// ############################################################################
void Sample_Range(std::vector<POINT_ERROR>& new_points,const std::vector<POINT_ERROR>& old_points,const EXACT_TRIANGLE& t,EXACT delta)
{
    //TODO: Allow options to change how points are prioritized
    std::set<POINT_ERROR,std::function<bool(const POINT_ERROR&,const POINT_ERROR&)>> point_errors(
            [] (const POINT_ERROR& l,const POINT_ERROR& r){return l.distance<r.distance;});

    for(POINT_ERROR point_error:old_points){
        RANGE<VECTOR<EXACT,3>> r;
        Generate_Range(point_error.point,delta,r);
        for(int i=1;i<=SAMPLES;i++){
            EXACT_POINT sample_point;
            rng.Fill_Uniform(sample_point,r);
            ERROR_POINT p(sample_point);
            for(int i=0;i<3;i++)
                p[i].qd=p[i].dbl+((point_error.point[i].dbl<point_error.point[i].qd)?1:-1)*fabsq(point_error.point[i].dbl-point_error.point[i].qd);
            POINT_ERROR temp(point_error.point,p,Compute_Collision(p,t));
            point_errors.insert(temp);
            maximum_error=(temp.distance.error() > maximum_error.distance.error()) ? temp : maximum_error;
            minimum_distance=(temp.distance<minimum_distance.distance)?temp:minimum_distance;
            if(temp.distance.dbl<temp.distance.error())
                maximum_error_failure=(temp.distance.error()>maximum_error_failure.distance.error())?temp:maximum_error_failure;}}

    int i=0;
    for(auto it=point_errors.begin(); it!=point_errors.end()&&i<CUT_TO_TOP;it++,i++)
        new_points.push_back(*it);
}

void Write_Experiment_Header()
{
    // Output Experiment Information
    LOG::printf("-----Experiment Information-----\n");
    char fmt[64]="  %-25s%P\n";
    LOG::printf(fmt,"Seed:",SEED);
    LOG::printf(fmt,"SAMPLES:",SAMPLES);
    LOG::printf(fmt,"CUT TO TOP:",CUT_TO_TOP);
    LOG::printf("--------------------------------\n\n");
}

void Write_Point_Information(std::string iteration,double length,const POINT_ERROR& pe)
{
    LOG::printf("Iteration %P Bounding Box Length: %P:\n\tParent: %P\n\tPoint: %P\n\tDistance: %0.20P\n\tError: %0.20P\n",
            iteration,length,pe.parent,pe.point,pe.distance.dbl,pe.distance.error());
}

int main(int argc,char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-Samples",&SAMPLES,"SAMPLES","Number of points sampled");
    parse_args.Add("-Cut_To_Top",&CUT_TO_TOP,"CUT_TO_TOP","Search around the first N sample points");
    parse_args.Add("-Sort_Option",&OPTION,"OPTION","1-Minimum Distance  2-Maximum Error");
    parse_args.Add("-Lower_bound",&lower_bound,"LOWER_BOUND","Bottom-left corner of triangle vertex bounding box");
    parse_args.Add("-Upper_bound",&upper_bound,"UPPER_BOUND","Top-right corner of triangle vertex bounding box");
    parse_args.Parse();
    Write_Experiment_Header();

    // Generate Base Triangle
    EXACT_TRIANGLE baseTriangle;
    Generate_Triangle(baseTriangle);
    LOG::printf("Base Triangle: %P\n",baseTriangle);

    std::vector<POINT_ERROR> points;
    points.push_back(POINT_ERROR(ERROR_POINT(0,0,0),ERROR_POINT(0,0,0),0));

    unsigned int iteration=0;
    EXACT accuracy=upper_bound-lower_bound;
    while(accuracy>points[0].distance.error()){
        std::vector<POINT_ERROR> tempPoints;
        Sample_Range(tempPoints,points,baseTriangle,accuracy);
        points=tempPoints;
        accuracy/=2;
        Write_Point_Information(std::to_string(++iteration),accuracy,points[0]);}

    Write_Point_Information("Minimum Distance",accuracy,minimum_distance);
    Write_Point_Information("Maximum Error",accuracy,maximum_error);
    Write_Point_Information("Maximum Error Failure",accuracy,maximum_error_failure);
    return 0;
}

