#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include "DOUBLE_ERROR.h"
#include <ctime>
#include <functional>

using namespace PhysBAM;

typedef double EXACT;
typedef DOUBLE_ERROR ERROR;
typedef VECTOR<EXACT,2> TV;
typedef VECTOR<DOUBLE_ERROR,2> ERROR_POINT;
typedef VECTOR<EXACT,2> EXACT_POINT;
typedef TRIANGLE_2D<EXACT> EXACT_TRIANGLE;

struct POINT_ERROR
{
    ERROR_POINT point;  // Sampled point
    ERROR distance;     // Distance from point to triangle

    POINT_ERROR(const ERROR_POINT& b,ERROR c):point(b),distance(c){}
};

struct COLLISION_RESULT
{
    bool isCollision;
    VECTOR<ERROR,3> point;
    COLLISION_RESULT():isCollision(false),point(VECTOR<DOUBLE_ERROR,3>(0,0,0)){}
    COLLISION_RESULT(bool bl, PhysBAM::VECTOR<DOUBLE_ERROR,3> v) : isCollision(bl),point(v){}
    COLLISION_RESULT(const COLLISION_RESULT& c):isCollision(c.isCollision),point(c.point){}
};

EXACT TOLERANCE_S = std::sqrt(ERROR::EPSILON);
EXACT TOLERANCE_L = 2*(1+5*ERROR::EPSILON)/(1-7*TOLERANCE_S);
EXACT TOLERANCE_T = TOLERANCE_S*TOLERANCE_L;
EXACT TOLERANCE_K = 21*ERROR::EPSILON*std::pow(TOLERANCE_L,2);

RANDOM_NUMBERS<EXACT> rng;;

//##############################################################################
// Function - Generate_Triangle
// Generates an exact random triangle with vertices in the range [lb,ub] and a
// distance of at most one from the origin
//##############################################################################
void Generate_Triangle(EXACT_TRIANGLE& t, EXACT magnitude)
{
    for(EXACT_POINT& vertex:t.X){
        TV direction;
        rng.Fill_Uniform(direction,-1,1);
        direction.Normalize();
        vertex=direction*rng.Get_Uniform_Number(0,magnitude);}
}

//##############################################################################
// Function - Generate_Point_On_Triangle
// Using exact barycentric coordinates, produces a random point on the triangle
//#############################################################################
void Generate_Point_On_Triangle(EXACT_POINT& p, const EXACT_TRIANGLE& t)
{
    VECTOR<EXACT,3> b;
    b[0]=rng.Get_Uniform_Number(0,1);
    b[1]=rng.Get_Uniform_Number(0,1-b[0]);
    b[2]=1-b[0]-b[1];
    p=t.Point_From_Barycentric_Coordinates(b);
}
/*
//####################################
// Function - Compute_VV
//####################################        
ERROR Compute_VV(const ERROR_POINT& p,const EXACT_POINT& v)
{
    ERROR_POINT p1(v);
    return sqrt((p-p1).Magnitude_Squared());
}

//####################################
// Function - Compute_VE
//####################################        
ERROR Compute_VE(const ERROR_POINT& p,const EXACT_POINT& v1,const EXACT_POINT& v2)
{
    ERROR_POINT p1(v1),p2(v2);
    return abs((p2-p1).Cross(p-p1).Magnitude()/(p2-p1).Magnitude());
}
*/
//####################################
// Function - Compute_VF
//####################################        
void Compute_VF(COLLISION_RESULT&c, const ERROR_POINT& p,const TRIANGLE_2D<EXACT>& t)
{
    ERROR_POINT p0(t.X[0]),p1(t.X[1]),p2(t.X[2]);
    ERROR area_a=((p1-p).Cross(p2-p)).Sum();
    ERROR area_b=((p-p0).Cross(p2-p0)).Sum();
    ERROR area_c=((p1-p0).Cross(p-p0)).Sum();
    if(abs(area_a)<TOLERANCE_K||abs(area_b)<TOLERANCE_K||abs(area_c)<TOLERANCE_K||
            ((area_a<0)!=(area_b<0))||((area_a<0)!=(area_c<0))||((area_b<0)!=(area_c<0))) c=COLLISION_RESULT(false,VECTOR<ERROR,3>(0,0,0));
    else{
        ERROR sum=area_a+area_b+area_c;
        c=COLLISION_RESULT(true,VECTOR<DOUBLE_ERROR,3>(area_a/sum,area_b/sum,area_c/sum));}
}

int main(int argc,char* argv[])
{
    unsigned int SEED=time(NULL);
    unsigned int SAMPLES=100000;
    EXACT RADIUS=1.0;
    bool PRINT_EXPRESSIONS=false;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-SEED",&SEED,"Seed","...");
    parse_args.Add("-SAMPLES",&SAMPLES,"Samples","Experiment iterations");
    parse_args.Add("-RADIUS",&RADIUS,"Radius","Maximum distance of each point from the origin.");
    parse_args.Add("-PRINT:",&PRINT_EXPRESSIONS,"Print Expression","Outputs prefix notation for each barycentric point calculation");
    parse_args.Parse();

    LOG::printf("-----Experiment Information-----\n");
    char fmt[64]="  %-25s%P\n";
    LOG::printf(fmt,"SEED:",SEED);
    LOG::printf(fmt,"SAMPLES:",SAMPLES);
    LOG::printf(fmt,"RADIUS:",RADIUS);
    LOG::printf("--------------------------------\n\n");

    rng.Set_Seed(SEED);
    // Generate Base Triangle
    EXACT_TRIANGLE baseTriangle;
    Generate_Triangle(baseTriangle,RADIUS);
    LOG::printf("Base Triangle: %P\n",baseTriangle);

    quad maximum_error=0;
    for(int i=0;i<SAMPLES;i++){
        EXACT_POINT p;
        Generate_Point_On_Triangle(p,baseTriangle);
        ERROR_POINT point(p);
        COLLISION_RESULT c;
        Compute_VF(c,point,baseTriangle);
        if(PRINT_EXPRESSIONS)for(ERROR e:c.point)LOG::printf("%P\n",e.val);
        if(c.isCollision)
            for(ERROR e:c.point)
                maximum_error=std::max(maximum_error,e.error());}
    LOG::printf("Maximum Error: %P\n",qdtostr(maximum_error));

    return 0;
}
