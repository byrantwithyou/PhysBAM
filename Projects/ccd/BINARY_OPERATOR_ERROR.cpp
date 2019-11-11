#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include "DOUBLE_ERROR.h"
#include <algorithm>
#include <ctime>
#include <functional>
#include <string>

using namespace PhysBAM;

typedef double EXACT;
typedef DOUBLE_ERROR ERROR;

int main(int argc,char* argv[])
{
    unsigned int SEED=time(NULL);
    unsigned int SAMPLES=100000;
    EXACT LOWER_BOUND=-1.0;
    EXACT UPPER_BOUND=1.0;
    char OPERATION[]="ADD";
    bool DEBUG=false;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-SEED",&SEED,"SEED","Random number generation seed");
    parse_args.Add("-SAMPLES",&SAMPLES,"SAMPLES","Experiment iteraction count");
    parse_args.Add("-LOWERBOUND",&LOWER_BOUND,"Lower Bound","...");
    parse_args.Add("-UPPERBOUND",&UPPER_BOUND,"Upper Bound","...");
    parse_args.Add("-OPERATION",&OPERATION,"Operation","ADD SUB MUL DIV");
    parse_args.Add("-DEBUG:",&DEBUG,"Debug","Output results of each iteration");
    parse_args.Parse();
 
    LOG::printf("--------------------------------\n\n");
    char fmt[64]=" %-25s%P\n";
    LOG::printf(fmt,"Seed:",SEED);
    LOG::printf(fmt,"Samples:",SAMPLES);
    LOG::printf(fmt,"Lower Bound:",LOWER_BOUND);
    LOG::printf(fmt,"Upper Bound:",UPPER_BOUND);
    LOG::printf(fmt,"Operation:",OPERATION);
    LOG::printf("--------------------------------\n");

    std::function<DOUBLE_ERROR(const DOUBLE_ERROR&,const DOUBLE_ERROR&)> operation;

    RANDOM_NUMBERS<EXACT> rng(SEED);
    if(strcmp(OPERATION,"ADD")==0) operation=[](const DOUBLE_ERROR& e1,const DOUBLE_ERROR& e2)->DOUBLE_ERROR{return e1+e2;};
    if(strcmp(OPERATION,"SUB")==0) operation=[](const DOUBLE_ERROR& e1,const DOUBLE_ERROR& e2)->DOUBLE_ERROR{return e1-e2;};
    if(strcmp(OPERATION,"MUL")==0) operation=[](const DOUBLE_ERROR& e1,const DOUBLE_ERROR& e2)->DOUBLE_ERROR{return e1*e2;};
    if(strcmp(OPERATION,"DIV")==0) operation=[](const DOUBLE_ERROR& e1,const DOUBLE_ERROR& e2)->DOUBLE_ERROR{return e1/e2;};

    __float128 maximum_error = 0;
    ERROR result;
    for(int i=1;i<=SAMPLES;i++){
        //d1 = 0.30500603283022598156;
        //d2 = 0.03257104180524249271;
        ERROR e1(rng.Get_Uniform_Number(LOWER_BOUND,UPPER_BOUND));
        ERROR e2(rng.Get_Uniform_Number(LOWER_BOUND,UPPER_BOUND));
        result=operation(e1,e2);
        if(DEBUG)LOG::printf("%0.20P %0.20P %0.20P %P\n",e1.dbl,e2.dbl,result.dbl,qdtostr(result.error()));
        maximum_error=std::max(maximum_error,result.error());}

    LOG::printf("Maximum Error: %P\n", qdtostr(maximum_error));
    return 0;
}

