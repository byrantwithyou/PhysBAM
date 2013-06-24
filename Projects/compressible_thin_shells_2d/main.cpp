#include <Tools/Parsing/PARSE_ARGS.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include "Strawman_Example/STRAWMAN_DRIVER.h"
#include "Strawman_Example/STRAWMAN_EXAMPLE.h"

using namespace PhysBAM;
int main(int argc,char* argv[]){
    bool run_strawman_example=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-strawman",&run_strawman_example,"Use strawman test");
    parse_args.Parse();

    typedef float T;
    if(run_strawman_example){
        STRAWMAN_EXAMPLE<VECTOR<T,2> > example;
        PARSE_ARGS parse_args(argc,argv);
        example.Parse(parse_args);
        STRAWMAN_DRIVER<VECTOR<T,2> > driver(example);
        example.Initialize();
        driver.Execute_Main_Program();}
    else
        LOG::cout<<"No example specified!"<<std::endl;
}
