#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "Strawman_Example/STRAWMAN_DRIVER.h"
#include "Strawman_Example/STRAWMAN_EXAMPLE.h"

using namespace PhysBAM;
int main(int argc,char* argv[]){
    bool run_strawman_example=PARSE_ARGS::Find_And_Remove("-strawman",argc,argv);

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
