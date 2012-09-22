#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "Strawman_Example/STRAWMAN_EXAMPLE.h"

using namespace PhysBAM;
int main(int argc,char* argv[]){

    bool run_strawman_example=false;
    int resolution=100;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-strawman",&run_strawman_example,"Use strawman test");
    parse_args.Add("-resolution",&resolution,"res","grid resolution");
    parse_args.Parse();

    typedef float T;
    if(run_strawman_example){
        STRAWMAN_EXAMPLE<T,T> example(resolution);
        example.Initialize();example.Run();}
    else
        LOG::cout<<"No example specified!"<<std::endl;
}
