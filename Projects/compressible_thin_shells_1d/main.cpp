#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "Strawman_Example/STRAWMAN_EXAMPLE.h"

using namespace PhysBAM;
int main(int argc,char* argv[]){
    bool run_strawman_example=PARSE_ARGS::Find_And_Remove("-strawman",argc,argv);
    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-resolution",100);
    parse_args.Parse(argc,argv);

    typedef float T;
    if(run_strawman_example){
        STRAWMAN_EXAMPLE<T,T> example(parse_args.Get_Integer_Value("-resolution"));
        example.Initialize();example.Run();}
    else
        LOG::cout<<"No example specified!"<<std::endl;
}
