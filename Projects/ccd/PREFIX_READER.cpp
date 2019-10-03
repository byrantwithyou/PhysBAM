#include <Core/Log/LOG.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <sstream>
#include <string>
#include "DOUBLE_ERROR.h"

using namespace PhysBAM;

enum INSTRUCTION {ADD,SUBTRACT,DIVIDE,MULTIPLY,OPERAND,OTHER};

INSTRUCTION Get_Token(char* buf, int sz, std::istringstream& iss)
{
    int idx=0;
    while(iss.peek()!='('&&iss.peek()!=')'&&iss.peek()!=','&&idx<sz-1) buf[idx++]=iss.get();
    while((iss.peek()=='('||iss.peek()==')'||iss.peek()==',')&&idx<sz-1) iss.get();
    
    if(idx==sz)LOG::printf("Critical Error: Buffer overflow.\n\tRead: %P\n", std::string(buf));
    buf[idx++]='\0';

    if(strcmp(buf,"ADD")==0) return ADD;
    if(strcmp(buf,"SUB")==0) return SUBTRACT;
    if(strcmp(buf,"DIV")==0) return DIVIDE;
    if(strcmp(buf,"MUL")==0) return MULTIPLY;

    try{std::stod(std::string(buf));return OPERAND;}
    catch(const std::invalid_argument& ia) {LOG::printf("CRITICAL ERROR: Could not identify token: [%P]\n", std::string(buf));}

    return OTHER;
}

DOUBLE_ERROR Calculate(std::istringstream& iss)
{
    const int SIZE=128;
    char buf[SIZE];
    buf[SIZE-1]='\0';
    INSTRUCTION instr = Get_Token(buf,SIZE,iss);
    DOUBLE_ERROR ret;
    switch(instr){
        case ADD:ret=Calculate(iss)+Calculate(iss);break;
        case SUBTRACT:ret=Calculate(iss)-Calculate(iss);break;
        case DIVIDE:ret=Calculate(iss)/Calculate(iss);break;
        case MULTIPLY:ret=Calculate(iss)*Calculate(iss);break;
        case OPERAND:ret=std::stod(buf);break;
        case OTHER:LOG::printf("Missed Operator: %P\n", buf);
    }
    LOG::printf("%P %P\n", qdtostr(ret.error()), ret.val);
    return ret;
}

int main(int argc,char* argv[])
{
    std::string INPUT="";
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-INPUT",&INPUT,"Input String","String to parse");
    parse_args.Parse();

    //TODO: Allow for parsing of each instruction in a file
    std::istringstream iss(INPUT);
    Calculate(iss);
    return 0;
}

