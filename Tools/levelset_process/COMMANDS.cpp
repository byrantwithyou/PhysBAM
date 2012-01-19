#include <cstring>
#include "COMMANDS.h"
using namespace PhysBAM;
using namespace std;

namespace PhysBAM {

bool string_is_bool(const std::string& arg,bool& val)
{
    if(arg=="0") val = false;
    else if(arg=="1") val = true;
    else if(arg=="true") val = true;
    else if(arg=="false") val = false;
    else return false;
    return true;
}

bool string_is_int(const std::string& arg, int& val)
{
    for(std::string::size_type i=0;i<arg.length();i++)if(!(isdigit(arg[i])||arg[i]=='-')) return false;
    val=atoi(arg.c_str());
    return true;
}

bool string_is_float(const std::string& arg,float& val)
{
    for(std::string::size_type i=0;i<arg.length();i++)if(!(isdigit(arg[i])||arg[i]=='.'||arg[i]=='-')) return false;
    val=(float)atof(arg.c_str());
    return true;
}

bool string_is_double(const std::string& arg,double& val)
{
    for(std::string::size_type i=0;i<arg.length();i++)if(!(isdigit(arg[i])||arg[i]=='.'||arg[i]=='-')) return false;
    val=atof(arg.c_str());
    return true;
}

const int Bool = 1;
const int Int = 2;
const int Float = 3;
const int Double = 4;
const int String = 5;

void COMMAND_BASE::Print_Description(std::ostream& out)
{
    out<<"-"<<Name();
    for(int i=0;i<args.m;i++) {
        out << " " << args(i).name;
    }
    out<<": "<<Description()<<endl;
}

void COMMAND_BASE::Print_Detailed_Description(ostream& out)
{
    out<<"-"<<Name()<<":"<<endl;
    out<<"Args:"<<endl;
    for(int i=0;i<args.m;i++) {
        switch(args(i).type) {
        case Bool:
            out<<"\tbool "<<args(i).name<<" (default "<<(*(bool*)(args(i).data) ? "true":"false")<<")"<<endl;
            break;
        case Int:
            out<<"\tint "<<args(i).name<<" (default "<<(*(int*)(args(i).data))<<")"<<endl;
            break;
        case Float:
            out<<"\tfloat "<<args(i).name<<" (default "<<(*(float*)(args(i).data))<<")"<<endl;
            break;
        case Double:
            out<<"\tdouble "<<args(i).name<<" (default "<<(*(double*)(args(i).data))<<")"<<endl;
            break;
        case String:
            out<<"\tstring "<<args(i).name<<" (default "<<(*(std::string*)(args(i).data))<<")"<<endl;
            break;
        default:
            out<<"Error!"<<endl;
            break;
        }
    }
    if(Max_Inputs() < 0)
        out<<"Inputs: at least "<<Min_Inputs()<<endl;
    else if(Max_Inputs() == Min_Inputs())
        out<<"Inputs: "<<Min_Inputs()<<endl;
    else
        out<<"Inputs: between "<<Min_Inputs()<<" and "<<Max_Inputs()<<endl;
    out<<"Outputs: "<<Num_Outputs()<<endl;
    out<<"Description: "; Print_Description(out);
}
void COMMAND_BASE::Add_Bool_Arg(const std::string& name,bool& val) {
    args.Append(COMMAND_ARG(Bool,name,&val)); }
void COMMAND_BASE::Add_Int_Arg(const std::string& name,int& val) {
    args.Append(COMMAND_ARG(Int,name,&val)); }
void COMMAND_BASE::Add_Float_Arg(const std::string& name,float& val) {
    args.Append(COMMAND_ARG(Float,name,&val)); }
void COMMAND_BASE::Add_Double_Arg(const std::string& name,double& val) {
    args.Append(COMMAND_ARG(Double,name,&val)); }
void COMMAND_BASE::Add_String_Arg(const std::string& name,std::string& val) {
    args.Append(COMMAND_ARG(String,name,&val)); }

bool COMMAND_BASE::Process_Line()
{
    assert(line != NULL);
    if(line->m < Num_Args()) {
        cout<<"Not enough arguments to "<<Name()<<endl;
        return false;
    }
    int args_remaining = line->m;
    int j;
    for(j=1; j<=Num_Args(); j++,args_remaining--) {
        if(!Set_Arg(j,(*line)(j))) {
            cout<<"Error reading argument "<<Get_Arg_Name(j)<<"="<<(*line)(j)<<" to "<<Name()<<endl;
            return false;
        }
    }
    int num_inputs = args_remaining - Num_Outputs();
    if(Max_Inputs() >= 0) {
        if(num_inputs > Max_Inputs()) {
            cout << "The number of max inputs to -" <<Name()<<" was exceeded."<<endl;
            return false;
        }
    }
    if(num_inputs < Min_Inputs()) {
        cout << "Not enough inputs given to -" <<Name()<<endl;
        return false;
    }
    inputs.Resize(num_inputs);
    for(int i=1;i<=num_inputs;i++,j++,args_remaining--) {
        inputs(i)=(*line)(j);
    }
    assert(args_remaining == Num_Outputs());
    outputs.Resize(Num_Outputs());
    for(int i=1;i<=Num_Outputs();i++,j++,args_remaining--) {
        outputs(i)=(*line)(j);
    }
    assert(j == line->m+1);
    assert(args_remaining == 0);
    return true;
}

bool COMMAND_BASE::Set_Arg(int arg,const std::string& val)
{
    if(val=="#") return true;  //default arguments
    switch(args(arg).type) {
    case Bool: return string_is_bool(val,*(bool*)args(arg).data);
    case Int: return string_is_int(val,*(int*)args(arg).data);
    case Float: return string_is_float(val,*(float*)args(arg).data);
    case Double: return string_is_double(val,*(double*)args(arg).data);
    case String: *(std::string*)args(arg).data = val;
        return true;
    }
    return false;
}

void Command_Print_Usage(const std::string& appName)
{
    cout << "USAGE: " << appName << " [COMMAND] args inputs outputs ..." << endl;
    cout << "# as an argument denotes using the default argument." << endl;
}

void Command_Print(COMMAND_BASE* commands [], int num_commands)
{
    cout << "COMMANDS:"<<endl;
    for(int i=0; i<num_commands; i++)
        commands[i]->Print_Description(cout);
    cout<<"For more information, type -help commandname"<<endl;
}

struct COMMAND_HELP:public COMMAND_BASE
{
    COMMAND_BASE** commands;
    int num_commands;
    std::string item_name;

    COMMAND_HELP(COMMAND_BASE** commands_input, int num_commands_input)
        :commands(commands_input), num_commands(num_commands_input)
    {
        Add_String_Arg("command",item_name);
    }
    virtual std::string Name() { return "help"; }
    virtual std::string Description() { return "returns help on a command."; }
    virtual int Do() {
        bool found_command = false;
        for(int k=0; k<num_commands && !found_command; k++) {
            if(commands[k]->Name()==item_name) {
                found_command = true;
                commands[k]->Print_Detailed_Description(cout);
            }
        }
        if(!found_command) {
            cout<<"Unknown command given to help :"<<item_name<<endl;
            return -1;
        }
        return 0;
    }
};

int Do_Command(COMMAND_BASE* c, ARRAY<std::string>& inputs)
{
    c->line = &inputs;
    if(!c->Process_Line()) return -1;
    c->line = NULL;
    int res = c->Do();
    if(res <= 0) return res;
    inputs.Resize(1);
    return 1;
}

int Command_Reader(COMMAND_BASE* commands [],int num_commands,int argc,const char** argv)
{
    ARRAY<std::string> inputs;
    COMMAND_BASE* last_command = NULL;
    COMMAND_HELP help_command(commands,num_commands);

    int command_count = 0;
    bool command_chosen = false;
    for(int i=1; i<argc; i++) {
        double tmp;
        if(argv[i][0] == '-' && !string_is_double(argv[i],tmp)) { //command
            if(string_is_double(argv[i],tmp))
            command_chosen = false;
            //check for help queries first
            if(0 == strcmp("help",&argv[i][1])) {
                if(last_command != NULL) {
                    int res = Do_Command(last_command,inputs);
                    if(res <= 0) return res;
                }
                last_command = &help_command;
                command_count++;
                command_chosen = true;
            }
            //search through rest of commands
            for(int k=0; k<num_commands && !command_chosen; k++) {
                if(commands[k]->Name()==&argv[i][1]) {
                    //we've found a command
                    if(last_command != NULL) {
                        int res = Do_Command(last_command,inputs);
                        if(res <= 0) return res;
                    }
                    last_command = commands[k];
                    command_chosen = true;
                    command_count++;
                }
            }
            if(!command_chosen) {
                printf("Unknown command %s\n", argv[i]);
                Command_Print(commands, num_commands);
                return -1;
            }
        }
        else {  //add an argument
            inputs.Append(argv[i]);
            command_chosen = false;
        }
    }
    if(command_count == 0) {
        if(inputs.m > 0)
            cout<<"Some stray elements were given on the command line"<<endl;
        Command_Print(commands, num_commands);
    }
    else {
        assert(last_command != NULL);
        int res = Do_Command(last_command,inputs);
        if(res <= 0) return res;
    }
    return 0;
}

};
