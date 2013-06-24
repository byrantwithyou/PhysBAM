#include <Tools/Arrays/ARRAY.h>
#include <sstream>
namespace PhysBAM{

struct COMMAND_ARG
{
    COMMAND_ARG():type(None),data(0){}
    COMMAND_ARG(int type_input,const std::string& name_input,void* data_input)
        :type(type_input),name(name_input),data(data_input) {}
    enum { None,Bool,Int,Float,Double,String };
    int type;
    std::string name;
    void* data;
};

struct COMMAND_BASE
{
    ARRAY<std::string>* line;  //the input line is divided into args,inputs,outputs
    ARRAY<std::string> inputs;
    ARRAY<std::string> outputs;
    ARRAY<COMMAND_ARG> args;

    COMMAND_BASE():line(0){}
    virtual ~COMMAND_BASE(){}
    void Add_Bool_Arg(const std::string& name,bool& val);
    void Add_Int_Arg(const std::string& name,int& val);
    void Add_Float_Arg(const std::string& name,float& val);
    void Add_Double_Arg(const std::string& name,double& val);
    void Add_String_Arg(const std::string& name,std::string& val);
    std::string Get_Arg_Name(int arg) const { return args(arg).name; }

    virtual std::string Name()=0;
    virtual void Print_Description(std::ostream& out);
    virtual void Print_Detailed_Description(std::ostream& out);
    virtual std::string Description() { return "unknown function."; }
    virtual int Num_Args() { return args.m; }
    virtual int Min_Inputs() { return 0; } 
    virtual int Max_Inputs() { return 0; } //if < 0, possibly infinite }
    virtual int Num_Outputs() { return 0; } //can only be 0 or 1, for now
    virtual bool Process_Line();
    virtual bool Set_Arg(int arg,const std::string& val);
    virtual int Do() = 0;     //return > 0 if success, <=0 on failure
    
    std::string Get_Input(int i) const { return inputs(i); }
    std::string Get_Output(int i=1) { return outputs(i); }
    int Num_Inputs() const { return inputs.m; }
};



struct COMMAND_AUTO : public COMMAND_BASE
{
    const std::string name;
    const std::string desc;

    COMMAND_AUTO(const std::string& name_input,const std::string& desc_input="")
        :name(name_input),desc(desc_input){}
    virtual std::string Name() { return name; }
    virtual std::string Description() { if(!desc.empty()) return desc; else return COMMAND_BASE::Description(); }
};

//T is assumed to have istream >> defined
template <class T>
struct COMMAND_INPUT_VALUE : public COMMAND_AUTO
{
    T* target;

    COMMAND_INPUT_VALUE(const std::string& name_input,T* target_input)
        :COMMAND_AUTO(name,": sets a value."),target(target_input) {}
    virtual int Num_Args() { return 1; }
    virtual bool Set_Arg(int arg,const std::string& val) {
        std::istringstream str;
        str.str(val);
        str >> *target;
        if(!str)
            return false;
        return true;
    }
    virtual int Do() {
        return 1;
    }
};

template <class T>
struct COMMAND_SET_VALUE : public COMMAND_AUTO
{
    T* target;
    const T& val;

    COMMAND_SET_VALUE(const std::string& name_input,T* target_input,const T& val_input)
        :COMMAND_AUTO(name_input,": sets a value."),target(target_input),val(val_input) {}
    virtual int Do() {
        *target = val;
        return 1;
    }
};

void Command_Print_Usage(const std::string& appName);
void Command_Print(COMMAND_BASE* commands [],int numCommands);
int Command_Reader(COMMAND_BASE* commands [],int numCommands,int argc,const char** argv);

};
