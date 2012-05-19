//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Symbolics/PROGRAM.h>

using namespace PhysBAM;

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    typedef double T;
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-expr","","expression to parse");
    parse_args.Add_Double_Argument("-x",5,"x");
    parse_args.Add_Double_Argument("-y",5,"y");
    parse_args.Add_Double_Argument("-z",5,"z");
    parse_args.Parse(argc,argv);
    std::string expr=parse_args.Get_String_Value("-expr");

    PROGRAM<double> prog;
    prog.var_in.Append("x");
    prog.var_in.Append("y");
    prog.var_in.Append("z");
    prog.var_out.Append("r");
    prog.var_out.Append("t");
    prog.var_out.Append("y");
    prog.Parse(expr.c_str());
    prog.Print();
    prog.Finalize();
    prog.Print();
    PROGRAM_CONTEXT<T> context(prog);
    context.data_in(0)=parse_args.Get_Double_Value("-x");
    context.data_in(1)=parse_args.Get_Double_Value("-y");
    context.data_in(2)=parse_args.Get_Double_Value("-z");
    LOG::cout<<prog.var_in<<std::endl;
    LOG::cout<<prog.var_out<<std::endl;
    LOG::cout<<context.data_in<<std::endl;
    prog.Execute(context);
    LOG::cout<<context.data_in<<std::endl;
    LOG::cout<<context.data_out<<std::endl;

    return 0;
}
