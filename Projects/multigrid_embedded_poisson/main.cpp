//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <cassert>

#include <iostream>
#include <string>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception/exception.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include <Jeffrey_Utilities/BASIC_TIMER.h>

#include "main_impl.h"

//#####################################################################
// main
//#####################################################################
int main(int argc, char* argv[])
{
    namespace bpo = boost::program_options;

    PhysBAM::BASIC_TIMER timer;

    // Parse the dimension and the scalar type to determine which instantiation
    // of main_impl to use.
    unsigned int dimension;
    std::string scalar_str;
    {
        bpo::options_description opts_desc("");
        opts_desc.add_options()
            ( "dimension", bpo::value< unsigned int >(&dimension)->default_value(3),        "" )
            ( "scalar"   , bpo::value< std::string >(&scalar_str)->default_value("double"), "" );
        bpo::variables_map var_map;
        try {
            bpo::store(
                bpo::command_line_parser(argc, argv).options(opts_desc).allow_unregistered().run(),
                var_map
            );
        }
        catch(boost::exception& e) {
            std::cerr << "ERROR: Exception thrown while parsing command line:\n"
                      << boost::diagnostic_information(e)
                      << "(Run with --help for usage.)"
                      << std::endl;
            return 1;
        }
        bpo::notify(var_map);
    }

    int main_result = 0;

    if(dimension != 2 && dimension != 3) {
        std::cerr << "ERROR: Unsigned int parameter --dimension must be 2 or 3.\n"
                  << "(Run with --help for usage.)"
                  << std::endl;
        main_result = 1;
    }
    if(scalar_str != "float" && scalar_str != "double") {
        std::cerr << "ERROR: Unable to parse string parameter --scalar : " << scalar_str << '\n'
                  << "(Run with --help for usage.)"
                  << std::endl;
        main_result = 1;
    }

    if(main_result == 0)
        switch(dimension) {
        case 2:
            if(scalar_str == "float")
                main_result = PhysBAM::Embedded_Poisson_V2::main_impl< float, 2 >(argc, argv);
            else if(scalar_str == "double")
                main_result = PhysBAM::Embedded_Poisson_V2::main_impl< double, 2 >(argc, argv);
            break;
        case 3:
            if(scalar_str == "float")
                main_result = PhysBAM::Embedded_Poisson_V2::main_impl< float, 3 >(argc, argv);
            else if(scalar_str == "double")
                main_result = PhysBAM::Embedded_Poisson_V2::main_impl< double, 3 >(argc, argv);
            break;
        default:;
        }

    std::cout << "[Total execution time: " << timer.Elapsed() << " s]"
              << std::endl;

    return main_result;
}
