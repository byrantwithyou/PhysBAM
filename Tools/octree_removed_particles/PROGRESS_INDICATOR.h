//#####################################################################
// Copyright 2002-2005, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_BLENDER
//##################################################################### 
#ifndef __PROGRESS_INDICATOR__
#define __PROGRESS_INDICATOR__

namespace PhysBAM {

class PROGRESS_INDICATOR
{
public:
    int total;
    int done;
    int percent_done;

    PROGRESS_INDICATOR(int total_input) : total(total_input),done(0),percent_done(0)
    {}

    bool Progress(int by=1)
    {
        done+=by;
        int new_percent_done=100*done/total;
        if(new_percent_done>percent_done){
            percent_done=new_percent_done;
            std::cout<<percent_done<<"% "<<std::flush;
            if(percent_done==100) std::cout<<std::endl;
            return true;
        }
        return false;
    }
};

}
#endif
