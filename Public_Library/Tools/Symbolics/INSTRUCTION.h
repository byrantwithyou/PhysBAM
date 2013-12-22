//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
// Class INSTRUCTION
//#####################################################################
#ifndef __INSTRUCTION__
#define __INSTRUCTION__

namespace PhysBAM{

enum op_type {
    op_nop,op_copy,op_add,op_sub,op_mul,op_div,op_mod,op_neg,op_inv,
    op_sqrt,op_exp,op_ln,op_pow,
    op_sin,op_cos,op_asin,op_acos,op_atan,op_atan2,
    op_lt,op_le,op_gt,op_ge,op_eq,op_ne,op_not,op_or,op_and,
    op_br_z,op_br_nz,op_jmp,op_label,op_phi,
    op_last
};

const int mem_shift=28;
const int mem_mask=0xF0000000;
const int mem_reg=0x00000000;
const int mem_const=0x10000000;
const int mem_in=0x20000000;
const int mem_out=0x30000000;

struct INSTRUCTION
{
    op_type type;
    int dest,src0,src1;
};
}
#endif
