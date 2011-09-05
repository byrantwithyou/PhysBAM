// Parts borrowed from http://www.math.utah.edu/~beebe/software/ieee/ofl.c
#include <iostream>
#include <assert.h>
#include <fenv.h>
#include <float.h>
#include <fpu_control.h>
#include <limits.h>
#include <math.h>
#include <signal.h>

// TWO MODES:
// TRAP_EXCEPTIONS -- use signal handler to catch (first) floating point exception and exit
// !TRAP_EXCEPTIONS -- print floating point flags after each operation
//#define TRAP_EXCEPTIONS

struct sigaction* old_action=0;

typedef unsigned long fpflag_t;

#define FP_INVALID_OPERATION        0x0001
#define FP_DENORMALIZED_OPERAND        0x0002
#define FP_ZERO_DIVIDE            0x0004
#define FP_OVERFLOW            0x0008
#define FP_UNDERFLOW            0x0010
#define FP_PRECISION            0x0020
#define FP_ALL_EXCEPTIONS        0x003f /* OR of all the above */

void fpe_callback(int sig_number,siginfo_t *info,void *data)
{
    assert(sig_number==SIGFPE);
    std::cerr << "Floating point exception: reason " << info->si_code << " = \"" <<
        (info->si_code==FPE_INTDIV?"integer divide by zero":info->si_code==FPE_INTOVF?"integer overflow":
        info->si_code==FPE_FLTDIV?"FP divide by zero":info->si_code==FPE_FLTOVF?"FP overflow":
        info->si_code==FPE_FLTUND?"FP underflow":info->si_code==FPE_FLTRES?"FP inexact result":
        info->si_code==FPE_FLTINV?"FP invalid operation":info->si_code==FPE_FLTSUB?"subscript out of range":"unknown")
        << "\", from address 0x" << std::hex << (unsigned int)info->si_addr << std::endl;

    // seems we have no choice but to exit here because after a trapped FPE we can't reliably continue execution (e.g. see
    // bottom of sigaction man page)
    exit(1);
}

void clear_fp_flags()
{
    __asm __volatile ("fclex"); /* clear exception flags to avoid interrupt on next fldcw! */
}

fpflag_t get_fp_flags()
{
    static short int flags;
    __asm __volatile ("fstsw %0" : : "m" (flags));
    /* Return only the exception flag bits, ignoring all others */
    return ((fpflag_t)(flags & FP_ALL_EXCEPTIONS));
}

void print_fp_flags(fpflag_t mask)
{
    if (mask != (fpflag_t)0)
    {
        std::cerr << "  [";
        if (mask & FP_INVALID_OPERATION)    std::cerr << " FP_INVALID_OPERATION ";
        if (mask & FP_DENORMALIZED_OPERAND) std::cerr << " FP_DENORMALIZED_OPERAND ";
        if (mask & FP_ZERO_DIVIDE)          std::cerr << " FP_ZERO_DIVIDE ";
        if (mask & FP_OVERFLOW)             std::cerr << " FP_OVERFLOW ";
        if (mask & FP_UNDERFLOW)            std::cerr << " FP_UNDERFLOW ";
        if (mask & FP_PRECISION)            std::cerr << " FP_PRECISION ";
        std::cerr << "]";
    }
}

void print_and_clear()
{
#ifndef TRAP_EXCEPTIONS
    print_fp_flags(get_fp_flags());
    clear_fp_flags();
#endif
    std::cerr << std::endl;
}

int main()
{
#ifdef TRAP_EXCEPTIONS
    old_action=new struct sigaction;
    memset(old_action,0,sizeof(*old_action));
    sigaction(SIGFPE,0,old_action);

    struct sigaction action;
    memset(&action,0,sizeof(action));
    action.sa_sigaction=fpe_callback;
    sigemptyset(&action.sa_mask);
    action.sa_flags=SA_SIGINFO;
    if(sigaction(SIGFPE,&action,0)) std::cerr << "Could not register FPE signal handler" << std::endl;

    fpu_control_t foo= _FPU_DEFAULT &~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
    _FPU_SETCW(foo);
    feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

    float inf=std::numeric_limits<float>::infinity();
    float max=std::numeric_limits<float>::max();
    float nan=std::numeric_limits<float>::quiet_NaN();
    float one=1;
    float zero=0;

    clear_fp_flags();
    std::cerr << "1/0 = " << one/zero;print_and_clear();
    std::cerr << "-1/0 = " << -one/zero;print_and_clear();
    std::cerr << "0/0 = " << zero/zero;print_and_clear();
    std::cerr << "log(0) = " << log(zero);print_and_clear();
    std::cerr << "sqrt(-1) = " << sqrt(-one);print_and_clear();
    std::cerr << std::endl;

    std::cerr << "max = " << max;print_and_clear();
    std::cerr << "max+1 = " << max+1;print_and_clear();
    std::cerr << "max-1 = " << max-1;print_and_clear();
    std::cerr << "2*max = " << 2*max;print_and_clear();
    std::cerr << ".5*max = " << 0.5*max;print_and_clear();
    std::cerr << "0*max = " << 0*max;print_and_clear();
    std::cerr << "1/max = " << 1/max;print_and_clear();
    std::cerr << std::endl;

    std::cerr << "inf = " << inf;print_and_clear();
    std::cerr << "inf+1 = " << inf+1;print_and_clear();
    std::cerr << "inf-1 = " << inf-1;print_and_clear();
    std::cerr << "2*inf = " << 2*inf;print_and_clear();
    std::cerr << ".5*inf = " << 0.5*inf;print_and_clear();
    std::cerr << "0*inf = " << 0*inf;print_and_clear();
    std::cerr << "1/inf = " << 1/inf;print_and_clear();
    std::cerr << std::endl;

    std::cerr << "nan = " << nan << std::endl;
    std::cerr << "nan+1 = " << nan+1 << std::endl;
    std::cerr << "nan-1 = " << nan-1 << std::endl;
    std::cerr << "2*nan = " << 2*nan << std::endl;
    std::cerr << ".5*nan = " << 0.5*nan << std::endl;
    std::cerr << "0*nan = " << 0*nan << std::endl;
    std::cerr << "1/nan = " << 1/nan << std::endl;
    std::cerr << std::endl;
}
