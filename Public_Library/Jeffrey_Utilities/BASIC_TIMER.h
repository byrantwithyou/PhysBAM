//#####################################################################
// Copyright 2010, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//
// PhysBAM::BASIC_TIMER has a similar interface to boost::timer.  It
// should provide more precise time measurements, at the expense of
// more platform-dependent behavior.
//
// Note: If using QueryPerformanceCounter/QueryPerformanceFrequency
// (the default on Windows), you'll probably want to #define NOMINMAX.
// Note: If using clock_gettime (the default on Linux), you'll have to
// add -lrt to your link command.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_BASIC_TIMER_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_BASIC_TIMER_HPP

#include <exception>

#include <boost/config.hpp>
#include <boost/exception/exception.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/throw_exception.hpp>
#include <boost/type_traits/is_same.hpp>

#if BOOST_VERSION > 103800
#include <boost/exception/errinfo_api_function.hpp>
#include <boost/exception/errinfo_errno.hpp>
#endif // #if BOOST_VERSION > 103800



#if !defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_REALTIME ) || \
    !defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_MONOTONIC ) || \
    !defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_PROCESS ) || \
    !defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_THREAD ) || \
    !defined( PHYSBAM_BASIC_TIMER_USE_POSIX_GETTIMEOFDAY ) || \
    !defined( PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE ) || \
    !defined( PHYSBAM_BASIC_TIMER_USE_CTIME )

// No "clock" chosen, so determine one which is available.

#if   defined( BOOST_HAS_CLOCK_GETTIME )
#if   defined( _POSIX_MONOTONIC_CLOCK )
#define PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_MONOTONIC
#elif defined( _POSIX_CPUTIME ) // #if defined( _POSIX_* )
#define PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_PROCESS
#else // #if defined( _POSIX_* )
#define PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_REALTIME
#endif  // #if defined( _POSIX_* )
#elif defined( BOOST_HAS_GETTIMEOFDAY ) // #if defined( BOOST_* )
#define PHYSBAM_BASIC_TIMER_USE_POSIX_GETTIMEOFDAY
#elif defined( BOOST_WINDOWS ) // #if defined( BOOST_* )
#define PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE
#else // #if defined( BOOST_* )
#define PHYSBAM_BASIC_TIMER_USE_CTIME
#endif // #if defined( BOOST_* )

#elif defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_REALTIME ) + \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_MONOTONIC ) + \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_PROCESS ) + \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_THREAD ) + \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_GETTIMEOFDAY ) + \
      defined( PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE ) + \
      defined( PHYSBAM_BASIC_TIMER_USE_CTIME ) \
      > 1

#error Must define at most one of PHYSBAM_BASIC_TIMER_USE_XXX.

#endif // #if !defined( PHYSBAM_BASIC_TIMER_USE_* ) || ...



#if   defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_REALTIME ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_MONOTONIC ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_PROCESS ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_THREAD )
#include <errno.h>
#include <time.h>
#elif defined( PHYSBAM_BASIC_TIMER_USE_POSIX_GETTIMEOFDAY ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
#include <errno.h>
#include <sys/time.h>
#elif defined( PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
#include <Windows.h>
#ifdef BOOST_HAS_MS_INT64
BOOST_MPL_ASSERT((boost::is_same< LONGLONG, __int64 >));
#endif // #ifdef BOOST_HAS_MS_INT64
#elif defined( PHYSBAM_BASIC_TIMER_USE_CTIME ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
#include <ctime>
#endif // #if defined( PHYSBAM_BASIC_TIMER_USE_* )



namespace PhysBAM
{

//#####################################################################
// struct BASIC_TIMER
//#####################################################################
struct BASIC_TIMER
{
    BASIC_TIMER();

    void Restart();
    double Elapsed() const;

private:
#if   defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_REALTIME ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_MONOTONIC ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_PROCESS ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_THREAD )
    typedef timespec TIME_TYPE;
#elif defined( PHYSBAM_BASIC_TIMER_USE_POSIX_GETTIMEOFDAY ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
    typedef timeval TIME_TYPE;
#elif defined( PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
    typedef LARGE_INTEGER TIME_TYPE;
#elif defined( PHYSBAM_BASIC_TIMER_USE_CTIME ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
    typedef std::clock_t TIME_TYPE;
#endif // #if defined( PHYSBAM_BASIC_TIMER_USE_* )

    static void Get_Current_Time(TIME_TYPE& t);
#ifdef PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE
    static double Freq();
    static double Freq_Impl();
#endif // #ifdef PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE
    TIME_TYPE m_start_time;
};

//#####################################################################
// struct BASIC_TIMER_ERROR
//#####################################################################
struct BASIC_TIMER_ERROR
    : virtual std::exception,
      virtual boost::exception
{
    const char* what() const throw()
    { return "PhysBAM::BASIC_TIMER_ERROR"; }
};

//#####################################################################
// BASIC_TIMER member function implementations
//#####################################################################
inline
BASIC_TIMER::
BASIC_TIMER()
{ Restart(); }

inline void
BASIC_TIMER::
Restart()
{ Get_Current_Time(m_start_time); }

inline double
BASIC_TIMER::
Elapsed() const
{
    TIME_TYPE current_time;
    Get_Current_Time(current_time);
#if   defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_REALTIME ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_MONOTONIC ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_PROCESS ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_THREAD )
    return (current_time.tv_sec - m_start_time.tv_sec) + 1.0e-9 * (current_time.tv_nsec - m_start_time.tv_nsec);
#elif defined( PHYSBAM_BASIC_TIMER_USE_POSIX_GETTIMEOFDAY ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
    return (current_time.tv_sec - m_start_time.tv_sec) + 1.0e-6 * (current_time.tv_usec - m_start_time.tv_usec);
#elif defined( PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
#ifdef BOOST_HAS_MS_INT64
    return (current_time.QuadPart - m_start_time.QuadPart) / Freq();
#else // #ifdef BOOST_HAS_MS_INT64
    static const double double_0x100000000 = static_cast< double >(1UL << 16) * static_cast< double >(1UL << 16);
    return (
        current_time.LowPart >= m_start_time.LowPart ?
        double_0x100000000 * (current_time.HighPart - m_start_time.HighPart) + (current_time.LowPart - m_start_time.LowPart):
        double_0x100000000 * (current_time.HighPart - m_start_time.HighPart) - (m_start_time.LowPart - current_time.LowPart)
    ) / Freq();
#endif // #ifdef BOOST_HAS_MS_INT64
#elif defined( PHYSBAM_BASIC_TIMER_USE_CTIME ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
    return static_cast< double >(current_time - m_start_time) / CLOCKS_PER_SEC;
#endif // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
}

inline void
BASIC_TIMER::
Get_Current_Time(TIME_TYPE& t)
{
#if   defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_REALTIME ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_MONOTONIC ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_PROCESS ) || \
      defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_THREAD )
#if   defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_REALTIME )
    const int e = clock_gettime(CLOCK_REALTIME, &t);
#elif defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_MONOTONIC ) // #if defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_* )
    const int e = clock_gettime(CLOCK_MONOTONIC, &t);
#elif defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_PROCESS ) // #if defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_* )
    const int e = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
#elif defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_GETTIME_THREAD ) // #if defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_* )
    const int e = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
#endif // #if defined( PHYSBAM_BASIC_TIMER_USE_POSIX_CLOCK_* )
    if(e == -1)
#if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR()
            << boost::errinfo_api_function("int clock_gettime(clockid_t,timespec*)")
            << boost::errinfo_errno(errno)
        );
#else // #if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR());
#endif // #if BOOST_VERSION > 103800
#elif defined( PHYSBAM_BASIC_TIMER_USE_POSIX_GETTIMEOFDAY ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
    const int e = gettimeofday(&t, NULL);
    if(e == -1)
#if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR()
            << boost::errinfo_api_function("int gettimeofday(timeval*,timezone*)")
            << boost::errinfo_errno(errno)
        );
#else // #if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR());
#endif // #if BOOST_VERSION > 103800
#elif defined( PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
    const BOOL b = QueryPerformanceCounter(&t);
    if(!b)
#if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR()
            << boost::errinfo_api_function("BOOL QueryPerformanceCounter(LARGE_INTEGER*)")
        );
#else // #if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR());
#endif // #if BOOST_VERSION > 103800
#elif defined( PHYSBAM_BASIC_TIMER_USE_CTIME ) // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
    t = std::clock();
    if(t == -1)
#if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR()
            << boost::errinfo_api_function("std::clock_t std::clock()")
        );
#else // #if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR());
#endif // #if BOOST_VERSION > 103800
#endif // #if defined( PHYSBAM_BASIC_TIMER_USE_* )
}

#ifdef PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE

inline double
BASIC_TIMER::
Freq()
{
    static const double result = Freq_Impl();
    return result;
}

inline double
BASIC_TIMER::
Freq_Impl()
{
    LARGE_INTEGER queried_freq;
    const BOOL b = QueryPerformanceFrequency(&queried_freq);
    if(!b)
#if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR()
            << boost::errinfo_api_function("BOOL QueryPerformanceFrequency(LARGE_INTEGER*)")
        );
#else // #if BOOST_VERSION > 103800
        BOOST_THROW_EXCEPTION(BASIC_TIMER_ERROR());
#endif // #if BOOST_VERSION > 103800
#ifdef BOOST_HAS_MS_INT64
    return static_cast< double >(queried_freq.QuadPart);
#else // #ifdef BOOST_HAS_MS_INT64
    const double double_0x100000000 = static_cast< double >(1UL << 16) * static_cast< double >(1UL << 16);
    return double_0x100000000 * queried_freq.HighPart + queried_freq.LowPart;
#endif // #ifdef BOOST_HAS_MS_INT64
}

#endif // #ifdef PHYSBAM_BASIC_TIMER_USE_WINDOWS_QUERY_PERFORMANCE

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_BASIC_TIMER_HPP
