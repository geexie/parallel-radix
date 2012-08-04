#include <timing.h>

#if defined WIN32 || defined _WIN32 || defined WINCE
# include <Windows.h>

LARGE_INTEGER timerFreq_;
LARGE_INTEGER counterAtStart_;

#include <tchar.h>
#if defined _MSC_VER
  #if _MSC_VER >= 1400
    #include <intrin.h>
  #elif defined _M_IX86
    static void __cpuid(int* cpuid_data, int)
    {
        __asm
        {
            push ebx
            push edi
            mov edi, cpuid_data
            mov eax, 1
            cpuid
            mov [edi], eax
            mov [edi + 4], ebx
            mov [edi + 8], ecx
            mov [edi + 12], edx
            pop edi
            pop ebx
        }
    }
  #endif
#endif
#else
# include <pthread.h>
# include <sys/time.h>
# include <time.h>

# if defined __MACH__ && defined __APPLE__
#  include <mach/mach.h>
#  include <mach/mach_time.h>
# endif

#endif

#include <stdarg.h>

#if defined __linux__ || defined __APPLE__
# include <unistd.h>
# include <stdio.h>
# include <sys/types.h>
# if defined ANDROID
#  include <sys/sysconf.h>
# else
#  include <sys/sysctl.h>
# endif
#endif


long long radix::getTickCount(void)
{
#if defined WIN32 || defined _WIN32 || defined WINCE
    LARGE_INTEGER counter;
    QueryPerformanceCounter( &counter );
    return (long long)counter.QuadPart;
#elif defined __linux || defined __linux__
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return (long long)tp.tv_sec*1000000000 + tp.tv_nsec;
#elif defined __MACH__ && defined __APPLE__
    return (int64)mach_absolute_time();
#else
    struct timeval tv;
    struct timezone tz;
    gettimeofday( &tv, &tz );
    return (int64)tv.tv_sec*1000000 + tv.tv_usec;
#endif
}

double radix::getTickFrequency(void)
{
#if defined WIN32 || defined _WIN32 || defined WINCE
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return (double)freq.QuadPart;
#elif defined __linux || defined __linux__
    return 1e9;
#elif defined __MACH__ && defined __APPLE__
    static double freq = 0;
    if( freq == 0 )
    {
        mach_timebase_info_data_t sTimebaseInfo;
        mach_timebase_info(&sTimebaseInfo);
        freq = sTimebaseInfo.denom*1e9/sTimebaseInfo.numer;
    }
    return freq;
#else
    return 1e6;
#endif
}
