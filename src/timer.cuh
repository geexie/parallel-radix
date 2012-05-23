#ifndef __tnns_timer_h_
#define __tnns_timer_h_

class Timer {
public:
    Timer()
    {
        cudaEventCreate(&_start);
        cudaEventCreate(&_stop);
    }

    ~Timer()
    {
        cudaEventDestroy( _start );
        cudaEventDestroy( _stop );
    }

    void go()
    {
        cudaEventRecord( _start, 0 );
    }

    //returns time after start event to current point in Ms
    float measure()
    {
        cudaEventRecord( _stop, 0 );
        cudaEventSynchronize( _stop );
        float time;
        cudaEventElapsedTime( &time, _start, _stop );
        return time;
    }

private:
    cudaEvent_t _start;
    cudaEvent_t _stop;
};

#endif