#ifndef TICTOC_H_
#define TICTOC_H_

#include <ctime>
#include <cstdlib>
#include <chrono>

class TicToc
{
public:
    TicToc()
    {
        Tic();
    }

    void Tic()
    {
        start = std::chrono::system_clock::now();
    }

    double Toc()
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> duration_time = end -start;
        return duration_time.count() * 1000;
    }
private:
    std::chrono::time_point<std::chrono::system_clock> start, end;
};

#endif