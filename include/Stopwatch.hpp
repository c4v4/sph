#ifndef SPH_INCLUDE_STOPWATCH_HPP_
#define SPH_INCLUDE_STOPWATCH_HPP_

#include <chrono>

class Stopwatch {
public:
    Stopwatch() : begin(std::chrono::steady_clock::now()), last(begin) { }

    double seconds_from_begin() { return std::chrono::duration<double>(now() - begin).count(); }
    double millisec_from_begin() { return std::chrono::duration<double, std::milli>(now() - begin).count(); }
    double microsec_from_begin() { return std::chrono::duration<double, std::micro>(now() - begin).count(); }
    double nanosec_from_begin() { return std::chrono::duration<double, std::nano>(now() - begin).count(); }

    double seconds_lap() {
        auto old_last = last;
        last = now();
        return std::chrono::duration<double>(last - old_last).count();
    }

    double millisec_lap() {
        auto old_last = last;
        last = now();
        return std::chrono::duration<double, std::milli>(last - old_last).count();
    }

    double microsec_lap() {
        auto old_last = last;
        last = now();
        return std::chrono::duration<double, std::micro>(last - old_last).count();
    }

    double nanosec_lap() {
        auto old_last = last;
        last = now();
        return std::chrono::duration<double, std::nano>(last - old_last).count();
    }

protected:
    std::chrono::steady_clock::time_point now() { return std::chrono::steady_clock::now(); }

private:
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point last;
};


class Timer : public Stopwatch {

public:
    Timer() : end(now()) { }
    Timer(double second_tlim) : end(now() + std::chrono::microseconds(static_cast<uint64_t>(second_tlim * 1E6))) { }

    double seconds_until_end() { return std::chrono::duration<double>(end - now()).count(); }
    double millisec_until_end() { return std::chrono::duration<double, std::milli>(end - now()).count(); }
    double microsec_until_end() { return std::chrono::duration<double, std::micro>(end - now()).count(); }
    double nanosec_until_end() { return std::chrono::duration<double, std::nano>(end - now()).count(); }

    bool exceeded_tlim() { return now() >= end; }

private:
    std::chrono::steady_clock::time_point end;
};

#endif