#pragma once
#include<chrono>
#include<string>
#include<vector>
#include<utility>
#include"../Types/Types.h"
#include"../Parameter_Space/Parameter_Space.h"

namespace ps = Parameter_Space;

namespace Time_Measurement
{

class Clock
{
    private:
    int my_rank;
    std::chrono::steady_clock::time_point program_start;
    std::chrono::steady_clock::time_point last_measurement;
    std::vector<std::pair<std::string, double>> m_measurements; // (task_name, duration_s)
    double m_total_s{0.};

    public:
    Clock() = default;
    Clock( const int rank );

    void measure( const std::string& task );
    void finalize();

    const std::vector<std::pair<std::string, double>>& get_measurements() const { return m_measurements; }
    double get_total_s() const { return m_total_s; }
};

class Simple_Estimator
{
    private:
    int my_rank;
    size_t num_iterations;
    std::chrono::steady_clock::time_point start_time;
    std::string task_name;

    public:
    Simple_Estimator() = default;
    Simple_Estimator( const int rank, const size_t num_iterations, std::string task_name );

    void enter_loop();
    void leave_loop();
    void estimate( const int iteration );

};



}