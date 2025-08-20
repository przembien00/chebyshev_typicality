#include<chrono>
#include<string>
#include"Time_Measurement.h"
#include"../cpp_libs/Print_Routines.h"

namespace ps = Parameter_Space;
namespace print = Print_Routines;

namespace Time_Measurement
{

Clock::Clock( const int rank ):
my_rank(rank)
{
    program_start = std::chrono::steady_clock::now();
    last_measurement = program_start;
}

void Clock::measure( const std::string& task )
{
    auto now = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(now - last_measurement).count();
    if( duration < 1 )
    {
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_measurement).count();
        print::print_R0(my_rank,"\033[1;36m" + task + " done, took " + std::to_string(duration) + " ms\033[0m\n");
        last_measurement = now;
        return;
    }
    last_measurement = now;
    print::print_R0(my_rank, "\033[1;36m" + task + " done, took " + std::to_string(duration) + " s\033[0m\n");
}

void Clock::finalize()
{
    auto end = std::chrono::steady_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end - program_start).count();
    if( total_duration < 1 )
    {
        total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - program_start).count();
        print::print_R0(my_rank, "Total program duration: " + std::to_string(total_duration) + " ms\n");
        return;
    }
    print::print_R0(my_rank, "Total program duration: " + std::to_string(total_duration) + " s\n");
}

Simple_Estimator::Simple_Estimator( const int rank, const size_t num_iterations, std::string task_name ):
my_rank(rank),
num_iterations(num_iterations),
task_name(task_name)
{}

void Simple_Estimator::enter_loop()
{
    start_time = std::chrono::steady_clock::now();
}

void Simple_Estimator::estimate( const int iteration )
{
    if( iteration == 0 )
    {
        auto now = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count();
        auto estimated_duration = duration * num_iterations / 1000; // in seconds
        print::print_R0(my_rank, "Single iteration duration: " + std::to_string(duration) + " ms\n");
        print::print_R0(my_rank,"Estimated duration of " + task_name + ": " + std::to_string(estimated_duration) + " s\n");
        print::print_R0(my_rank, "---------- \033[1mProgress of " + task_name + "\033[0m ----------\n");
    }
    RealType progress = static_cast<RealType>(iteration + 1) / static_cast<RealType>(num_iterations);
    int bar_width = 46; // Width of the progress bar
    print::print_R0(my_rank, "\r[");
    int filled = static_cast<int>(bar_width * progress);
    for (int i = 0; i < filled; ++i)
    {
        print::print_R0(my_rank, "=");
    }
    for (int i = filled; i < bar_width; ++i)
    {
        print::print_R0(my_rank, " ");
    }
    print::print_R0(my_rank, "] " + std::to_string(static_cast<int>(progress * 100)) + "%");
    std::cout.flush(); // Carriage return to overwrite the line
}

void Simple_Estimator::leave_loop()
{
    print::print_R0(my_rank, "\n-----------------------------------------------------\n");
}

}