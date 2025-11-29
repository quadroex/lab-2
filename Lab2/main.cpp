// Compiler: MSVC
// Standart: C++20

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <execution>
#include <random>
#include <chrono>
#include <format>
#include <cmath>
#include <thread>
#include <iomanip>
#include <functional>

template <class F>
double time_count(F&& f)
{
    auto time1 = std::chrono::high_resolution_clock::now();
    f();
    auto time2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_time = time2 - time1;

    std::cout << std::fixed << std::setprecision(3) << ms_time.count() << " ms";

    return ms_time.count();
}

std::vector<double> generate_data(size_t size)
{
    std::vector<double> data(size);

    std::mt19937 gen;
    std::uniform_real_distribution<> dist(-1488.0, 1488.0);

    std::generate(data.begin(), data.end(), [&]() { return dist(gen); });

    return data;
}

// finishes unprocessed sequence remains
void worker_thread(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end, double& result)
{
    if (start >= end || std::distance(start, end) < 2)
    {
        result = 0.0;
        return;
    }

    result = std::transform_reduce(std::execution::seq,
        start + 1, end,
        start,
        0.0,
        std::plus<>(),
        [](double curr, double prev) { return std::abs(curr - prev); }
    );
}

template <typename Policy>
void run_library_algo(Policy policy, const std::vector<double>& data, const std::string& name)
{
    size_t n = data.size();
    if (n < 2) return;

    std::cout << "For " << name << " time = ";

    time_count([&]()
        {
            std::vector<double> diffs(n);
            std::adjacent_difference(policy, data.begin(), data.end(), diffs.begin());

            double sum = std::transform_reduce
            (
                policy,
                diffs.begin() + 1, diffs.end(),
                0.0,
                std::plus<>(),
                [](double val) { return std::abs(val); }
            );

        });
    std::cout << std::endl;
}

double run_custom_algo(const std::vector<double>& data, int K, size_t n) {
    if (n < 2) return 0.0;

    size_t total_diffs = n - 1;
    size_t threads_diffs = total_diffs / K;
    size_t remained = total_diffs % K;

    std::vector<std::thread> threads;
    std::vector<double> partial_results(K, 0.0);
    size_t start_index = 0;

    double duration = time_count([&]() {
        for (int i = 0; i < K; ++i) {
            size_t count_diffs_for_thread = threads_diffs + (i < remained ? 1 : 0);
            if (count_diffs_for_thread == 0) continue;

            auto it_begin = data.begin() + start_index;
            auto it_end = it_begin + count_diffs_for_thread + 1;

            threads.emplace_back(worker_thread, it_begin, it_end, std::ref(partial_results[i]));
            start_index += count_diffs_for_thread;
        }

        for (auto& t : threads) {
            if (t.joinable()) t.join();
        }

        double total_sum = std::reduce(std::execution::seq, partial_results.begin(), partial_results.end()); // result
        });

    return duration;
}

int main() {
    int hw_threads = std::thread::hardware_concurrency();

    std::cout << "Hardware concurrency: " << hw_threads << " threads.\n";
    std::vector<size_t> data_sizes = { 1000000, 10000000, 100000000 };

    for (size_t size : data_sizes) {
        std::cout << "\nData size: " << size << " elements.\n";
        auto sequence = generate_data(size);
        std::cout << "\nLibrary algorithms:\n";

        run_library_algo(std::execution::seq, sequence, "std::seq");
        run_library_algo(std::execution::par, sequence, "std::par");
        run_library_algo(std::execution::unseq, sequence, "std::unseq");
        run_library_algo(std::execution::par_unseq, sequence, "std::par_unseq");

        std::cout << "\nCustom algorithm:\n";
        std::cout << "Threads numbder (K) accomplished in Time (in ms)\n";

        int max_k = hw_threads * 3;
        double time;
        double min_time;
        int best_k = 1;

        for (int k = 1; k <= max_k; ++k)
        {
            std::cout << "For K = " << k << " time = ";
            time = run_custom_algo(sequence, k, size);
            std::cout << std::endl;

            if (k == 1 || time < min_time)
            {
                min_time = time;
                best_k = k;
            }
        }

        std::cout << "\nBest K = " << best_k << " (" << min_time << " ms). K / CPU threads: "  << (double)best_k / hw_threads << "\n\n";
    }

    std::cout << "Done.";
    return 0;
}