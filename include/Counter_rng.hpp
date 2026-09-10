/**
 * @file Counter_rng.hpp
 * @brief Counter-based random numbers, keyed by position in the run.
 *
 * What a draw returns depends only on the key it is asked for -- the run's
 * seed, the iteration, and whatever identifies the thing being drawn for --
 * never on how many draws came before it. That buys three things a sequential
 * generator cannot:
 *
 * - **Thread safety without a lock.** There is no shared state to race on.
 *   The Brownian step used to call a shared std::normal_distribution on a
 *   shared std::mt19937 from inside `#pragma omp parallel for`, which is a
 *   data race on the generator: the noise was neither reproducible nor
 *   guaranteed to still be Gaussian, and an equilibrium fluctuation spectrum
 *   is only ever as good as the noise that drives it.
 * - **Reproducibility independent of scheduling.** A vertex draws the same
 *   number whatever order the OpenMP team reaches it in, and whatever number
 *   of threads.
 * - **Independent streams for free.** The Metropolis flip sweep and the
 *   Brownian displacement can share a seed and still never collide, because
 *   they tag their keys differently.
 *
 * These were local to Dynamic_model.cpp until the flip sweep needed the same
 * primitives. Nothing about them changed in the move; the periodic workload is
 * byte-identical across it.
 */

#pragma once

#include <cmath>
#include <cstdint>

namespace slimed
{

/// SplitMix64 -- one avalanche round on a 64-bit counter.
inline std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/// A uniform on (0, 1) -- never exactly 0, so a log() of it stays finite.
inline double uniform_open01(std::uint64_t bits)
{
    return (static_cast<double>(bits >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

/**
 * @brief One standard normal, keyed by (run, iteration, vertex, axis).
 *
 * Box-Muller on two uniforms drawn from the same key. The cosine branch only,
 * so one call is one number and the key alone determines it.
 */
inline double standard_normal(std::uint64_t stepKey, std::uint64_t vertex, std::uint64_t axis)
{
    const std::uint64_t key = splitmix64(stepKey + splitmix64(vertex * 4ULL + axis));
    const double u1 = uniform_open01(splitmix64(key));
    const double u2 = uniform_open01(splitmix64(key ^ 0xD1B54A32D192ED03ULL));
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

/// A uniform integer in [0, bound), or 0 when the range is empty.
inline int uniform_below(std::uint64_t key, int bound)
{
    if (bound <= 0)
    {
        return 0;
    }
    const int drawn = static_cast<int>(uniform_open01(splitmix64(key)) * bound);
    // uniform_open01() is strictly below 1, but the multiply is not exact.
    return (drawn >= bound) ? bound - 1 : drawn;
}

/**
 * @brief A Poisson draw with mean @p lambda, by Knuth's method.
 *
 * Multiplies uniforms until the product falls below exp(-lambda), which takes
 * about `lambda + 1` iterations. That is the right trade at the rates this is
 * used at -- an edge-flip sweep draws a handful of attempts per step -- and it
 * needs no tables and no state.
 *
 * exp(-lambda) underflows to zero somewhere past lambda = 745, which would
 * spin forever, so the mean is capped. A cap that ever binds means the flip
 * rate or the time step is far outside the regime this was written for, and
 * Mesh::edge_flip_sweep() says so rather than quietly truncating.
 */
inline int poisson(std::uint64_t key, double lambda)
{
    if (!(lambda > 0.0))
    {
        return 0;
    }
    if (lambda > 700.0)
    {
        lambda = 700.0;
    }
    const double threshold = std::exp(-lambda);
    double product = 1.0;
    int count = 0;
    while (true)
    {
        product *= uniform_open01(splitmix64(key + static_cast<std::uint64_t>(count) * 0x9E37ULL));
        if (product <= threshold)
        {
            return count;
        }
        count++;
    }
}

} // namespace slimed
