#include <usvg/bernstein_basis.hpp>

#include <cmath>

namespace usvg
{
    double BernsteinBasis::sample(double t, int i, int n) noexcept
    {
        return binomial(n, i) * safePow(1 - t, n - i) * safePow(t, i);
    }

    double BernsteinBasis::sample_dt(double t, int i, int n) noexcept
    {
        return binomial(n, i) * (i * safePow(1 - t, n - i) * safePow(t, i - 1) - (n - i) * safePow(1 - t, n - i - 1) * safePow(t, i));
    }

    double BernsteinBasis::sample_dtt(double t, int i, int n) noexcept
    {
        return binomial(n, i) * ((n - i - 1) * (n - i) * safePow(1 - t, n - i - 2) * safePow(t, i) - 2 * i * (n - i) * safePow(1 - t, n - i - 1) * safePow(t, i - 1) + (i - 1) * i * safePow(1 - t, n - i) * safePow(t, i - 2));
    }

    constexpr double BernsteinBasis::safePow(double t, int i) noexcept
    {
        if (i < 0)
            return 0;
        if (t == 0 && i == 0)
            return 1;
        return std::pow(t, i);
    }

    constexpr int BernsteinBasis::binomial(int n, int i) noexcept
    {
        if (i > n || i < 0 || n < 0)
            return 0;
        if (i == 0 || i == n)
            return 1;
        return binomial(n - 1, i - 1) + binomial(n - 1, i);
    }
}
