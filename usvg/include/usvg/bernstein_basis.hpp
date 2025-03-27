#pragma once

namespace usvg
{
    /**
     * @brief Class that evaluates Bernstein basis polynomials.
     */
    class BernsteinBasis
    {
    public:
        /**
         * @brief Delete constructor.
         */
        BernsteinBasis() = delete;

        /**
         * @brief Evaluates the Bernstein basis polynomial B_i^n(t) = (n i) * (1-t)^(n-i) * t^i
         * @param t Parameter value in [0,1]
         * @param i Index of the basis.
         * @param n Number of basis functions.
         * @return Value of the basis function.
         */
        [[nodiscard]] static double sample(double t, int i, int n) noexcept;

        /**
         * @brief Evaluates the first-order derivative of the Bernstein polynomial.
         * @param t Parameter value in [0,1]
         * @param i Index of the basis.
         * @param n Number of basis functions.
         * @return First derivative of the basis function.
         */
        [[nodiscard]] static double sample_dt(double t, int i, int n) noexcept;

        /**
         * @brief Evaluates the second-order derivative of the Bernstein polynomial.
         * @param t Parameter value in [0,1]
         * @param i Index of the basis.
         * @param n Number of basis functions.
         * @return Second derivative of the basis function.
         */
        [[nodiscard]] static double sample_dtt(double t, int i, int n) noexcept;

    private:
        /**
         * @brief Computes t^i with 0^0=1
         * @param t Base
         * @param i Exponent
         * @return t^i
         */
        static constexpr double safePow(double t, int i) noexcept;

        /**
         * @brief Computes the binomial coefficient recursively at compile time. Returns the number of ways to select a sequence of i distinct objects, retaining the order of selection, from a set of n objects. Note that 0<=i<=n, otherwise the function returns 0.
         * @param n Number of objects to select from.
         * @param i Number of objects to select.
         * @return Binomial coefficient.
         */
        [[nodiscard]] static constexpr int binomial(int n, int i) noexcept;
    };
}
