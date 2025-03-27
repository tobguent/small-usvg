#pragma once

#include <stdint.h>

namespace usvg
{
    /**
     * @brief Types of boundary conditions.
     */
    enum EBoundaryCondition : uint8_t
    {
        /**
         * @brief Dirichlet boundary condition.
         */
        Dirichlet,
        /**
         * @brief Neumann boundary condition.
         */
        Neumann
    };
}
