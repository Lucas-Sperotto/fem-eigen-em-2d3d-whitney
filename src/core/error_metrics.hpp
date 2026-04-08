/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/error_metrics.hpp                                        */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Funcoes auxiliares para metricas de erro numerico.              */
/*****************************************************************************/

#pragma once

#include <cmath>
#include <limits>

namespace error_metrics
{

inline double absolute_relative_error_percent(double reference_value, double computed_value)
{
    if (!std::isfinite(reference_value) || !std::isfinite(computed_value))
        return std::numeric_limits<double>::quiet_NaN();

    const double reference_abs = std::abs(reference_value);
    if (reference_abs <= 1.0e-30)
    {
        if (std::abs(computed_value) <= 1.0e-30)
            return 0.0;
        return std::numeric_limits<double>::infinity();
    }

    return 100.0 * std::abs(reference_value - computed_value) / reference_abs;
}

} // namespace error_metrics
