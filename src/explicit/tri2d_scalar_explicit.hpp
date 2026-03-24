/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/explicit/tri2d_scalar_explicit.hpp                            */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Formas fechadas (closed-form) para o elemento escalar 2D        */
/* triangular linear da Secao 2.1.                                            */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Eq. (31) e (32).    */
/*****************************************************************************/
/* Observacao: Este cabecalho materializa as formas explicitas usadas como    */
/* referencia analitica local para o backend closed-form.                     */
/*****************************************************************************/

#pragma once

#include <array>

namespace explicit_tri2d
{

/******************************************************************************/
/* FUNCAO: tri2d_scalar_closed_form_eq_30_33                                  */
/* DESCRICAO: Calcula as matrizes elementares escalares do triangulo linear   */
/* usando as formas fechadas do artigo: a rigidez local segue a Eq. (31)      */
/* (usando a definição da equação 30) e a massa consistente segue a Eq. (33). */
/* Essas duas formas tambem sao reaproveitadas nos blocos escalares das Eq.   */
/* (123) e (125) da Secao 2.2.3 e, por rearranjo, nas Eq. (140) e (142) da    */
/* Secao 2.2.4.                                                               */
/* ENTRADA: area: double; b: const std::array<double, 3> &; c: const          */
/* std::array<double, 3> &; Sel: double[3][3]; Tel: double[3][3].             */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void tri2d_scalar_closed_form_eq_30_33(
    double area,
    const std::array<double, 3> &b,
    const std::array<double, 3> &c,
    double Sel[3][3],
    double Tel[3][3])
{
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            // Eq. (30) in Eq. (31): grad-grad escalar no triangulo linear.
            Sel[i][j] = (b[i] * b[j] + c[i] * c[j]) / (4.0 * area);
            // Eq. (33): massa consistente escalar do triangulo linear.
            Tel[i][j] = (area / 12.0) * ((i == j) ? 2.0 : 1.0);
        }
    }
}

} // namespace explicit_tri2d
