/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/fem_scalar.hpp                                           */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo escalar 2D (malha, montagem e identificacao modal).      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.1.          */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "assembly_backend.hpp"
#include "mesh2d.hpp"
#include "dense.hpp"
#include "explicit/tri2d_scalar_explicit.hpp"
#include <array>
#include <cmath>

struct TriGeom
{
    double A;                // area do triangulo
    std::array<double, 3> b; // coeficientes b_i das funcoes de forma lineares
    std::array<double, 3> c; // coeficientes c_i das funcoes de forma lineares
};

// Cubatura triangular simetrica de ordem 2.
// Como os integrandos do elemento P1 escalar sao no maximo quadraticos, esta
// regra integra exatamente os termos de rigidez e massa.
constexpr std::array<std::array<double, 3>, 3> kTriGaussP2 = {{
    {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
    {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
    {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
}};

/******************************************************************************/
/* FUNCAO: tri_geom                                                           */
/* DESCRICAO: Calcula area e coeficientes geometricos do triangulo usados na  */
/* montagem escalar.                                                          */
/* ENTRADA: m: const Mesh2D &; t: const Tri &.                                */
/* SAIDA: TriGeom.                                                            */
/******************************************************************************/
inline TriGeom tri_geom(const Mesh2D &m, const Tri &t)
{
    const auto &n1 = m.nodes[t.v[0]];
    const auto &n2 = m.nodes[t.v[1]];
    const auto &n3 = m.nodes[t.v[2]];

    double x1 = n1.x, y1 = n1.y;
    double x2 = n2.x, y2 = n2.y;
    double x3 = n3.x, y3 = n3.y;

    // A = 0.5 * det
    double det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    double A = 0.5 * std::abs(det);

    // bi = yj - yk ; ci = xk - xj (cíclico)
    std::array<double, 3> b = {(y2 - y3), (y3 - y1), (y1 - y2)};
    std::array<double, 3> c = {(x3 - x2), (x1 - x3), (x2 - x1)};

    return {A, b, c};
}

/******************************************************************************/
/* FUNCAO: tri_shape_gradients_scalar                                         */
/* DESCRICAO: Converte os coeficientes geometricos b_i e c_i nos gradientes   */
/* constantes das funcoes de forma nodais do triangulo linear.                */
/* ENTRADA: g: const TriGeom &; dndx: std::array<double, 3> &; dndy:          */
/* std::array<double, 3> &.                                                   */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void tri_shape_gradients_scalar(
    const TriGeom &g,
    std::array<double, 3> &dndx,
    std::array<double, 3> &dndy)
{
    for (int i = 0; i < 3; i++)
    {
        dndx[i] = g.b[i] / (2.0 * g.A);
        dndy[i] = g.c[i] / (2.0 * g.A);
    }
}

/******************************************************************************/
/* FUNCAO: element_mats_scalar_closed_form                                    */
/* DESCRICAO: Calcula as matrizes elementares pela forma fechada da Secao 2.1 */
/* usando diretamente as Eq. (30) e (33) do artigo.                           */
/* ENTRADA: g: const TriGeom &; Sel: double[3][3]; Tel: double[3][3].         */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void element_mats_scalar_closed_form(
    const TriGeom &g,
    double Sel[3][3],
    double Tel[3][3])
{
    explicit_tri2d::tri2d_scalar_closed_form_eq_30_33(g.A, g.b, g.c, Sel, Tel);
}

/******************************************************************************/
/* FUNCAO: element_mats_scalar_gauss                                          */
/* DESCRICAO: Calcula as matrizes elementares via cubatura triangular de 3    */
/* pontos. Para o elemento P1 escalar, esta regra e exata porque integra      */
/* exatamente termos constantes e quadraticos.                                */
/* ENTRADA: g: const TriGeom &; Sel: double[3][3]; Tel: double[3][3].         */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void element_mats_scalar_gauss(
    const TriGeom &g,
    double Sel[3][3],
    double Tel[3][3])
{
    std::array<double, 3> dndx{};
    std::array<double, 3> dndy{};
    tri_shape_gradients_scalar(g, dndx, dndy);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Sel[i][j] = 0.0;
            Tel[i][j] = 0.0;
        }
    }

    for (const auto &lam : kTriGaussP2)
    {
        constexpr double weight = 1.0 / 3.0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                const double grad_dot = dndx[i] * dndx[j] + dndy[i] * dndy[j];
                const double shape_prod = lam[i] * lam[j];
                Sel[i][j] += g.A * weight * grad_dot;
                Tel[i][j] += g.A * weight * shape_prod;
            }
        }
    }
}

/******************************************************************************/
/* FUNCAO: element_mats_scalar                                                */
/* DESCRICAO: Despacha o calculo das matrizes elementares escalares para o    */
/* backend solicitado, permitindo comparar closed-form e quadratura Gauss.    */
/* ENTRADA: g: const TriGeom &; backend: ElementAssemblyBackend; Sel:         */
/* double[3][3]; Tel: double[3][3].                                           */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void element_mats_scalar(
    const TriGeom &g,
    ElementAssemblyBackend backend,
    double Sel[3][3],
    double Tel[3][3])
{
    if (backend == ElementAssemblyBackend::ClosedForm)
    {
        element_mats_scalar_closed_form(g, Sel, Tel);
        return;
    }

    element_mats_scalar_gauss(g, Sel, Tel);
}

/******************************************************************************/
/* FUNCAO: element_mats_scalar                                                */
/* DESCRICAO: Mantem compatibilidade com o codigo legado, preservando a forma */
/* fechada como comportamento padrao quando nenhum backend e informado.       */
/* ENTRADA: g: const TriGeom &; Sel: double[3][3]; Tel: double[3][3].         */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void element_mats_scalar(const TriGeom &g, double Sel[3][3], double Tel[3][3])
{
    element_mats_scalar(g, ElementAssemblyBackend::ClosedForm, Sel, Tel);
}
