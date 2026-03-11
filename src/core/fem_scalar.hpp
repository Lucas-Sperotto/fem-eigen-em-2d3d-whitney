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
#include "mesh2d.hpp"
#include "dense.hpp"
#include <array>
#include <cmath>

struct TriGeom
{
    double A;                // area do triangulo
    std::array<double, 3> b; // coeficientes b_i das funcoes de forma lineares
    std::array<double, 3> c; // coeficientes c_i das funcoes de forma lineares
};

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
/* FUNCAO: element_mats_scalar                                                */
/* DESCRICAO: Calcula as matrizes elementares de rigidez e massa do elemento triangular escalar (Secao 2.1). */
/* ENTRADA: g: const TriGeom &; Sel: double; Tel: double.                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void element_mats_scalar(const TriGeom &g, double Sel[3][3], double Tel[3][3])
{
    // Sel_ij = \int grad(ai)·grad(aj) dA
    // grad(ai) = [bi, ci] / (2A)  => Sel_ij = (bi*bj+ci*cj)/(4A)
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Sel[i][j] = (g.b[i] * g.b[j] + g.c[i] * g.c[j]) / (4.0 * g.A);
        }
    }

    // Tel (massa consistente) = (A/12) * [[2,1,1],[1,2,1],[1,1,2]]
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Tel[i][j] = (g.A / 12.0) * ((i == j) ? 2.0 : 1.0);
        }
    }
}
