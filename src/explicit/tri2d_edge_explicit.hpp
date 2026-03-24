/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/explicit/tri2d_edge_explicit.hpp                              */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Formas fechadas (closed-form) do elemento de aresta triangular  */
/* 2D da Secao 2.2.1.                                                         */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Eq. (57)-(77).      */
/*****************************************************************************/
/* Observacao: Este cabecalho organiza os coeficientes A_m, B_m, C_m, D_m e   */
/* os termos It1..It5 para reuso futuro em 2.2.3 e 2.2.4.                     */
/*****************************************************************************/

#pragma once

#include "edge/edge_basis.hpp"

#include <array>
#include <cmath>

namespace explicit_tri2d
{

struct Tri2DEdgeClosedFormData
{
    double area = 0.0;
    double xtri = 0.0;
    double ytri = 0.0;
    double mean_x2 = 0.0; // (1/A) * int_T x^2 dA
    double mean_y2 = 0.0; // (1/A) * int_T y^2 dA
    std::array<double, 3> a{};
    std::array<double, 3> b{};
    std::array<double, 3> c{};
    std::array<double, 3> x{};
    std::array<double, 3> y{};
    std::array<double, 3> edge_len{};
    std::array<double, 3> Acoef{};
    std::array<double, 3> Bcoef{};
    std::array<double, 3> Ccoef{};
    std::array<double, 3> Dcoef{};
};

/******************************************************************************/
/* FUNCAO: tri2d_edge_closed_form_data                                        */
/* DESCRICAO: Constroi os coeficientes geometricos necessarios para as formas */
/* fechadas do elemento de aresta 2D. Em particular, materializa as Eq. (57)  */
/* a (60) e os momentos geometricos usados nas Eq. (73) a (77). Estes dados   */
/* sao a base dos blocos explicitos das Eq. (66), (67), (120)-(125) e, por    */
/* reaproveitamento, tambem das Eq. (137)-(142).                              */
/* ENTRADA: tg: const TriGeomEdge &.                                          */
/* SAIDA: Tri2DEdgeClosedFormData.                                            */
/******************************************************************************/
inline Tri2DEdgeClosedFormData tri2d_edge_closed_form_data(const TriGeomEdge &tg)
{
    Tri2DEdgeClosedFormData data;
    data.area = tg.g.A;
    data.edge_len = tg.L;

    for (int i = 0; i < 3; ++i)
    {
        data.x[i] = tg.X[i].x;
        data.y[i] = tg.X[i].y;
        data.b[i] = tg.g.b[i];
        data.c[i] = tg.g.c[i];
    }

    // Coeficientes a_i da coordenada simplex:
    // lambda_i = (a_i + b_i x + c_i y) / (2A).
    data.a[0] = data.x[1] * data.y[2] - data.x[2] * data.y[1];
    data.a[1] = data.x[2] * data.y[0] - data.x[0] * data.y[2];
    data.a[2] = data.x[0] * data.y[1] - data.x[1] * data.y[0];

    data.xtri = (data.x[0] + data.x[1] + data.x[2]) / 3.0;
    data.ytri = (data.y[0] + data.y[1] + data.y[2]) / 3.0;

    const double sum_x2 =
        data.x[0] * data.x[0] +
        data.x[1] * data.x[1] +
        data.x[2] * data.x[2];
    const double sum_y2 =
        data.y[0] * data.y[0] +
        data.y[1] * data.y[1] +
        data.y[2] * data.y[2];

    // Eq. (76)-(77): termos normalizados por A, isto e, (1/A) * int x^2 dA
    // e (1/A) * int y^2 dA.
    data.mean_x2 = (sum_x2 + 9.0 * data.xtri * data.xtri) / 12.0;
    data.mean_y2 = (sum_y2 + 9.0 * data.ytri * data.ytri) / 12.0;

    const std::array<std::array<int, 2>, 3> edge_nodes = {{
        {{0, 1}},
        {{1, 2}},
        {{2, 0}},
    }};

    for (int m = 0; m < 3; ++m)
    {
        const int i = edge_nodes[m][0];
        const int j = edge_nodes[m][1];
        // Eq. (57): A_m = a_i b_j - a_j b_i.
        data.Acoef[m] = data.a[i] * data.b[j] - data.a[j] * data.b[i];
        // Eq. (58): B_m = c_i b_j - c_j b_i.
        data.Bcoef[m] = data.c[i] * data.b[j] - data.c[j] * data.b[i];
        // Eq. (59): C_m = a_i c_j - a_j c_i.
        data.Ccoef[m] = data.a[i] * data.c[j] - data.a[j] * data.c[i];
        // Eq. (60): D_m = b_i c_j - b_j c_i = -B_m.
        data.Dcoef[m] = data.b[i] * data.c[j] - data.b[j] * data.c[i];
    }

    return data;
}

/******************************************************************************/
/* FUNCAO: tri2d_edge_it_terms_eq_73_77                                       */
/* DESCRICAO: Calcula os cinco termos auxiliares It1..It5 que aparecem na Eq. */
/* (67), ja reduzidos pelas formulas de integracao Eq. (73) a (77). Esses     */
/* termos tambem entram, por composicao, nas Eq. (120), (124), (137) e (141). */
/* ENTRADA: data: const Tri2DEdgeClosedFormData &; m: int; n: int.            */
/* SAIDA: std::array<double, 5>.                                              */
/******************************************************************************/
inline std::array<double, 5> tri2d_edge_it_terms_eq_73_77(
    const Tri2DEdgeClosedFormData &data,
    int m,
    int n)
{
    return {{
        // Eq. (73): It1.
        data.Acoef[m] * data.Acoef[n] + data.Ccoef[m] * data.Ccoef[n],
        // Eq. (74): It2.
        (data.Ccoef[m] * data.Dcoef[n] + data.Ccoef[n] * data.Dcoef[m]) * data.xtri,
        // Eq. (75): It3.
        (data.Acoef[m] * data.Bcoef[n] + data.Acoef[n] * data.Bcoef[m]) * data.ytri,
        // Eq. (76): It4.
        data.Bcoef[m] * data.Bcoef[n] * data.mean_y2,
        // Eq. (77): It5.
        data.Dcoef[m] * data.Dcoef[n] * data.mean_x2,
    }};
}

/******************************************************************************/
/* FUNCAO: tri2d_edge_curlcurl_entry_eq_66                                    */
/* DESCRICAO: Calcula uma entrada local do bloco curl-curl pela forma fechada */
/* da Eq. (66), sem aplicar ainda o fator material 1/mu_r. Este mesmo nucleo  */
/* e reutilizado em Sel(tt) da Eq. (120) e em Sel(tt) da Eq. (137).           */
/* ENTRADA: data: const Tri2DEdgeClosedFormData &; m: int; n: int.            */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double tri2d_edge_curlcurl_entry_eq_66(
    const Tri2DEdgeClosedFormData &data,
    int m,
    int n)
{
    const double scale =
        (data.edge_len[m] * data.edge_len[n]) /
        (4.0 * data.area * data.area * data.area);
    // Eq. (66): integral curl-curl no triangulo.
    return scale * data.Dcoef[m] * data.Dcoef[n];
}

/******************************************************************************/
/* FUNCAO: tri2d_edge_mass_entry_eq_67                                        */
/* DESCRICAO: Calcula uma entrada local do bloco de massa vetorial pela forma */
/* fechada da Eq. (67), sem aplicar ainda o fator material eps_r. Esta forma  */
/* tambem e reutilizada em Tel(tt) da Eq. (124), em Sel(tt) da Eq. (137) e em */
/* Tel(tt) da Eq. (141), com os fatores globais apropriados.                  */
/* ENTRADA: data: const Tri2DEdgeClosedFormData &; m: int; n: int.            */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double tri2d_edge_mass_entry_eq_67(
    const Tri2DEdgeClosedFormData &data,
    int m,
    int n)
{
    const auto it = tri2d_edge_it_terms_eq_73_77(data, m, n);
    const double scale =
        (data.edge_len[m] * data.edge_len[n]) /
        (16.0 * data.area * data.area * data.area);
    // Eq. (67): integral W_m . W_n no triangulo.
    return scale * (it[0] + it[1] + it[2] + it[3] + it[4]);
}

/******************************************************************************/
/* FUNCAO: tri2d_edge_closed_form_eq_66_67                                    */
/* DESCRICAO: Monta os blocos locais de rigidez e massa do elemento de aresta */
/* triangular 2D usando diretamente as Eq. (66) e (67). Estes blocos servem   */
/* de base para os sistemas acoplados das Eq. (120)-(125) e Eq. (137)-(142).  */
/* ENTRADA: tg: const TriGeomEdge &; inv_mu_r: double; eps_r: double; Sel:    */
/* double[3][3]; Tel: double[3][3].                                           */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void tri2d_edge_closed_form_eq_66_67(
    const TriGeomEdge &tg,
    double inv_mu_r,
    double eps_r,
    double Sel[3][3],
    double Tel[3][3])
{
    const Tri2DEdgeClosedFormData data = tri2d_edge_closed_form_data(tg);
    for (int m = 0; m < 3; ++m)
    {
        for (int n = 0; n < 3; ++n)
        {
            // Eq. (66) com ponderacao material 1/mu_r.
            Sel[m][n] = inv_mu_r * tri2d_edge_curlcurl_entry_eq_66(data, m, n);
            // Eq. (67) com ponderacao material eps_r.
            Tel[m][n] = eps_r * tri2d_edge_mass_entry_eq_67(data, m, n);
        }
    }
}

} // namespace explicit_tri2d
