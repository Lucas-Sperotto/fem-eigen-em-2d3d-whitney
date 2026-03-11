/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec2/helmvec23_shared.hpp                                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema acoplado vetorial+escalar para obter k0 dado beta.      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.3, Eq.    */
/* (108)-(109), Fig. 11, Tabela 8.                                            */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/lapack_eig.hpp"
#include "core/mesh2d.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace helmvec23
{
/******************************************************************************/
/* FUNCAO: unique_sorted                                                      */
/* DESCRICAO: Remove duplicatas e ordena valores para estabilizar listas de candidatos espectrais. */
/* ENTRADA: v: std::vector<double>; tol: double.                              */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> unique_sorted(std::vector<double> v, double tol = 1e-8)
{
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end(), [tol](double a, double b)
                        { return std::abs(a - b) < tol; }),
            v.end());
    return v;
}

// Mantem apenas raizes reais positivas do PAV real generalizado.
/******************************************************************************/
/* FUNCAO: collect_positive_real_roots                                        */
/* DESCRICAO: Filtra autovalores reais positivos e converte para raizes       */
/* fisicas do problema (k0 ou beta), conforme Eq. 108-109 e Eq. 126-127.     */
/* ENTRADA: res: const GenEigGeneralResult &; imag_tol: double; lambda_min:   */
/* double.                                                                    */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> collect_positive_real_roots(
    const GenEigGeneralResult &res,
    double imag_tol = 1e-7,
    double lambda_min = 1e-10)
{
    std::vector<double> out;
    out.reserve((size_t)res.n);
    for (int i = 0; i < res.n; ++i)
    {
        if (!std::isfinite(res.lambda_re[i]))
            continue;
        if (std::abs(res.lambda_im[i]) > imag_tol)
            continue;
        if (res.lambda_re[i] <= lambda_min)
            continue;
        out.push_back(std::sqrt(res.lambda_re[i]));
    }
    return unique_sorted(std::move(out));
}

// Perfil de material com interface horizontal:
//   eps = eps_below se y_centroid < y_split, senao eps_above.
/******************************************************************************/
/* FUNCAO: eps_step_y                                                         */
/* DESCRICAO: Define perfil de permissividade por degrau na direcao y para casos parcialmente preenchidos. */
/* Este perfil reproduz a geometria estratificada da Figura 11 e Figura 12.   */
/* ENTRADA: mesh: const Mesh2D &; y_split: double; eps_below: double;         */
/* eps_above: double.                                                         */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> eps_step_y(
    const Mesh2D &mesh,
    double y_split,
    double eps_below,
    double eps_above)
{
    std::vector<double> eps(mesh.tris.size(), eps_above);
    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const Tri &t = mesh.tris[(size_t)tid];
        const double yc =
            (mesh.nodes[(size_t)t.v[0]].y +
             mesh.nodes[(size_t)t.v[1]].y +
             mesh.nodes[(size_t)t.v[2]].y) /
            3.0;
        eps[(size_t)tid] = (yc < y_split) ? eps_below : eps_above;
    }
    return eps;
}

// Perfil de material com interface vertical:
//   eps = eps_left se x_centroid < x_split, senao eps_right.
/******************************************************************************/
/* FUNCAO: eps_step_x                                                         */
/* DESCRICAO: Define perfil de permissividade por degrau na direcao x para casos parcialmente preenchidos. */
/* Este perfil reproduz a geometria com interface vertical da Figura 13.      */
/* ENTRADA: mesh: const Mesh2D &; x_split: double; eps_left: double;          */
/* eps_right: double.                                                         */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> eps_step_x(
    const Mesh2D &mesh,
    double x_split,
    double eps_left,
    double eps_right)
{
    std::vector<double> eps(mesh.tris.size(), eps_right);
    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const Tri &t = mesh.tris[(size_t)tid];
        const double xc =
            (mesh.nodes[(size_t)t.v[0]].x +
             mesh.nodes[(size_t)t.v[1]].x +
             mesh.nodes[(size_t)t.v[2]].x) /
            3.0;
        eps[(size_t)tid] = (xc < x_split) ? eps_left : eps_right;
    }
    return eps;
}

/******************************************************************************/
/* FUNCAO: pick_closest_unused                                                */
/* DESCRICAO: Seleciona o candidato nao utilizado mais proximo de um valor de */
/* referencia.                                                                */
/* ENTRADA: target: double; cands: const std::vector<double> &; used:         */
/* std::vector<char> &.                                                       */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double pick_closest_unused(
    double target,
    const std::vector<double> &cands,
    std::vector<char> &used)
{
    int best = -1;
    double best_err = std::numeric_limits<double>::infinity();
    for (int i = 0; i < (int)cands.size(); ++i)
    {
        if (used[(size_t)i])
            continue;
        const double e = std::abs(cands[(size_t)i] - target);
        if (e < best_err)
        {
            best_err = e;
            best = i;
        }
    }
    if (best < 0)
        return std::numeric_limits<double>::quiet_NaN();
    used[(size_t)best] = 1;
    return cands[(size_t)best];
}

} // namespace helmvec23
