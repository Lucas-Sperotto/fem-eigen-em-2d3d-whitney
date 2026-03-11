/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec/edge_mode_post.hpp                                    */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Executavel da formulacao vetorial transversal com elementos de  */
/* aresta.                                                                    */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.1.        */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/mesh2d.hpp"
#include "edge/edge_assembly.hpp"
#include "edge/edge_basis.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace helmvec_post
{
/******************************************************************************/
/* FUNCAO: print_positive_kc                                                  */
/* DESCRICAO: Filtra autovalores fisicos e imprime os primeiros valores de kc */
/* para inspecao rapida do espectro.                                          */
/* ENTRADA: lambdas: const std::vector<double> &; how_many: int; tol: double. */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_positive_kc(
    const std::vector<double> &lambdas,
    int how_many,
    double tol = 1e-9)
{
    int printed = 0;
    for (double lam : lambdas)
    {
        if (lam < tol)
            continue;
        const double kc = std::sqrt(lam);
        std::cout << std::setw(2) << (printed + 1)
                  << "  kc=" << std::fixed << std::setprecision(6) << kc
                  << "\n";
        if (++printed >= how_many)
            break;
    }
}

/******************************************************************************/
/* FUNCAO: pick_first_positive_mode                                           */
/* DESCRICAO: Localiza o primeiro modo fisico acima da tolerancia para evitar */
/* modos nulos/espurios.                                                      */
/* ENTRADA: lambdas: const std::vector<double> &; tol: double.                */
/* SAIDA: int.                                                                */
/******************************************************************************/
inline int pick_first_positive_mode(
    const std::vector<double> &lambdas,
    double tol = 1e-9)
{
    for (int i = 0; i < (int)lambdas.size(); ++i)
    {
        if (lambdas[(size_t)i] > tol)
            return i;
    }
    return -1;
}

// Reconstrucao de um modo transversal de aresta nos centroides dos triangulos:
//   Ft(xc) = sum_m e_m * W_m(xc), m=0..2
// onde e_m ja incorpora o sinal de orientacao local/global da aresta.
/******************************************************************************/
/* FUNCAO: reconstruct_cell_field_from_edge_mode                              */
/* DESCRICAO: Reconstrói o campo transversal por elemento a partir dos        */
/* coeficientes de aresta.                                                    */
/* ENTRADA: mesh: const Mesh2D &; sys: const EdgeSystem &; zcol: const        */
/* std::vector<double> &; mode_idx: int; cell_vx: std::vector<double> &;      */
/* cell_vy: std::vector<double> &; normalize: bool.                           */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void reconstruct_cell_field_from_edge_mode(
    const Mesh2D &mesh,
    const EdgeSystem &sys,
    const std::vector<double> &zcol,
    int mode_idx,
    std::vector<double> &cell_vx,
    std::vector<double> &cell_vy,
    bool normalize = true)
{
    std::array<double, 3> lam = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    cell_vx.assign(mesh.tris.size(), 0.0);
    cell_vy.assign(mesh.tris.size(), 0.0);

    auto idx_col = [&](int i, int j)
    { return (size_t)j * (size_t)sys.ed.ndof + (size_t)i; };

    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const Tri &tri = mesh.tris[(size_t)tid];
        const TriGeomEdge tg = tri_geom_edge(mesh, tri);
        const TriEdges &te = sys.ed.tri_edges[(size_t)tid];

        double e_loc[3] = {0.0, 0.0, 0.0};
        for (int m = 0; m < 3; ++m)
        {
            const int eid = te.e[m];
            const int sgn = te.sgn[m];
            const int dof = sys.ed.edge_to_dof[(size_t)eid];
            const double val = (dof >= 0) ? zcol[idx_col(dof, mode_idx)] : 0.0;
            e_loc[m] = (double)sgn * val;
        }

        double fx = 0.0;
        double fy = 0.0;
        for (int m = 0; m < 3; ++m)
        {
            const Vec2 Wm = whitney_W_local(m, tg, lam);
            fx += e_loc[m] * Wm.x;
            fy += e_loc[m] * Wm.y;
        }
        cell_vx[(size_t)tid] = fx;
        cell_vy[(size_t)tid] = fy;
    }

    if (!normalize)
        return;

    double vmax = 0.0;
    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const double v = std::sqrt(cell_vx[(size_t)tid] * cell_vx[(size_t)tid] +
                                   cell_vy[(size_t)tid] * cell_vy[(size_t)tid]);
        vmax = std::max(vmax, v);
    }
    if (vmax > 0.0)
    {
        for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
        {
            cell_vx[(size_t)tid] /= vmax;
            cell_vy[(size_t)tid] /= vmax;
        }
    }
}
} // namespace helmvec_post
