/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/scalar_mode_post.hpp                                   */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Pos-processamento didatico da formulacao escalar 2D para kc,    */
/* incluindo rastreabilidade entre potencial escalar, gradiente e campos      */
/* eletricos transversais das Eq. (35)-(41).                                  */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.1, Tabelas  */
/* 1-3.                                                                       */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/helm10_scalar_system.hpp"
#include "core/mesh2d.hpp"
#include "meshfree/efgmi_2d.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace helm10_post
{
// Imprime kc = sqrt(lambda) para os primeiros autovalores positivos/fisicos.
/******************************************************************************/
/* FUNCAO: print_positive_kc                                                  */
/* DESCRICAO: Filtra autovalores fisicos e imprime os primeiros valores de kc */
/* para inspecao rapida do espectro.                                          */
/* ENTRADA: lambdas: const std::vector<double> &; how_many: int;              */
/* drop_zero_mode: bool; tol: double.                                         */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_positive_kc(
    const std::vector<double> &lambdas,
    int how_many,
    bool drop_zero_mode,
    double tol = 1e-9)
{
    int printed = 0;
    for (double lam : lambdas)
    {
        if (lam < 0.0)
            continue;
        if (drop_zero_mode && lam < tol)
            continue;
        const double kc = std::sqrt(lam);
        std::cout << std::setw(2) << (printed + 1)
                  << "  kc=" << std::setprecision(6) << std::fixed << kc
                  << "\n";
        if (++printed >= how_many)
            break;
    }
}

// Retorna o primeiro indice com lambda > tol.
// Util para TE quando lambda=0 corresponde ao modo constante.
/******************************************************************************/
/* FUNCAO: pick_first_positive_mode                                           */
/* DESCRICAO: Localiza o primeiro modo fisico acima da tolerancia para evitar */
/* modos nulos/espurios.                                                      */
/* ENTRADA: lambdas: const std::vector<double> &; drop_zero_mode: bool; tol:  */
/* double.                                                                    */
/* SAIDA: int.                                                                */
/******************************************************************************/
inline int pick_first_positive_mode(
    const std::vector<double> &lambdas,
    bool drop_zero_mode,
    double tol = 1e-9)
{
    for (int i = 0; i < (int)lambdas.size(); ++i)
    {
        if (lambdas[(size_t)i] < 0.0)
            continue;
        if (drop_zero_mode && lambdas[(size_t)i] < tol)
            continue;
        return i;
    }
    return -1;
}

// Expande um vetor modal (armazenado em coluna-major LAPACK) para todos os nos
// da malha, reinserindo nos de Dirichlet eliminados com valor zero.
/******************************************************************************/
/* FUNCAO: extract_mode_nodal_from_Z                                          */
/* DESCRICAO: Expande o autovetor para campo nodal completo, reinserindo DOFs */
/* eliminados quando necessario.                                              */
/* ENTRADA: mesh: const Mesh2D &; sys: const ScalarSystem &; zcol: const      */
/* std::vector<double> &; mode_idx: int; normalize: bool.                     */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> extract_mode_nodal_from_Z(
    const Mesh2D &mesh,
    const ScalarSystem &sys,
    const std::vector<double> &zcol,
    int mode_idx,
    bool normalize = true)
{
    std::vector<double> phi(mesh.nodes.size(), 0.0);
    const int n = sys.ndof;
    auto idx_col = [&](int i, int j)
    { return (size_t)j * (size_t)n + (size_t)i; };

    for (size_t ni = 0; ni < mesh.nodes.size(); ++ni)
    {
        const int dof = sys.dof_map[ni];
        phi[ni] = (dof >= 0) ? zcol[idx_col(dof, mode_idx)] : 0.0;
    }

    if (sys.backend == ElementAssemblyBackend::EfgmiConsistent)
    {
        const auto ctx = efgmi2d::make_context(mesh);
        std::vector<double> sampled(mesh.nodes.size(), 0.0);
        for (size_t ni = 0; ni < mesh.nodes.size(); ++ni)
        {
            const Node2D &node = mesh.nodes[ni];
            sampled[ni] = efgmi2d::evaluate_scalar_field(ctx, phi, node.x, node.y, true).value;
        }
        phi = std::move(sampled);
    }

    if (!normalize)
        return phi;

    double max_abs = 0.0;
    for (double v : phi)
        max_abs = std::max(max_abs, std::abs(v));
    if (max_abs > 0.0)
    {
        for (double &v : phi)
            v /= max_abs;
    }
    return phi;
}
} // namespace helm10_post
