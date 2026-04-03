/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec3/main_helmvec3_rect.cpp                               */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema acoplado vetorial+escalar para obter beta dado k0.      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.4, Eq.    */
/* (126)-(127), Fig. 12-13, Tabelas 9-10.                                     */
/*****************************************************************************/
/* Observacao: Este driver desempenha o papel do programa HELMVEC3 do         */
/* apendice em FORTRAN. A montagem global correspondente a Eq. (136) ocorre em*/
/* src/helmvec2/helmvec2_coupled_system.cpp. Comentarios priorizam didatica,  */
/* rastreabilidade e validacao.                                               */
/*****************************************************************************/

#include "article/tp3485_systems.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d.hpp"
#include "core/mesh2d_rect.hpp"
#include "core/timing_utils.hpp"
#include "explicit/tri2d_coupled_explicit.hpp"
#include "helmvec2/helmvec23_shared.hpp"
#include "helmvec2/coupled_cli_options.hpp"
#include "helmvec2/helmvec2_coupled_system.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
struct Table10Block
{
    double d_over_a = 0.0;
    std::vector<double> analytic_beta_over_k0;
    std::vector<double> helmvec3_beta_over_k0;
};

/******************************************************************************/
/* FUNCAO: print_block_3x3                                                    */
/* DESCRICAO: Imprime uma matriz 3x3 em formato compacto para depuracao       */
/* didatica dos blocos locais da formulacao acoplada.                         */
/* ENTRADA: name: const char *; M: const double[3][3].                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void print_block_3x3(const char *name, const double M[3][3])
{
    std::cout << "  " << name << " =\n";
    for (int i = 0; i < 3; ++i)
    {
        std::cout << "   ";
        for (int j = 0; j < 3; ++j)
            std::cout << " " << M[i][j];
        std::cout << "\n";
    }
}

/******************************************************************************/
/* FUNCAO: print_first_triangle_closed_form_debug                             */
/* DESCRICAO: Imprime os blocos locais closed-form do primeiro triangulo para */
/* a formulacao de 2.2.4, mostrando tanto as Eq. (137)-(142) quanto a forma  */
/* rearranjada usada no sistema global Eq. (136).                             */
/* ENTRADA: mesh: const Mesh2D &; k0: double; eps: const std::vector<double>  */
/* &; mu: const std::vector<double> &.                                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void print_first_triangle_closed_form_debug(
    const Mesh2D &mesh,
    double k0,
    const std::vector<double> &eps,
    const std::vector<double> &mu)
{
    if (mesh.tris.empty())
        return;

    const Tri &t = mesh.tris.front();
    const TriGeomEdge tg = tri_geom_edge(mesh, t);
    const auto blk142 = explicit_tri2d::tri2d_beta_closed_form_eq_137_142(
        tg,
        k0,
        eps.front(),
        mu.front());
    const auto blk136 = explicit_tri2d::tri2d_beta_rearranged_closed_form_eq_136(
        tg,
        k0,
        eps.front(),
        mu.front());

    std::cout << "\n[debug] primeiro triangulo: blocos locais closed-form de 2.2.4\n";
    std::cout << "  k0_amostra=" << k0 << "\n";
    std::cout << "  artigo Eq. (137)-(142):\n";
    print_block_3x3("Sel_tt", blk142.Sel_tt);
    print_block_3x3("Tel_tz", blk142.Tel_tz);
    print_block_3x3("Tel_zt", blk142.Tel_zt);
    print_block_3x3("Sel_zz", blk142.Sel_zz);
    print_block_3x3("Tel_tt", blk142.Tel_tt);
    print_block_3x3("Tel_zz", blk142.Tel_zz);

    std::cout << "  rearranjo usado no codigo para Eq. (136):\n";
    print_block_3x3("P_tt", blk136.P_tt);
    print_block_3x3("P_zz", blk136.P_zz);
    print_block_3x3("Q_tt", blk136.Q_tt);
    print_block_3x3("Q_tz", blk136.Q_tz);
    print_block_3x3("Q_zt", blk136.Q_zt);
    print_block_3x3("Q_zz", blk136.Q_zz);
}

/******************************************************************************/
/* FUNCAO: beta_ratio_candidates_from_k0                                      */
/* DESCRICAO: Resolve o sistema de beta para um k0 fixo e filtra candidatos   */
/* fisicos de beta/k0 para construir ramos de dispersao da Secao 2.2.4        */
/* (Eq. 126-127).                                                             */
/* ENTRADA: mesh: const Mesh2D &; eps: const std::vector<double> &; mu: const */
/* std::vector<double> &; k0: double; backend: ElementAssemblyBackend.        */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<double> beta_ratio_candidates_from_k0(
    const Mesh2D &mesh,
    const std::vector<double> &eps,
    const std::vector<double> &mu,
    double k0,
    ElementAssemblyBackend backend,
    timing::Breakdown *perf = nullptr)
{
    timing::Stopwatch stage;
    auto sys = tp3485::build_eq136_helmvec3_system_E(mesh, k0, eps, mu, backend);
    if (perf != nullptr)
        perf->assembly_ms += stage.elapsed_ms();
    stage.reset();
    auto eig = generalized_eigs_real_vec(sys.P, sys.Q);
    if (perf != nullptr)
        perf->solve_ms += stage.elapsed_ms();
    auto beta = helmvec23::collect_positive_real_roots(eig, 1e-4);

    double eps_max = 1.0;
    for (double e : eps)
        eps_max = std::max(eps_max, e);
    const double ratio_max = std::sqrt(eps_max) + 0.25;

    std::vector<double> ratios;
    for (double b : beta)
    {
        const double r = b / k0;
        if (r <= 1e-6 || r > ratio_max)
            continue;
        ratios.push_back(r);
    }
    return helmvec23::unique_sorted(std::move(ratios));
}

/******************************************************************************/
/* FUNCAO: match_ratio_to_reference                                           */
/* DESCRICAO: Associa resultados FEM a referencias analiticas usando criterio */
/* de proximidade em beta/k0 para comparar com Tabela 9 e Tabela 10.          */
/* ENTRADA: mesh: const Mesh2D &; eps: const std::vector<double> &; mu: const */
/* std::vector<double> &; br_over_lambda: const std::vector<double> &; b:     */
/* double; ref_ratio: const std::vector<double> &; debug_candidates: bool;    */
/* backend: ElementAssemblyBackend.                                           */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<double> match_ratio_to_reference(
    const Mesh2D &mesh,
    const std::vector<double> &eps,
    const std::vector<double> &mu,
    const std::vector<double> &br_over_lambda,
    double b,
    const std::vector<double> &ref_ratio,
    bool debug_candidates,
    ElementAssemblyBackend backend,
    timing::Breakdown *perf = nullptr)
{
    std::vector<double> out;
    out.reserve(br_over_lambda.size());

    for (int i = 0; i < (int)br_over_lambda.size(); ++i)
    {
        if (!std::isfinite(ref_ratio[i]))
        {
            out.push_back(std::numeric_limits<double>::quiet_NaN());
            continue;
        }

        const double s = br_over_lambda[i];
        const double k0 = 2.0 * M_PI * s / b;
        auto cands = beta_ratio_candidates_from_k0(mesh, eps, mu, k0, backend, perf);
        if (cands.empty())
        {
            out.push_back(std::numeric_limits<double>::quiet_NaN());
            continue;
        }

        int best = 0;
        double best_err = std::abs(cands[0] - ref_ratio[i]);
        for (int j = 1; j < (int)cands.size(); ++j)
        {
            const double e = std::abs(cands[j] - ref_ratio[i]);
            if (e < best_err)
            {
                best_err = e;
                best = j;
            }
        }
        out.push_back(cands[best]);

        if (debug_candidates)
        {
            std::cout << "  [debug] s=" << s << " cands:";
            for (double c : cands)
                std::cout << " " << c;
            std::cout << "\n";
        }
    }

    return out;
}

/******************************************************************************/
/* FUNCAO: trace_ratio_branch                                                 */
/* DESCRICAO: Rastreia uma rama de dispersao beta/k0 ao longo da variacao de parametro do problema. */
/* Esta rotina implementa continuidade modal numerica entre pontos de         */
/* amostragem consecutivos em b_r/lambda_0.                                   */
/* ENTRADA: mesh: const Mesh2D &; eps: const std::vector<double> &; mu: const */
/* std::vector<double> &; br_over_lambda: const std::vector<double> &; b:     */
/* double; seed_ratio: double; backend: ElementAssemblyBackend.               */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<double> trace_ratio_branch(
    const Mesh2D &mesh,
    const std::vector<double> &eps,
    const std::vector<double> &mu,
    const std::vector<double> &br_over_lambda,
    double b,
    double seed_ratio,
    ElementAssemblyBackend backend,
    timing::Breakdown *perf = nullptr)
{
    std::vector<double> out;
    out.reserve(br_over_lambda.size());

    double prev = seed_ratio;
    for (double s : br_over_lambda)
    {
        const double k0 = 2.0 * M_PI * s / b;
        auto cands = beta_ratio_candidates_from_k0(mesh, eps, mu, k0, backend, perf);
        if (cands.empty())
        {
            out.push_back(std::numeric_limits<double>::quiet_NaN());
            continue;
        }

        int best = 0;
        double best_err = std::abs(cands[0] - prev);
        for (int i = 1; i < (int)cands.size(); ++i)
        {
            const double e = std::abs(cands[i] - prev);
            if (e < best_err)
            {
                best_err = e;
                best = i;
            }
        }
        out.push_back(cands[best]);
        prev = cands[best];
    }
    return out;
}

} // namespace

/******************************************************************************/
/* FUNCAO: main                                                               */
/* DESCRICAO: Ponto de entrada do executavel: interpreta argumentos, prepara o*/
/* caso e executa o fluxo numerico principal.                                 */
/* ENTRADA: argc: int; argv: char **.                                         */
/* SAIDA: int.                                                                */
/******************************************************************************/
int main(int argc, char **argv)
{
    timing::Breakdown perf;
    timing::Stopwatch total_watch;
    // Parametros das Figuras 12/13 da secao 2.2.4.5.
    // Objetivo: obter beta para k0 conhecido e reconstruir curvas de dispersao.
    const double a = 1.0;
    const double b = 0.45;
    const double d12 = 0.5 * b;
    const double eps_fill = 2.45;

    // Configuracoes opcionais em tempo de execucao:
    // argv[1] = d/a para uma pre-visualizacao de continuacao da Figura 13.
    // argv[2], argv[3] = nx, ny for the rectangular mesh.
    // argv[4] = debug legado (0/1).
    double d13_over_a = 0.20;
    int nx = 10, ny = 5; // 100 triangulos, como nos exemplos do artigo.
    helmvec2::CoupledCliOptions cli;
    try
    {
        cli = helmvec2::parse_coupled_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./helmvec3_rect [d_over_a [nx ny [debug]]]"
                  << " [--backend closed-form|gauss]"
                  << " [--debug-local-blocks] [--debug-candidates]\n";
        return 2;
    }

    if (!cli.positionals.empty())
        d13_over_a = std::atof(cli.positionals[0].c_str());
    if (cli.positionals.size() >= 3)
    {
        nx = std::atoi(cli.positionals[1].c_str());
        ny = std::atoi(cli.positionals[2].c_str());
    }
    if (cli.positionals.size() >= 4)
    {
        const bool legacy_debug = (std::atoi(cli.positionals[3].c_str()) != 0);
        if (legacy_debug)
        {
            cli.debug_local_blocks = true;
            cli.debug_candidates = true;
        }
    }

    Mesh2D mesh = make_rect_mesh(a, b, nx, ny);
    std::vector<double> mu(mesh.tris.size(), 1.0);

    const std::vector<double> br_over_lambda_9 = {0.2, 0.3, 0.4, 0.5, 0.6};
    const std::vector<double> ref_ana_9 = {0.48, 1.00, 1.18, 1.26, 1.30};
    const std::vector<double> ref_hvec3_9 = {0.47, 1.01, 1.17, 1.28, 1.35};

    const std::vector<double> br_over_lambda_10 = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<Table10Block> table10 = {
        {0.0, {nan, 0.03, 0.52, 0.70, 0.79, 0.83, 0.88}, {nan, 0.04, 0.56, 0.71, 0.78, 0.83, 0.87}},
        {0.167, {nan, 0.21, 0.60, 0.72, 0.82, 0.88, 0.91}, {nan, 0.18, 0.59, 0.74, 0.81, 0.87, 0.90}},
        {0.286, {nan, 0.51, 0.78, 0.90, 0.99, 1.03, 1.10}, {nan, 0.44, 0.74, 0.88, 1.05, 1.03, 1.09}},
        {0.375, {nan, 0.68, 0.91, 1.05, 1.13, 1.20, 1.25}, {nan, 0.66, 0.90, 1.03, 1.11, 1.18, 1.23}},
        {0.5, {0.40, 0.90, 1.10, 1.20, 1.25, 1.30, 1.35}, {0.42, 0.89, 1.09, 1.19, 1.24, 1.31, 1.35}},
        {0.6, {0.70, 1.02, 1.18, 1.23, 1.31, 1.38, 1.41}, {0.67, 1.03, 1.19, 1.27, 1.33, 1.37, 1.40}},
        {0.8, {0.90, 1.18, 1.29, 1.38, 1.41, 1.43, 1.44}, {0.91, 1.18, 1.30, 1.37, 1.42, 1.44, 1.47}},
    };

    auto eps12 = helmvec23::eps_step_y(mesh, d12, eps_fill, 1.0);
    if (cli.debug_local_blocks)
    {
        const double k0_sample = 2.0 * M_PI * br_over_lambda_9.front() / b;
        print_first_triangle_closed_form_debug(mesh, k0_sample, eps12, mu);
    }
    auto ratio9 = match_ratio_to_reference(
        mesh,
        eps12,
        mu,
        br_over_lambda_9,
        b,
        ref_ana_9,
        cli.debug_candidates,
        cli.backend,
        &perf);
    // Bloco acima corresponde ao caso da Figura 12 (interface horizontal).

    std::cout << "[2.2.4] beta from given k0 (Figure 12)\n";
    std::cout << "a=" << a << " b=" << b << " d=" << d12 << " eps_fill=" << eps_fill
              << " nx=" << nx << " ny=" << ny << " tris=" << mesh.tris.size()
              << " backend=" << element_assembly_backend_name(cli.backend) << "\n";
    std::cout << "br/lambda0  beta/k0(FEM)  Analytic(ref)  HELMVEC3(ref)\n";
    for (int i = 0; i < (int)br_over_lambda_9.size(); ++i)
    {
        std::cout << br_over_lambda_9[i]
                  << "  " << ratio9[i]
                  << "  " << ref_ana_9[i]
                  << "  " << ref_hvec3_9[i] << "\n";
    }

    // Pre-visualizacao opcional de continuacao para um valor de d/a (depuracao visual).
    auto eps13_single = helmvec23::eps_step_x(mesh, d13_over_a * a, eps_fill, 1.0);
    auto ratio13_preview = trace_ratio_branch(
        mesh,
        eps13_single,
        mu,
        br_over_lambda_9,
        b,
        0.5,
        cli.backend,
        &perf);
    // Bloco acima reproduz a logica da Figura 13 (interface vertical).
    std::cout << "\n[2.2.4] beta from given k0 (Figure 13, single d/a preview)\n";
    std::cout << "d/a=" << d13_over_a << "\n";
    std::cout << "br/lambda0  beta/k0(FEM branch)\n";
    for (int i = 0; i < (int)br_over_lambda_9.size(); ++i)
    {
        std::cout << br_over_lambda_9[i] << "  " << ratio13_preview[i] << "\n";
    }

    // Validacao completa alinhada com a Tabela 10.
    std::cout << "\n[2.2.4] Figure 13 / Table 10 validation\n";
    std::cout << "d/a  br/lambda0  beta/k0(FEM matched)  Analytical(ref)  HELMVEC3(ref)\n";
    for (const auto &blk : table10)
    {
        if (cli.debug_candidates)
            std::cout << "  [debug] Figure13 block d/a=" << blk.d_over_a << "\n";

        // Para cada d/a, calcula beta/k0 e faz casamento por proximidade
        // com a familia analitica correspondente da Tabela 10.
        auto eps13 = helmvec23::eps_step_x(mesh, blk.d_over_a * a, eps_fill, 1.0);
        auto fem = match_ratio_to_reference(
            mesh,
            eps13,
            mu,
            br_over_lambda_10,
            b,
            blk.analytic_beta_over_k0,
            cli.debug_candidates,
            cli.backend,
            &perf);

        for (int i = 0; i < (int)br_over_lambda_10.size(); ++i)
        {
            if (!std::isfinite(blk.analytic_beta_over_k0[i]) || !std::isfinite(blk.helmvec3_beta_over_k0[i]))
                continue;

            std::cout << blk.d_over_a
                      << "  " << br_over_lambda_10[i]
                      << "  " << fem[i]
                      << "  " << blk.analytic_beta_over_k0[i]
                      << "  " << blk.helmvec3_beta_over_k0[i] << "\n";
        }
    }

    perf.total_ms = total_watch.elapsed_ms();
    timing::print_breakdown("helmvec3_rect", perf);
    return 0;
}
