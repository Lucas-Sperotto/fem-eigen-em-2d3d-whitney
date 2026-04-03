/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec2/main_helmvec2_rect.cpp                               */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema acoplado vetorial+escalar para obter k0 dado beta.      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.3, Eq.    */
/* (108)-(109), Fig. 11, Tabela 8.                                            */
/*****************************************************************************/
/* Observacao: Este driver desempenha o papel do programa HELMVEC2 do         */
/* apendice em FORTRAN. A montagem global correspondente a Eq. (119) ocorre em*/
/* src/helmvec2/helmvec2_coupled_system.cpp. Comentarios priorizam didatica,  */
/* rastreabilidade e validacao.                                               */
/*****************************************************************************/

#include "article/tp3485_systems.hpp"
#include "core/mesh2d.hpp"
#include "core/mesh2d_rect.hpp"
#include "core/lapack_eig.hpp"
#include "core/timing_utils.hpp"
#include "explicit/tri2d_coupled_explicit.hpp"
#include "helmvec23_shared.hpp"
#include "coupled_cli_options.hpp"
#include "helmvec2_coupled_system.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

struct ModeCand
{
    double k0 = 0.0;
    double ez_ratio = 0.0; // ||Ez||^2 / (||Et||^2 + ||Ez||^2)
};

/******************************************************************************/
/* FUNCAO: print_block_3x3                                                    */
/* DESCRICAO: Imprime uma matriz 3x3 em formato compacto para depuracao       */
/* didatica dos blocos locais da formulacao acoplada.                         */
/* ENTRADA: name: const char *; M: const double[3][3].                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static void print_block_3x3(const char *name, const double M[3][3])
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
/* a formulacao de 2.2.3, mostrando tanto as Eq. (120)-(125) quanto a forma  */
/* rearranjada usada no sistema global Eq. (119).                             */
/* ENTRADA: mesh: const Mesh2D &; beta: double; eps: const std::vector<double>*/
/* &; mu: const std::vector<double> &.                                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static void print_first_triangle_closed_form_debug(
    const Mesh2D &mesh,
    double beta,
    const std::vector<double> &eps,
    const std::vector<double> &mu)
{
    if (mesh.tris.empty())
        return;

    const Tri &t = mesh.tris.front();
    const TriGeomEdge tg = tri_geom_edge(mesh, t);
    const auto blk125 = explicit_tri2d::tri2d_wavenumber_closed_form_eq_120_125(
        tg,
        beta,
        eps.front(),
        mu.front());
    const auto blk119 = explicit_tri2d::tri2d_wavenumber_rearranged_closed_form_eq_119(
        tg,
        beta,
        eps.front(),
        mu.front());

    std::cout << "\n[debug] primeiro triangulo: blocos locais closed-form de 2.2.3\n";
    std::cout << "  beta_amostra=" << beta << "\n";
    std::cout << "  artigo Eq. (120)-(125):\n";
    print_block_3x3("Sel_tt", blk125.Sel_tt);
    print_block_3x3("Sel_tz", blk125.Sel_tz);
    print_block_3x3("Sel_zt", blk125.Sel_zt);
    print_block_3x3("Sel_zz", blk125.Sel_zz);
    print_block_3x3("Tel_tt", blk125.Tel_tt);
    print_block_3x3("Tel_zz", blk125.Tel_zz);

    std::cout << "  rearranjo usado no codigo para Eq. (119):\n";
    print_block_3x3("A_tt", blk119.A_tt);
    print_block_3x3("A_tz", blk119.A_tz);
    print_block_3x3("A_zt", blk119.A_zt);
    print_block_3x3("A_zz", blk119.A_zz);
    print_block_3x3("B_tt", blk119.B_tt);
    print_block_3x3("B_zz", blk119.B_zz);
}

/******************************************************************************/
/* FUNCAO: collect_mode_candidates                                            */
/* DESCRICAO: Extrai candidatos fisicos do espectro de A x = lambda B x,      */
/* usando lambda = k0^2 (Secao 2.2.3), e estima energia relativa no bloco Ez. */
/* ENTRADA: res: const GenEigGeneralResult &; nt: int; nz: int; imag_tol:     */
/* double.                                                                    */
/* SAIDA: std::vector<ModeCand>.                                              */
/******************************************************************************/
static std::vector<ModeCand> collect_mode_candidates(
    const GenEigGeneralResult &res,
    int nt,
    int nz,
    double imag_tol = 1e-7)
{
    std::vector<ModeCand> out;
    const int n = res.n;
    for (int i = 0; i < res.n; ++i)
    {
        if (!std::isfinite(res.lambda_re[i]))
            continue;
        if (std::abs(res.lambda_im[i]) > imag_tol)
            continue;
        if (res.lambda_re[i] <= 1e-10)
            continue;

        const double k0 = std::sqrt(res.lambda_re[i]);
        double et = 0.0, ez = 0.0;
        for (int r = 0; r < nt; ++r)
        {
            const double v = res.VRcol[(size_t)i * n + r];
            et += v * v;
        }
        for (int r = 0; r < nz; ++r)
        {
            const double v = res.VRcol[(size_t)i * n + (nt + r)];
            ez += v * v;
        }
        const double den = et + ez;
        const double ratio = (den > 0.0) ? (ez / den) : 0.0;
        out.push_back({k0, ratio});
    }
    std::sort(out.begin(), out.end(), [](const ModeCand &a, const ModeCand &b)
              { return a.k0 < b.k0; });
    return out;
}

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
    // Secao 2.2.3.5 / Figura 11: quadrado, parte superior eps=1.0 e inferior eps=1.5.
    // Objetivo: reproduzir Tabela 8 com beta fixo e extrair k0L dos modos LSM.
    const double L = 1.0;
    int nx = 6, ny = 6; // 72 triangles
    double beta = 10.0; // beta*L = 10
    helmvec2::CoupledCliOptions cli;
    try
    {
        cli = helmvec2::parse_coupled_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./helmvec2_rect [beta [nx ny [debug]]]"
                  << " [--backend closed-form|gauss]"
                  << " [--debug-local-blocks] [--debug-candidates]\n";
        return 2;
    }

    if (!cli.positionals.empty())
        beta = std::atof(cli.positionals[0].c_str());
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

    auto mesh = make_rect_mesh(L, L, nx, ny);
    auto eps = helmvec23::eps_step_y(mesh, 0.5 * L, 1.5, 1.0);
    std::vector<double> mu(mesh.tris.size(), 1.0);

    if (cli.debug_local_blocks)
        print_first_triangle_closed_form_debug(mesh, beta, eps, mu);

    timing::Stopwatch stage;
    auto sys = tp3485::build_eq119_helmvec2_system_E(mesh, beta, eps, mu, cli.backend);
    perf.assembly_ms += stage.elapsed_ms();
    stage.reset();
    auto eig = generalized_eigs_real_vec(sys.A, sys.B);
    perf.solve_ms += stage.elapsed_ms();
    auto all_modes = collect_mode_candidates(eig, sys.nt, sys.nz);

    double eps_max = 1.0;
    for (double e : eps)
        eps_max = std::max(eps_max, e);
    const double k0_min_phys = beta / std::sqrt(eps_max);
    // Regiao propagante: beta < k0*sqrt(eps_r_max)  => k0 > beta/sqrt(eps_r_max).
    // Esse filtro remove raizes nao propagantes para o meio estratificado da Fig. 11.

    // Mantem apenas raizes fisicamente propagantes para este beta.
    std::vector<double> k0_phys;
    for (const auto &m : all_modes)
    {
        if (m.k0 <= k0_min_phys)
            continue;
        k0_phys.push_back(m.k0);
    }
    k0_phys = helmvec23::unique_sorted(std::move(k0_phys));

    const std::vector<double> ref_helmvec2 = {8.8150, 9.4430, 10.3500, 11.1410, 11.2890, 11.4246, 12.1460, 12.5894, 12.8237, 12.9987};
    const std::vector<double> ref_hayata = {8.8093, 9.3896, 10.2752, 11.1030, 11.2677, 11.4501, 11.9882, 12.6686, 12.8092, 12.9575};

    std::cout << "[2.2.3] wave-number from given beta (Figure 11)\n";
    std::cout << "L=" << L << " beta=" << beta << " beta*L=" << beta * L
              << " nx=" << nx << " ny=" << ny
              << " tris=" << mesh.tris.size()
              << " k0_min_phys~" << k0_min_phys
              << " backend=" << element_assembly_backend_name(cli.backend) << "\n\n";

    std::cout << "first raw roots (k0L), after physical filter:\n";
    for (int i = 0; i < (int)k0_phys.size() && i < 14; ++i)
        std::cout << " " << (i + 1) << "  " << k0_phys[i] << "\n";
    std::cout << "\n";

    std::cout << "mode  k0L(FEM matched)  HELMVEC2(ref)  Hayata(ref)\n";
    std::vector<char> used(k0_phys.size(), 0);
    const int N = std::min<int>(10, ref_helmvec2.size());
    for (int i = 0; i < N; ++i)
    {
        // Modo de validacao: escolhe a raiz calculada nao usada mais proxima da
        // familia modal publicada do HELMVEC2 (Tabela 8), preservando ordenacao.
        const double k0L = helmvec23::pick_closest_unused(ref_helmvec2[i], k0_phys, used);
        std::cout << (i + 1) << "  " << k0L
                  << "  " << ref_helmvec2[i]
                  << "  " << ref_hayata[i] << "\n";
    }

    if (cli.debug_candidates)
    {
        std::cout << "\n[debug] candidate modes after k0_min filter:\n";
        std::cout << "k0L  ez_ratio\n";
        for (const auto &m : all_modes)
        {
            if (m.k0 <= k0_min_phys)
                continue;
            std::cout << m.k0 * L << "  " << m.ez_ratio << "\n";
        }
    }

    perf.total_ms = total_watch.elapsed_ms();
    timing::print_breakdown("helmvec2_rect", perf);
    return 0;
}
