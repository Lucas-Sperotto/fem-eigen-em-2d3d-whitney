/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/main_mixed_rect.cpp                                  */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema misto vetorial+escalar para kc, separando blocos        */
/* transverso/longitudinal.                                                   */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.2, Eq.    */
/* (92).                                                                      */
/*****************************************************************************/
/* Observacao: Este driver desempenha o papel do programa HELMVEC1 do         */
/* apendice em FORTRAN. A montagem global correspondente a Eq. (92) ocorre em */
/* src/helmvec1/helmvec1_mixed_system.cpp. Comentarios priorizam didatica,    */
/* rastreabilidade e validacao.                                               */
/*****************************************************************************/

#include "article/tp3485_systems.hpp"
#include "core/lapack_eig.hpp"
#include "core/mesh2d_rect.hpp"
#include "core/timing_utils.hpp"
#include "helmvec1_mixed_system.hpp"
#include "mixed_cli_options.hpp"
#include "mixed_debug.hpp"
#include "mixed_mode_utils.hpp"
#include "mixed_rect_reference.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>

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
    // Malha padrao (2*nx*ny triangulos). Para comparacao em formato de tabela,
    // usamos uma discretizacao uniforme moderadamente refinada e expomos (nx, ny).
    int nx = 12, ny = 6;
    helmvec1::MixedCliOptions cli;
    try
    {
        cli = helmvec1::parse_mixed_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./mixed_rect [nx ny] [--backend closed-form|gauss]"
                  << " [--debug-local-blocks] [--debug-candidates]\n";
        return 2;
    }

    if (cli.positionals.size() >= 2)
    {
        nx = std::atoi(cli.positionals[0].c_str());
        ny = std::atoi(cli.positionals[1].c_str());
    }

    const double a = 1.0;
    const double b = 0.5;
    const double rho = 1e-10; // filtro em lambda = kc^2

    Mesh2D mesh = make_rect_mesh(a, b, nx, ny);
    std::vector<double> eps(mesh.tris.size(), 1.0);
    std::vector<double> mu(mesh.tris.size(), 1.0);

    if (cli.debug_local_blocks)
        helmvec1_debug::print_first_triangle_closed_form_debug(mesh, 1.0, 1.0);

    const int nprint = 10;
    const auto ana_te = analytic_rect_te(a, b, nprint);
    const auto ana_tm = analytic_rect_tm(a, b, nprint);

    // Formulacao E: x=[Et; Ez], Sx=(kc^2)Tx com blocos da Eq. (92).
    std::cout << "backend=" << element_assembly_backend_name(cli.backend) << "\n";

    timing::Stopwatch stage;
    auto sys_e = tp3485::build_eq92_helmvec1_system_E(mesh, eps, mu, cli.backend);
    perf.assembly_ms += stage.elapsed_ms();
    stage.reset();
    auto res_e = generalized_eigs_sym_vec(sys_e.S, sys_e.T);
    perf.solve_ms += stage.elapsed_ms();

    std::vector<double> k_te_edge_e, k_tm_scalar_e;
    split_modes_by_block_energy(res_e, sys_e.nt, sys_e.nz, rho, k_te_edge_e, k_tm_scalar_e);
    if (cli.debug_candidates)
    {
        helmvec1_debug::print_block_candidates_debug(
            "[E] candidatos TE (bloco de aresta):",
            k_te_edge_e,
            "[E] candidatos TM (bloco escalar):",
            k_tm_scalar_e);
    }

    print_rect_compare_table("[E-formulation] TE cutoffs (edge block)", ana_te, k_te_edge_e, rho);
    print_rect_compare_table("[E-formulation] TM cutoffs (scalar block)", ana_tm, k_tm_scalar_e, rho);

    // Formulacao H: operador constitutivo dual por troca eps/mu.
    stage.reset();
    auto sys_h = tp3485::build_eq92_helmvec1_system_H(mesh, eps, mu, cli.backend);
    perf.assembly_ms += stage.elapsed_ms();
    stage.reset();
    auto res_h = generalized_eigs_sym_vec(sys_h.S, sys_h.T);
    perf.solve_ms += stage.elapsed_ms();

    std::vector<double> k_tm_edge_h, k_te_scalar_h;
    split_modes_by_block_energy(res_h, sys_h.nt, sys_h.nz, rho, k_tm_edge_h, k_te_scalar_h);
    if (cli.debug_candidates)
    {
        helmvec1_debug::print_block_candidates_debug(
            "[H] candidatos TM (bloco de aresta):",
            k_tm_edge_h,
            "[H] candidatos TE (bloco escalar):",
            k_te_scalar_h);
    }

    // Na formulacao H, o bloco de aresta corresponde a familia TM.
    print_rect_compare_table("[H-formulation] TM cutoffs (edge block)", ana_tm, k_tm_edge_h, rho);
    // Na formulacao H, o bloco escalar corresponde a familia TE.
    print_rect_compare_table("[H-formulation] TE cutoffs (scalar block)", ana_te, k_te_scalar_h, rho);

    perf.total_ms = total_watch.elapsed_ms();
    timing::print_breakdown("mixed_rect", perf);
    return 0;
}
