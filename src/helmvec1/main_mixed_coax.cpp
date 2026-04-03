/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/main_mixed_coax.cpp                                  */
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
#include "core/mesh2d_coax.hpp"
#include "core/timing_utils.hpp"
#include "helmvec1_mixed_system.hpp"
#include "mixed_cli_options.hpp"
#include "mixed_debug.hpp"
#include "mixed_mode_utils.hpp"
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
    // Caso coaxial de referencia para a secao 2.2.2 (blocos escalar+vetorial homogeneos).
    const double r1 = 1.0;
    const double r2 = 4.0;
    int nr = 10, nt = 48;
    helmvec1::MixedCliOptions cli;
    try
    {
        cli = helmvec1::parse_mixed_cli_options(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro ao interpretar argumentos: " << e.what() << "\n";
        std::cerr << "Uso: ./mixed_coax [nr nt] [--backend closed-form|gauss]"
                  << " [--debug-local-blocks] [--debug-candidates]\n";
        return 2;
    }

    if (cli.positionals.size() >= 2)
    {
        nr = std::atoi(cli.positionals[0].c_str());
        nt = std::atoi(cli.positionals[1].c_str());
    }

    Mesh2D mesh = make_coax_mesh(r1, r2, nr, nt);
    std::cout << "Mixed coax: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size()
              << " r1=" << r1 << " r2=" << r2
              << " nr=" << nr << " nt=" << nt
              << " backend=" << element_assembly_backend_name(cli.backend) << "\n";

    std::vector<double> eps(mesh.tris.size(), 1.0);
    std::vector<double> mu(mesh.tris.size(), 1.0);
    const double rho = 1e-10;

    if (cli.debug_local_blocks)
        helmvec1_debug::print_first_triangle_closed_form_debug(mesh, 1.0, 1.0);

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

    std::cout << "\n[E-formulation | coax]\n";
    print_first_modes("TE (edge block), first 8 kc:", k_te_edge_e, 8);
    print_first_modes("TM (scalar block), first 8 kc:", k_tm_scalar_e, 8);

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

    std::cout << "\n[H-formulation | coax]\n";
    print_first_modes("TM (edge block), first 8 kc:", k_tm_edge_h, 8);
    print_first_modes("TE (scalar block), first 8 kc:", k_te_scalar_h, 8);

    perf.total_ms = total_watch.elapsed_ms();
    timing::print_breakdown("mixed_coax", perf);
    return 0;
}
