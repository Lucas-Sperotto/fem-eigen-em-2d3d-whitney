/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/main_mixed_coax.cpp                                  */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema misto vetorial+escalar para kc, separando blocos        */
/* transverso/longitudinal.                                                   */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.2, Eq.    */
/* (92).                                                                      */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "core/lapack_eig.hpp"
#include "core/mesh2d_coax.hpp"
#include "helmvec1_mixed_system.hpp"
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
    // Caso coaxial de referencia para a secao 2.2.2 (blocos escalar+vetorial homogeneos).
    const double r1 = 1.0;
    const double r2 = 4.0;
    int nr = 10, nt = 48;
    if (argc >= 3)
    {
        nr = std::atoi(argv[1]);
        nt = std::atoi(argv[2]);
    }

    Mesh2D mesh = make_coax_mesh(r1, r2, nr, nt);
    std::cout << "Mixed coax: nodes=" << mesh.nodes.size()
              << " tris=" << mesh.tris.size()
              << " r1=" << r1 << " r2=" << r2
              << " nr=" << nr << " nt=" << nt << "\n";

    std::vector<double> eps(mesh.tris.size(), 1.0);
    std::vector<double> mu(mesh.tris.size(), 1.0);
    const double rho = 1e-10;

    auto sys_e = build_system92_E(mesh, eps, mu);
    auto res_e = generalized_eigs_sym_vec(sys_e.S, sys_e.T);
    std::vector<double> k_te_edge_e, k_tm_scalar_e;
    split_modes_by_block_energy(res_e, sys_e.nt, sys_e.nz, rho, k_te_edge_e, k_tm_scalar_e);

    std::cout << "\n[E-formulation | coax]\n";
    print_first_modes("TE (edge block), first 8 kc:", k_te_edge_e, 8);
    print_first_modes("TM (scalar block), first 8 kc:", k_tm_scalar_e, 8);

    auto sys_h = build_system92_H(mesh, eps, mu);
    auto res_h = generalized_eigs_sym_vec(sys_h.S, sys_h.T);
    std::vector<double> k_tm_edge_h, k_te_scalar_h;
    split_modes_by_block_energy(res_h, sys_h.nt, sys_h.nz, rho, k_tm_edge_h, k_te_scalar_h);

    std::cout << "\n[H-formulation | coax]\n";
    print_first_modes("TM (edge block), first 8 kc:", k_tm_edge_h, 8);
    print_first_modes("TE (scalar block), first 8 kc:", k_te_scalar_h, 8);

    return 0;
}
