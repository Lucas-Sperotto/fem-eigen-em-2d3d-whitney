/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_debug.hpp                                      */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios didaticos de depuracao para a formulacao mista 2D   */
/* da Secao 2.2.2.                                                            */
/*****************************************************************************/

#pragma once

#include "core/fem_scalar.hpp"
#include "edge/edge_basis.hpp"
#include "explicit/tri2d_edge_explicit.hpp"
#include "explicit/tri2d_scalar_explicit.hpp"
#include "mixed_mode_utils.hpp"

#include <iostream>
#include <vector>

namespace helmvec1_debug
{

/******************************************************************************/
/* FUNCAO: print_block_3x3                                                    */
/* DESCRICAO: Imprime uma matriz 3x3 em formato compacto para depuracao local.*/
/* ENTRADA: name: const char *; M: const double[3][3].                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_block_3x3(const char *name, const double M[3][3])
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
/* DESCRICAO: Imprime os blocos locais do primeiro triangulo para a           */
/* formulacao mista da Eq. (92): bloco de aresta, bloco escalar e a ideia de  */
/* arranjo block-diagonal no cutoff.                                          */
/* ENTRADA: mesh: const Mesh2D &; eps_r: double; mu_r: double.                */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_first_triangle_closed_form_debug(
    const Mesh2D &mesh,
    double eps_r,
    double mu_r)
{
    if (mesh.tris.empty())
        return;

    const Tri &t = mesh.tris.front();
    const TriGeom g = tri_geom(mesh, t);
    const TriGeomEdge tg = tri_geom_edge(mesh, t);
    double Sel_edge[3][3] = {{0.0}};
    double Tel_edge[3][3] = {{0.0}};
    double Sel_scal[3][3] = {{0.0}};
    double Tel_scal[3][3] = {{0.0}};
    explicit_tri2d::tri2d_edge_closed_form_eq_66_67(tg, 1.0 / mu_r, eps_r, Sel_edge, Tel_edge);
    explicit_tri2d::tri2d_scalar_closed_form_eq_30_33(g.A, g.b, g.c, Sel_scal, Tel_scal);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            Sel_scal[i][j] *= (1.0 / mu_r);
            Tel_scal[i][j] *= eps_r;
        }

    std::cout << "\n[debug] primeiro triangulo: blocos locais closed-form 2D mistos\n";
    std::cout << "  area=" << g.A << " eps_r=" << eps_r << " mu_r=" << mu_r << "\n";
    std::cout << "  blocos locais que alimentam a Eq. (92):\n";
    print_block_3x3("St(edge)", Sel_edge);
    print_block_3x3("Tt(edge)", Tel_edge);
    print_block_3x3("Sz(scalar)", Sel_scal);
    print_block_3x3("Tz(scalar)", Tel_scal);
}

/******************************************************************************/
/* FUNCAO: print_block_candidates_debug                                       */
/* DESCRICAO: Imprime listas compactas de candidatos kc por bloco dominante   */
/* apos a separacao energetica do sistema misto.                              */
/* ENTRADA: title0: const char *; vals0: const std::vector<double> &; title1: */
/* const char *; vals1: const std::vector<double> &.                          */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_block_candidates_debug(
    const char *title0,
    const std::vector<double> &vals0,
    const char *title1,
    const std::vector<double> &vals1)
{
    std::cout << "\n[debug] candidatos por bloco dominante:\n";
    print_first_modes(title0, vals0, 20);
    print_first_modes(title1, vals1, 20);
}

} // namespace helmvec1_debug
