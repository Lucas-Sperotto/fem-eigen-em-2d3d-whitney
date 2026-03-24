/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec/edge_debug.hpp                                        */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios didaticos de depuracao para a formulacao de aresta  */
/* 2D da Secao 2.2.1.                                                         */
/*****************************************************************************/

#pragma once

#include "edge/edge_basis.hpp"
#include "explicit/tri2d_edge_explicit.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace helmvec_debug
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
/* DESCRICAO: Imprime os blocos locais vetoriais do primeiro triangulo,       */
/* ligados as Eq. (66) e (67), que alimentam o sistema global da Eq. (65).    */
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
    const TriGeomEdge tg = tri_geom_edge(mesh, t);
    double Sel[3][3] = {{0.0}};
    double Tel[3][3] = {{0.0}};
    explicit_tri2d::tri2d_edge_closed_form_eq_66_67(tg, 1.0 / mu_r, eps_r, Sel, Tel);

    std::cout << "\n[debug] primeiro triangulo: blocos locais closed-form 2D de aresta\n";
    std::cout << "  area=" << tg.g.A << " eps_r=" << eps_r << " mu_r=" << mu_r << "\n";
    std::cout << "  artigo Eq. (66) e Eq. (67):\n";
    print_block_3x3("Sel", Sel);
    print_block_3x3("Tel", Tel);
}

/******************************************************************************/
/* FUNCAO: print_positive_kc_candidates_debug                                 */
/* DESCRICAO: Imprime as primeiras raizes positivas kc antes do matching      */
/* modal/analitico.                                                           */
/* ENTRADA: lambda: const std::vector<double> &; zero_tol: double; max_print: */
/* int.                                                                       */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_positive_kc_candidates_debug(
    const std::vector<double> &lambda,
    double zero_tol,
    int max_print = 20)
{
    std::cout << "\n[debug] primeiros candidatos kc positivos antes do matching:\n";
    int shown = 0;
    for (double l : lambda)
    {
        if (l <= zero_tol)
            continue;
        std::cout << "  " << (shown + 1) << "  " << std::sqrt(l) << "\n";
        ++shown;
        if (shown >= max_print)
            break;
    }
}

} // namespace helmvec_debug
