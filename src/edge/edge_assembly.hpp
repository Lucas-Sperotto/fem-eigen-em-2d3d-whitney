/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge/edge_assembly.hpp                                        */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo 2D de elementos de aresta (DOFs, base de Whitney e       */
/* montagem).                                                                 */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.1.        */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/dense.hpp"
#include "core/mesh2d.hpp"
#include "edge_dofs.hpp"

// Sistema global de elementos de aresta para o problema generalizado:
//   S x = lambda T x
// com:
//   S_ij = integral_Omega (1/mu_r) curl(W_i) curl(W_j) dOmega
//   T_ij = integral_Omega eps_r (W_i . W_j) dOmega
// (formulacao transversal da Secao 2.2.1 do artigo).
struct EdgeSystem
{
  DenseMat S, T;
  EdgeDofs ed;
};

// Material homogeneo; usa a politica padrao de DOFs de aresta.
/******************************************************************************/
/* FUNCAO: build_helm10_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial transversal de aresta da secao 2.2.1.  */
/* Implementa Sx=lambdaTx da formulacao transversal (Secao 2.2.1).            */
/* ENTRADA: mesh: const Mesh2D &; eps_r: double; mu_r: double.                */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
EdgeSystem build_helm10_edge_system(const Mesh2D &mesh, double eps_r = 1.0, double mu_r = 1.0);

// Material homogeneo com politica explicita de condicao de contorno.
EdgeSystem build_helm10_edge_system(const Mesh2D &mesh, EdgeBC bc, double eps_r = 1.0, double mu_r = 1.0);

// Material nao homogeneo definido por triangulo.
EdgeSystem build_helm10_edge_system(
    const Mesh2D &mesh,
    EdgeBC bc,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri);
