/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge3d/edge3d_assembly.hpp                                    */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo 3D de elementos de aresta tetraedricos (Whitney 1-form) e*/
/* montagem.                                                                  */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 3.1, integrais*/
/* I1..I10.                                                                   */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/dense.hpp"
#include "core/mesh3d.hpp"
#include "core/sparse_sym.hpp"
#include "edge3d_dofs.hpp"
#include <vector>

struct EdgeSystem3D
{
  // Problema generalizado 3D:
  //   S x = lambda T x
  // com:
  //   S_ij = int_Omega (1/mu_r) curl(W_i).curl(W_j) dV
  //   T_ij = int_Omega eps_r (W_i.W_j) dV
  // (Secao 3.1; integrandos associados aos termos I1..I10).
  DenseMat S, T;
  EdgeDofs3D ed;
};

struct EdgeSystem3DSparse
{
  // Mesmo problema Sx=lambdaTx, mas armazenado no formato esparso simetrico
  // (apenas triangular inferior) para reduzir memoria em malhas grandes.
  SparseSymMat S, T;
  EdgeDofs3D ed;
};

/******************************************************************************/
/* FUNCAO: build_helm3d_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial 3D de aresta para autovalores de       */
/* cavidades. Implementa a formulacao vetorial 3D da Secao 3.1.               */
/* ENTRADA: mesh: const Mesh3D &; eps_r: double; mu_r: double.                */
/* SAIDA: EdgeSystem3D.                                                       */
/******************************************************************************/
EdgeSystem3D build_helm3d_edge_system(
    const Mesh3D &mesh,
    double eps_r = 1.0,
    double mu_r = 1.0);

/******************************************************************************/
/* FUNCAO: build_helm3d_edge_system                                           */
/* DESCRICAO: Variante homogenea com configuracao explicita de BC de aresta.  */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r: double; mu_r: double.  */
/* SAIDA: EdgeSystem3D.                                                       */
/******************************************************************************/
EdgeSystem3D build_helm3d_edge_system(
    const Mesh3D &mesh,
    Edge3DBC bc,
    double eps_r = 1.0,
    double mu_r = 1.0);

/******************************************************************************/
/* FUNCAO: build_helm3d_edge_system                                           */
/* DESCRICAO: Variante geral heterogenea por tetraedro para eps_r e mu_r.     */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r_tet: const              */
/* std::vector<double> &; mu_r_tet: const std::vector<double> &.              */
/* SAIDA: EdgeSystem3D.                                                       */
/******************************************************************************/
EdgeSystem3D build_helm3d_edge_system(
    const Mesh3D &mesh,
    Edge3DBC bc,
    const std::vector<double> &eps_r_tet,
    const std::vector<double> &mu_r_tet);

/******************************************************************************/
/* FUNCAO: build_helm3d_edge_system_sparse                                    */
/* DESCRICAO: Monta o sistema vetorial 3D em formato esparso simetrico        */
/* (triangulo inferior). Implementa a formulacao vetorial 3D da Secao 3.1.    */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r_tet: const              */
/* std::vector<double> &; mu_r_tet: const std::vector<double> &.              */
/* SAIDA: EdgeSystem3DSparse.                                                 */
/******************************************************************************/
EdgeSystem3DSparse build_helm3d_edge_system_sparse(
    const Mesh3D &mesh,
    Edge3DBC bc,
    const std::vector<double> &eps_r_tet,
    const std::vector<double> &mu_r_tet);
