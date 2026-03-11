/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/helmvec1_mixed_system.hpp                            */
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

#pragma once
#include "core/block_ops.hpp"
#include "core/dense.hpp"
#include "core/helm10_scalar_system.hpp"
#include "core/mesh2d.hpp"
#include "edge/edge_assembly.hpp"
#include <vector>

struct MixedSystem92
{
    // Sistema generalizado em blocos (Secao 2.2.2, Eq. 92):
    //   S x = (kc^2) T x, com x = [bloco_transversal_aresta; bloco_escalar_longitudinal]
    DenseMat S; // [St  0]
                // [ 0 Sz]
    DenseMat T; // [Tt  0]
                // [ 0 Tz]

    int nt = 0; // numero de DOFs de aresta (componente transversal)
    int nz = 0; // numero de DOFs escalares (componente longitudinal)

    // Sub-sistemas mantidos para diagnostico e pos-processamento modal.
    EdgeSystem edge;
    ScalarSystem scal;
};

// Formulacao em E (Et, Ez):
//   - BC de aresta: Et tangencial = 0 em PEC
//   - BC escalar: Ez = 0 em PEC
/******************************************************************************/
/* FUNCAO: build_system92_E                                                   */
/* DESCRICAO: Monta o sistema misto em blocos da formulacao E para obtencao de*/
/* kc. Corresponde ao sistema em blocos da Eq. (92), Secao 2.2.2.             */
/* ENTRADA: mesh: const Mesh2D &; eps_r_tri: const std::vector<double> &;     */
/* mu_r_tri: const std::vector<double> &.                                     */
/* SAIDA: MixedSystem92.                                                      */
/******************************************************************************/
MixedSystem92 build_system92_E(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri);

// Formulacao em H (Ht, Hz), operador dual por troca eps/mu:
//   - BC de aresta: mantem arestas de contorno (condicao natural)
//   - BC escalar: Neumann
MixedSystem92 build_system92_H(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri);
