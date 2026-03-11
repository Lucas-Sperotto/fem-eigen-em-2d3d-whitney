/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec2/helmvec2_coupled_system.hpp                          */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema acoplado vetorial+escalar para obter k0 dado beta.      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.3, Eq.    */
/* (108)-(109), Fig. 11, Tabela 8.                                            */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/dense.hpp"
#include "core/helm10_scalar_system.hpp"
#include "core/mesh2d.hpp"
#include "edge/edge_assembly.hpp"
#include <vector>

struct CoupledWaveNumberSystem
{
    // Secao 2.2.3 (obtencao de k0 para beta fixo):
    //   A x = (k0^2) B x      (Eq. 108-109 apos discretizacao)
    // com x = [Et; Ez].
    // Blocos:
    //   A = [S_tt(beta)  S_tz(beta);
    //        S_zt(beta)  S_zz(beta)]
    //   B = [T_tt        0;
    //        0           T_zz(beta)]
    // O operador acoplado em geral nao e simetrico.
    DenseMat A;
    DenseMat B;
    int nt = 0; // numero de DOFs de aresta no bloco Et
    int nz = 0; // numero de DOFs nodais no bloco Ez

    // Sub-blocos mantidos para diagnostico/pos-processamento.
    EdgeSystem edge;
    ScalarSystem scal;
};

struct CoupledBetaSystem
{
    // Secao 2.2.4 (obtencao de beta para k0 fixo):
    //   P x = (beta^2) Q x    (Eq. 126-127 apos rearranjo matricial)
    // com x = [Et; Ez].
    // Blocos:
    //   P = [P_tt(k0)  0;
    //        0         P_zz(k0)]
    //   Q = [Q_tt      Q_tz;
    //        Q_zt      Q_zz]
    DenseMat P;
    DenseMat Q;
    int nt = 0; // numero de DOFs de aresta no bloco Et
    int nz = 0; // numero de DOFs nodais no bloco Ez

    // Sub-blocos mantidos para diagnostico/pos-processamento.
    EdgeSystem edge;
    ScalarSystem scal;
};

// Formulacao em E (Et, Ez), usada nas secoes 2.2.3 e 2.2.4:
//   - BC de aresta: Et tangencial = 0 em PEC
//   - BC escalar: Ez = 0 em PEC
/******************************************************************************/
/* FUNCAO: build_coupled_wavenumber_system_E                                  */
/* DESCRICAO: Monta o sistema acoplado A x = k0^2 B x para k0 dado beta.      */
/* Corresponde ao problema da Secao 2.2.3 (Eq. 108-109), com formulacao em E  */
/* usando os blocos Et/Ez.                                                    */
/* ENTRADA: mesh: const Mesh2D &; beta: double; eps_r_tri: const              */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &.              */
/* SAIDA: CoupledWaveNumberSystem.                                            */
/******************************************************************************/
CoupledWaveNumberSystem build_coupled_wavenumber_system_E(
    const Mesh2D &mesh,
    double beta,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri);

/******************************************************************************/
/* FUNCAO: build_coupled_beta_system_E                                        */
/* DESCRICAO: Monta o sistema acoplado P x = beta^2 Q x para beta dado k0.    */
/* Corresponde ao problema da Secao 2.2.4 (Eq. 126-127), com formulacao em E  */
/* e acoplamento entre Et e Ez.                                               */
/* ENTRADA: mesh: const Mesh2D &; k0: double; eps_r_tri: const                */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &.              */
/* SAIDA: CoupledBetaSystem.                                                  */
/******************************************************************************/
CoupledBetaSystem build_coupled_beta_system_E(
    const Mesh2D &mesh,
    double k0,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri);
