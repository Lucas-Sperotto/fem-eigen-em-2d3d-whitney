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
#include "core/assembly_backend.hpp"
#include "core/dense.hpp"
#include "core/helm10_scalar_system.hpp"
#include "core/mesh2d.hpp"
#include "edge/edge_assembly.hpp"
#include <vector>

struct CoupledWaveNumberSystem
{
    // Secao 2.2.3 (obtencao de k0 para beta fixo):
    //   A x = (k0^2) B x
    // Na leitura do artigo, esta estrutura representa o sistema global montado
    // da Eq. (119), obtido a partir da discretizacao das Eq. (108)-(109) e da
    // assembleia dos blocos locais Eq. (120)-(125).
    // com x = [Et; Ez].
    // Blocos:
    //   A = [S_tt(beta)  S_tz(beta);
    //        S_zt(beta)  S_zz(beta)]
    //   B = [T_tt        0;
    //        0           T_zz(beta)]
    // O operador acoplado em geral nao e simetrico.
    // Este bloco corresponde ao papel do programa HELMVEC2 no apendice
    // em FORTRAN.
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
    //   P x = (beta^2) Q x
    // Na leitura do artigo, esta estrutura representa o sistema global montado
    // da Eq. (136), obtido apos o rearranjo das Eq. (126)-(127) e a assembleia
    // dos blocos locais Eq. (137)-(142).
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
    // Este bloco corresponde ao papel do programa HELMVEC3 no apendice
    // em FORTRAN.
    EdgeSystem edge;
    ScalarSystem scal;
};

// Formulacao em E (Et, Ez), usada nas secoes 2.2.3 e 2.2.4:
//   - BC de aresta: Et tangencial = 0 em PEC
//   - BC escalar: Ez = 0 em PEC
/******************************************************************************/
/* FUNCAO: build_coupled_wavenumber_system_E                                  */
/* DESCRICAO: Monta o sistema acoplado A x = k0^2 B x para k0 dado beta.      */
/* Corresponde ao problema da Secao 2.2.3. No nivel de sistema global, este   */
/* e o ponto onde o codigo monta a Eq. (119), usando formulacao em E com      */
/* blocos Et/Ez. O parametro backend permite escolher entre a montagem por     */
/* quadratura e a versao closed-form local das Eq. (120)-(125).               */
/* ENTRADA: mesh: const Mesh2D &; beta: double; eps_r_tri: const              */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: CoupledWaveNumberSystem.                                            */
/******************************************************************************/
CoupledWaveNumberSystem build_coupled_wavenumber_system_E(
    const Mesh2D &mesh,
    double beta,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend = ElementAssemblyBackend::GaussianQuadrature);

/******************************************************************************/
/* FUNCAO: build_coupled_beta_system_E                                        */
/* DESCRICAO: Monta o sistema acoplado P x = beta^2 Q x para beta dado k0.    */
/* Corresponde ao problema da Secao 2.2.4. No nivel de sistema global, este   */
/* e o ponto onde o codigo monta a Eq. (136), com formulacao em E e           */
/* acoplamento entre Et e Ez. O parametro backend permite escolher entre a    */
/* montagem por quadratura e os blocos constituintes closed-form.             */
/* ENTRADA: mesh: const Mesh2D &; k0: double; eps_r_tri: const                */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: CoupledBetaSystem.                                                  */
/******************************************************************************/
CoupledBetaSystem build_coupled_beta_system_E(
    const Mesh2D &mesh,
    double k0,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend = ElementAssemblyBackend::GaussianQuadrature);
