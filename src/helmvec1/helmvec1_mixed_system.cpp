/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/helmvec1_mixed_system.cpp                            */
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

#include "helmvec1_mixed_system.hpp"

namespace
{
/******************************************************************************/
/* FUNCAO: finalize_block_diagonal_system                                     */
/* DESCRICAO: Consolida dimensoes e monta os blocos diagonais globais (S,T)   */
/* do sistema misto vetorial+escalar da Eq. (92).                             */
/* ENTRADA: ms: MixedSystem92 &.                                              */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void finalize_block_diagonal_system(MixedSystem92 &ms)
{
    ms.nt = ms.edge.ed.ndof;
    ms.nz = ms.scal.ndof;
    ms.S = block_diag(ms.edge.S, ms.scal.S);
    ms.T = block_diag(ms.edge.T, ms.scal.T);
}
} // namespace

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
    const std::vector<double> &mu_r_tri)
{
    MixedSystem92 ms;

    // Formulacao E da Secao 2.2.2:
    //   St ~ (1/mu) * curl-curl
    //   Tt ~ eps * massa vetorial
    // com Et tangencial nulo no contorno PEC.
    ms.edge = build_helm10_edge_system(
        mesh,
        EdgeBC::TE_PEC_TangentialZero,
        eps_r_tri,
        mu_r_tri);

    // Bloco longitudinal escalar para Ez com Dirichlet homogenea em PEC.
    ms.scal = build_helm10_scalar_system(
        mesh,
        ScalarBC::TM_Dirichlet,
        eps_r_tri,
        mu_r_tri);

    finalize_block_diagonal_system(ms);
    return ms;
}

/******************************************************************************/
/* FUNCAO: build_system92_H                                                   */
/* DESCRICAO: Monta o sistema misto em blocos da formulacao H para obtencao de*/
/* kc. Corresponde ao sistema em blocos da Eq. (92), Secao 2.2.2.             */
/* ENTRADA: mesh: const Mesh2D &; eps_r_tri: const std::vector<double> &;     */
/* mu_r_tri: const std::vector<double> &.                                     */
/* SAIDA: MixedSystem92.                                                      */
/******************************************************************************/
MixedSystem92 build_system92_H(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri)
{
    MixedSystem92 ms;

    // Formulacao H (dual da formulacao E):
    //   S ~ (1/eps) * curl-curl
    //   T ~ mu * massa vetorial
    // Reuso do mesmo montador trocando os vetores de material:
    //   eps_proxy <- mu_r_tri
    //   mu_proxy  <- eps_r_tri
    ms.edge = build_helm10_edge_system(
        mesh,
        EdgeBC::TM_PEC_NormalZero,
        /*eps_proxy*/ mu_r_tri,
        /*mu_proxy */ eps_r_tri);

    // Bloco escalar dual (Hz) com condicao natural de Neumann.
    ms.scal = build_helm10_scalar_system(
        mesh,
        ScalarBC::TE_Neumann,
        /*eps_proxy*/ mu_r_tri,
        /*mu_proxy */ eps_r_tri);

    finalize_block_diagonal_system(ms);
    return ms;
}
