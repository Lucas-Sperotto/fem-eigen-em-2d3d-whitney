/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/helmvec1_mixed_system.cpp                            */
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
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "helmvec1_mixed_system.hpp"

namespace
{
/******************************************************************************/
/* FUNCAO: load_named_eq92_blocks_from_subsystems                             */
/* DESCRICAO: Copia, com nomes explicitos da Eq. (92), os blocos locais ja    */
/* montados dos subproblemas de aresta e escalar para dentro de MixedSystem92.*/
/* ENTRADA: ms: MixedSystem92 &.                                              */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void load_named_eq92_blocks_from_subsystems(MixedSystem92 &ms)
{
    ms.nt = ms.edge.ed.ndof;
    ms.nz = ms.scal.ndof;
    ms.St = ms.edge.S;
    ms.Tt = ms.edge.T;
    ms.Sz = ms.scal.S;
    ms.Tz = ms.scal.T;
}

/******************************************************************************/
/* FUNCAO: assemble_eq92_global_from_named_blocks                             */
/* DESCRICAO: Monta explicitamente a Eq. (92) a partir dos blocos nomeados    */
/* St, Tt, Sz e Tz, preservando a leitura didatica do sistema block-diagonal. */
/* ENTRADA: ms: MixedSystem92 &.                                              */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void assemble_eq92_global_from_named_blocks(MixedSystem92 &ms)
{
    ms.S = block_diag(ms.St, ms.Sz);
    ms.T = block_diag(ms.Tt, ms.Tz);
}
} // namespace

/******************************************************************************/
/* FUNCAO: build_system92_E                                                   */
/* DESCRICAO: Monta o sistema misto em blocos da formulacao E para obtencao de*/
/* kc. Corresponde ao sistema em blocos da Eq. (92), Secao 2.2.2.             */
/* ENTRADA: mesh: const Mesh2D &; eps_r_tri: const std::vector<double> &;     */
/* mu_r_tri: const std::vector<double> &; backend:                            */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: MixedSystem92.                                                      */
/******************************************************************************/
MixedSystem92 build_system92_E(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend)
{
    MixedSystem92 ms;

    // Formulacao E da Secao 2.2.2:
    //   Eq. (88): Sel(t) = (1/mu_r) * curl-curl
    //            (na pratica, reaproveita Eq. (66) no bloco edge)
    //   Eq. (90): Tel(t) = eps_r * massa vetorial
    //            (na pratica, reaproveita Eq. (67) no bloco edge)
    // com Et tangencial nulo no contorno PEC.
    ms.edge = build_helm10_edge_system(
        mesh,
        EdgeBC::TE_PEC_TangentialZero,
        eps_r_tri,
        mu_r_tri,
        backend);

    // Bloco longitudinal escalar para Ez com Dirichlet homogenea em PEC:
    //   Eq. (89): Sel(z) = (1/mu_r) * grad-grad
    //   Eq. (91): Tel(z) = eps_r * massa escalar
    // com reaproveitamento das Eq. (31) e (33) do bloco escalar.
    ms.scal = build_helm10_scalar_system(
        mesh,
        ScalarBC::TM_Dirichlet,
        eps_r_tri,
        mu_r_tri,
        backend);

    load_named_eq92_blocks_from_subsystems(ms);
    assemble_eq92_global_from_named_blocks(ms);
    return ms;
}

/******************************************************************************/
/* FUNCAO: build_system92_H                                                   */
/* DESCRICAO: Monta o sistema misto em blocos da formulacao H para obtencao de*/
/* kc. Corresponde ao sistema em blocos da Eq. (92), Secao 2.2.2.             */
/* ENTRADA: mesh: const Mesh2D &; eps_r_tri: const std::vector<double> &;     */
/* mu_r_tri: const std::vector<double> &; backend:                            */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: MixedSystem92.                                                      */
/******************************************************************************/
MixedSystem92 build_system92_H(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend)
{
    MixedSystem92 ms;

    // Formulacao H (dual da formulacao E):
    //   Eq. (88)-(91) com troca constitutiva E<->H:
    //   S ~ (1/eps_r) * curl-curl / grad-grad
    //   T ~ mu_r * massa vetorial / massa escalar
    // Reuso do mesmo montador trocando os vetores de material:
    //   eps_proxy <- mu_r_tri
    //   mu_proxy  <- eps_r_tri
    ms.edge = build_helm10_edge_system(
        mesh,
        EdgeBC::TM_PEC_NormalZero,
        /*eps_proxy*/ mu_r_tri,
        /*mu_proxy */ eps_r_tri,
        backend);

    // Bloco escalar dual (Hz) com condicao natural de Neumann.
    ms.scal = build_helm10_scalar_system(
        mesh,
        ScalarBC::TE_Neumann,
        /*eps_proxy*/ mu_r_tri,
        /*mu_proxy */ eps_r_tri,
        backend);

    load_named_eq92_blocks_from_subsystems(ms);
    assemble_eq92_global_from_named_blocks(ms);
    return ms;
}
