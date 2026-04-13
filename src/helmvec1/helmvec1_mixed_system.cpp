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
inline ElementAssemblyBackend scalar_backend_for_mixed(ElementAssemblyBackend backend)
{
    return backend;
}

inline ElementAssemblyBackend vector_backend_for_mixed(ElementAssemblyBackend backend)
{
    if (backend == ElementAssemblyBackend::EfgmiConsistent)
        return ElementAssemblyBackend::ClosedForm;
    return backend;
}

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
    ms.edge = build_helm10_edge_system(
        mesh,
        EdgeBC::TE_PEC_TangentialZero,
        eps_r_tri,
        mu_r_tri,
        vector_backend_for_mixed(backend));
    ms.scal = build_helm10_scalar_system(
        mesh,
        ScalarBC::TM_Dirichlet,
        eps_r_tri,
        mu_r_tri,
        scalar_backend_for_mixed(backend));
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
    ms.edge = build_helm10_edge_system(
        mesh,
        EdgeBC::TM_PEC_NormalZero,
        mu_r_tri,
        eps_r_tri,
        vector_backend_for_mixed(backend));
    ms.scal = build_helm10_scalar_system(
        mesh,
        ScalarBC::TE_Neumann,
        mu_r_tri,
        eps_r_tri,
        scalar_backend_for_mixed(backend));
    load_named_eq92_blocks_from_subsystems(ms);
    assemble_eq92_global_from_named_blocks(ms);
    return ms;
}
