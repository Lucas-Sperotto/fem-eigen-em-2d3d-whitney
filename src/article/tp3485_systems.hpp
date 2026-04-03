/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/article/tp3485_systems.hpp                                    */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Fachada didatica nomeada pelas equacoes globais do artigo.      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Eq. (43), (65),     */
/* (92), (119), (136) e (178).                                                */
/*****************************************************************************/
/* Observacao: Esta camada existe para tornar a rastreabilidade matematica    */
/* explicita no proprio ponto de entrada do codigo C++.                       */
/*****************************************************************************/

#pragma once

#include "core/assembly_backend.hpp"
#include "core/helm10_scalar_system.hpp"
#include "edge/edge_assembly.hpp"
#include "edge3d/edge3d_assembly.hpp"
#include "helmvec1/helmvec1_mixed_system.hpp"
#include "helmvec2/helmvec2_coupled_system.hpp"
#include <vector>

namespace tp3485
{

/******************************************************************************/
/* FUNCAO: build_eq43_helm10_system                                           */
/* DESCRICAO: Fachada didatica da Eq. (43), sistema global escalar associado  */
/* ao programa HELM10 do artigo.                                              */
/* ENTRADA: mesh: const Mesh2D &; bc: ScalarBC; backend:                      */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: ScalarSystem.                                                       */
/******************************************************************************/
inline ScalarSystem build_eq43_helm10_system(
    const Mesh2D &mesh,
    ScalarBC bc,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_helm10_scalar_system(mesh, bc, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq43_helm10_system                                           */
/* DESCRICAO: Variante heterogenea da Eq. (43), com materiais por triangulo.  */
/* ENTRADA: mesh: const Mesh2D &; bc: ScalarBC; eps_r_tri: const              */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: ScalarSystem.                                                       */
/******************************************************************************/
inline ScalarSystem build_eq43_helm10_system(
    const Mesh2D &mesh,
    ScalarBC bc,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_helm10_scalar_system(mesh, bc, eps_r_tri, mu_r_tri, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq65_helmvec_system                                          */
/* DESCRICAO: Fachada didatica da Eq. (65), sistema global vetorial           */
/* transversal associado ao programa HELMVEC do artigo.                       */
/* ENTRADA: mesh: const Mesh2D &; bc: EdgeBC; eps_r_tri: const                */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
inline EdgeSystem build_eq65_helmvec_system(
    const Mesh2D &mesh,
    EdgeBC bc,
    double eps_r = 1.0,
    double mu_r = 1.0,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_helm10_edge_system(mesh, bc, eps_r, mu_r, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq65_helmvec_system                                          */
/* DESCRICAO: Variante heterogenea da Eq. (65), com materiais por triangulo.  */
/* ENTRADA: mesh: const Mesh2D &; bc: EdgeBC; eps_r_tri: const                */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
inline EdgeSystem build_eq65_helmvec_system(
    const Mesh2D &mesh,
    EdgeBC bc,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_helm10_edge_system(mesh, bc, eps_r_tri, mu_r_tri, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq92_helmvec1_system_E                                       */
/* DESCRICAO: Fachada didatica da Eq. (92) na formulacao E, associada ao      */
/* programa HELMVEC1 no cutoff.                                               */
/* ENTRADA: mesh: const Mesh2D &; eps_r_tri: const std::vector<double> &;     */
/* mu_r_tri: const std::vector<double> &; backend: ElementAssemblyBackend.    */
/* SAIDA: MixedSystem92.                                                      */
/******************************************************************************/
inline MixedSystem92 build_eq92_helmvec1_system_E(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_system92_E(mesh, eps_r_tri, mu_r_tri, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq92_helmvec1_system_H                                       */
/* DESCRICAO: Fachada didatica da Eq. (92) na formulacao H, associada ao      */
/* programa HELMVEC1 no cutoff.                                               */
/* ENTRADA: mesh: const Mesh2D &; eps_r_tri: const std::vector<double> &;     */
/* mu_r_tri: const std::vector<double> &; backend: ElementAssemblyBackend.    */
/* SAIDA: MixedSystem92.                                                      */
/******************************************************************************/
inline MixedSystem92 build_eq92_helmvec1_system_H(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_system92_H(mesh, eps_r_tri, mu_r_tri, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq119_helmvec2_system_E                                      */
/* DESCRICAO: Fachada didatica da Eq. (119), sistema acoplado em k0^2 da      */
/* Secao 2.2.3, associado ao programa HELMVEC2.                               */
/* ENTRADA: mesh: const Mesh2D &; beta: double; eps_r_tri: const              */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: CoupledWaveNumberSystem.                                            */
/******************************************************************************/
inline CoupledWaveNumberSystem build_eq119_helmvec2_system_E(
    const Mesh2D &mesh,
    double beta,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_coupled_wavenumber_system_E(mesh, beta, eps_r_tri, mu_r_tri, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq136_helmvec3_system_E                                      */
/* DESCRICAO: Fachada didatica da Eq. (136), sistema acoplado em beta^2 da    */
/* Secao 2.2.4, associado ao programa HELMVEC3.                               */
/* ENTRADA: mesh: const Mesh2D &; k0: double; eps_r_tri: const                */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: CoupledBetaSystem.                                                  */
/******************************************************************************/
inline CoupledBetaSystem build_eq136_helmvec3_system_E(
    const Mesh2D &mesh,
    double k0,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_coupled_beta_system_E(mesh, k0, eps_r_tri, mu_r_tri, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq178_fem3d_system_dense                                     */
/* DESCRICAO: Fachada didatica da Eq. (178) no fluxo denso, correspondente ao */
/* papel historico do FEM3D0.                                                 */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r_tet: const              */
/* std::vector<double> &; mu_r_tet: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: EdgeSystem3D.                                                       */
/******************************************************************************/
inline EdgeSystem3D build_eq178_fem3d_system_dense(
    const Mesh3D &mesh,
    Edge3DBC bc,
    double eps_r = 1.0,
    double mu_r = 1.0,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_helm3d_edge_system(mesh, bc, eps_r, mu_r, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq178_fem3d_system_dense                                     */
/* DESCRICAO: Variante heterogenea da Eq. (178) no fluxo denso.               */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r_tet: const              */
/* std::vector<double> &; mu_r_tet: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: EdgeSystem3D.                                                       */
/******************************************************************************/
inline EdgeSystem3D build_eq178_fem3d_system_dense(
    const Mesh3D &mesh,
    Edge3DBC bc,
    const std::vector<double> &eps_r_tet,
    const std::vector<double> &mu_r_tet,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_helm3d_edge_system(mesh, bc, eps_r_tet, mu_r_tet, backend);
}

/******************************************************************************/
/* FUNCAO: build_eq178_fem3d_system_sparse                                    */
/* DESCRICAO: Fachada didatica da Eq. (178) no fluxo esparso, correspondente  */
/* ao papel historico do FEM3D1.                                              */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r_tet: const              */
/* std::vector<double> &; mu_r_tet: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: EdgeSystem3DSparse.                                                 */
/******************************************************************************/
inline EdgeSystem3DSparse build_eq178_fem3d_system_sparse(
    const Mesh3D &mesh,
    Edge3DBC bc,
    double eps_r = 1.0,
    double mu_r = 1.0,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_helm3d_edge_system_sparse(
        mesh,
        bc,
        std::vector<double>(mesh.tets.size(), eps_r),
        std::vector<double>(mesh.tets.size(), mu_r),
        backend);
}

/******************************************************************************/
/* FUNCAO: build_eq178_fem3d_system_sparse                                    */
/* DESCRICAO: Variante heterogenea da Eq. (178) no fluxo esparso.             */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r_tet: const              */
/* std::vector<double> &; mu_r_tet: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: EdgeSystem3DSparse.                                                 */
/******************************************************************************/
inline EdgeSystem3DSparse build_eq178_fem3d_system_sparse(
    const Mesh3D &mesh,
    Edge3DBC bc,
    const std::vector<double> &eps_r_tet,
    const std::vector<double> &mu_r_tet,
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm)
{
    return build_helm3d_edge_system_sparse(mesh, bc, eps_r_tet, mu_r_tet, backend);
}

} // namespace tp3485
