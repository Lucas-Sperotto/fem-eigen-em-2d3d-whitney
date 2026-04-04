/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/field_reconstruction.hpp                              */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Reconstrucao didatica dos campos transversais a partir do       */
/* potencial longitudinal escalar da Secao 2.1.                               */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Eq. (35)-(41).     */
/*****************************************************************************/
/* Observacao: Este modulo deixa explicita a interpretacao fisica da          */
/* variavel escalar do HELM10:                                                */
/* - TE: potencial longitudinal associado a Hz;                               */
/* - TM: potencial longitudinal associado a Ez.                               */
/* A partir desse potencial, o gradiente transversal e convertido em Ex e Ey. */
/*****************************************************************************/

#pragma once

#include "core/mesh2d.hpp"

#include <string>
#include <vector>

namespace helm10::field_reconstruction
{

enum class LongitudinalScalarKind
{
    TE_Hz,
    TM_Ez
};

struct HomogeneousMedium
{
    double frequency_hz = 0.0;
    double omega_rad_s = 0.0;
    double eps_r = 1.0;
    double mu_r = 1.0;
    double epsilon = 0.0;
    double mu = 0.0;
    double k = 0.0;
};

struct ReconstructedField2D
{
    LongitudinalScalarKind kind = LongitudinalScalarKind::TE_Hz;
    std::vector<double> psi;
    std::vector<double> dpsi_dx;
    std::vector<double> dpsi_dy;
    std::vector<double> ex;
    std::vector<double> ey;
    std::vector<double> ex_without_ztm;
    std::vector<double> ey_without_ztm;
    std::vector<double> ex_with_ztm;
    std::vector<double> ey_with_ztm;
    double kc = 0.0;
    double k = 0.0;
    double beta = 0.0;
    double ztm = 0.0;
    bool below_cutoff = false;
    bool tm_has_scaled_fields = false;
    std::string longitudinal_label;
    std::string field_equations;
    std::string status_label;
    std::string status_message;
};

/******************************************************************************/
/* FUNCAO: vacuum_permittivity                                                */
/* DESCRICAO: Retorna a permissividade eletrica do vacuo em SI.               */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: double.                                                             */
/******************************************************************************/
double vacuum_permittivity();

/******************************************************************************/
/* FUNCAO: vacuum_permeability                                                */
/* DESCRICAO: Retorna a permeabilidade magnetica do vacuo em SI.              */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: double.                                                             */
/******************************************************************************/
double vacuum_permeability();

/******************************************************************************/
/* FUNCAO: make_homogeneous_medium_from_relative                              */
/* DESCRICAO: Constrói os parametros fisicos de um meio homogeneo a partir de */
/* frequencia e constantes relativas, para uso na reconstrucao TM.            */
/* ENTRADA: frequency_hz: double; eps_r: double; mu_r: double.                */
/* SAIDA: HomogeneousMedium.                                                  */
/******************************************************************************/
HomogeneousMedium make_homogeneous_medium_from_relative(
    double frequency_hz,
    double eps_r = 1.0,
    double mu_r = 1.0);

/******************************************************************************/
/* FUNCAO: make_homogeneous_medium_from_k                                     */
/* DESCRICAO: Constrói o meio homogeneo a partir do numero de onda do meio,   */
/* deixando a frequencia consistente com `k = omega * sqrt(mu * eps)`.        */
/* ENTRADA: k_medium: double; eps_r: double; mu_r: double.                    */
/* SAIDA: HomogeneousMedium.                                                  */
/******************************************************************************/
HomogeneousMedium make_homogeneous_medium_from_k(
    double k_medium,
    double eps_r = 1.0,
    double mu_r = 1.0);

/******************************************************************************/
/* FUNCAO: make_homogeneous_medium_above_kc                                   */
/* DESCRICAO: Escolhe automaticamente uma frequencia tal que o numero de onda */
/* do meio satisfaca k > kc_max, garantindo propagacao para os modos          */
/* solicitados com uma pequena margem de seguranca.                           */
/* ENTRADA: kc_max: double; eps_r: double; mu_r: double; safety_factor: double.*/
/* SAIDA: HomogeneousMedium.                                                  */
/******************************************************************************/
HomogeneousMedium make_homogeneous_medium_above_kc(
    double kc_max,
    double eps_r = 1.0,
    double mu_r = 1.0,
    double safety_factor = 1.05);

/******************************************************************************/
/* FUNCAO: make_homogeneous_medium_for_tm_reference_ztm                       */
/* DESCRICAO: Escolhe automaticamente um meio cuja frequencia satisfaz        */
/* `Ztm = ztm_target` para um modo TM de referencia com cutoff `kc_ref`.      */
/* Isso e util para uma politica didatica em que o modo TM limitante fixa a   */
/* escala de impedancia e os demais modos passam a ser comparados a partir    */
/* dela.                                                                      */
/* ENTRADA: kc_ref: double; ztm_target: double; eps_r: double; mu_r: double.  */
/* SAIDA: HomogeneousMedium.                                                  */
/******************************************************************************/
HomogeneousMedium make_homogeneous_medium_for_tm_reference_ztm(
    double kc_ref,
    double ztm_target = 1.0,
    double eps_r = 1.0,
    double mu_r = 1.0);

/******************************************************************************/
/* FUNCAO: longitudinal_label                                                 */
/* DESCRICAO: Retorna o rotulo fisico da incognita longitudinal escalar.      */
/* ENTRADA: kind: LongitudinalScalarKind.                                     */
/* SAIDA: const char *.                                                       */
/******************************************************************************/
const char *longitudinal_label(LongitudinalScalarKind kind);

/******************************************************************************/
/* FUNCAO: field_equation_label                                               */
/* DESCRICAO: Retorna a faixa de equacoes do artigo usada na reconstrucao     */
/* transversal. No caso TM, a saida padrao do projeto salva o campo apenas    */
/* com gradiente e sinal, sem multiplicar por Ztm.                            */
/* ENTRADA: kind: LongitudinalScalarKind.                                     */
/* SAIDA: const char *.                                                       */
/******************************************************************************/
const char *field_equation_label(LongitudinalScalarKind kind);

/******************************************************************************/
/* FUNCAO: reconstruct_transverse_fields                                      */
/* DESCRICAO: Reconstrói dpsi/dx, dpsi/dy, Ex e Ey a partir do potencial      */
/* longitudinal escalar. Para TE usa as Eq. (38)-(39). Para TM, o padrao      */
/* didatico atual salva Ex e Ey apenas com gradiente e sinal, sem aplicar     */
/* o fator Ztm; beta e Ztm permanecem disponiveis como informacao auxiliar.   */
/* ENTRADA: mesh: const Mesh2D &; psi_nodal: const std::vector<double> &;     */
/* kind: LongitudinalScalarKind; kc: double; medium: const HomogeneousMedium &;*/
/* normalize_transverse_field: bool; area_tolerance: double.                  */
/* SAIDA: ReconstructedField2D.                                               */
/******************************************************************************/
ReconstructedField2D reconstruct_transverse_fields(
    const Mesh2D &mesh,
    const std::vector<double> &psi_nodal,
    LongitudinalScalarKind kind,
    double kc,
    const HomogeneousMedium &medium,
    bool normalize_transverse_field = false,
    double area_tolerance = 1e-14);

} // namespace helm10::field_reconstruction
