/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/field_reconstruction.cpp                              */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Implementacao da reconstrucao didatica dos campos transversais  */
/* do HELM10 a partir do potencial longitudinal escalar.                      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Eq. (35)-(41).     */
/*****************************************************************************/
/* Observacao: O objetivo aqui nao e apenas produzir numeros, mas mostrar     */
/* explicitamente a passagem: potencial longitudinal -> gradiente transversal */
/* -> componentes fisicas Ex e Ey.                                            */
/*****************************************************************************/

#include "helm10/field_reconstruction.hpp"

#include "core/fem_scalar.hpp"
#include "meshfree/efgmi_2d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace helm10::field_reconstruction
{
namespace
{

/******************************************************************************/
/* FUNCAO: normalize_vector_field                                             */
/* DESCRICAO: Normaliza um campo vetorial nodal pela maior magnitude presente.*/
/* ENTRADA: ex: std::vector<double> &; ey: std::vector<double> &.             */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void normalize_vector_field(
    std::vector<double> &ex,
    std::vector<double> &ey)
{
    double vmax = 0.0;
    for (size_t i = 0; i < ex.size(); ++i)
    {
        const double v = std::sqrt(ex[i] * ex[i] + ey[i] * ey[i]);
        vmax = std::max(vmax, v);
    }

    if (vmax <= 0.0)
        return;

    for (size_t i = 0; i < ex.size(); ++i)
    {
        ex[i] /= vmax;
        ey[i] /= vmax;
    }
}

/******************************************************************************/
/* FUNCAO: compute_checked_smoothed_gradient                                  */
/* DESCRICAO: Reconstrói o gradiente nodal suavizado do potencial escalar     */
/* usando as Eq. (36) e (37), com verificacao explicita de triangulos         */
/* degenerados antes de dividir por 2A.                                       */
/* ENTRADA: mesh: const Mesh2D &; psi_nodal: const std::vector<double> &;     */
/* dpsi_dx: std::vector<double> &; dpsi_dy: std::vector<double> &;            */
/* area_tolerance: double.                                                    */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void compute_checked_smoothed_gradient(
    const Mesh2D &mesh,
    const std::vector<double> &psi_nodal,
    std::vector<double> &dpsi_dx,
    std::vector<double> &dpsi_dy,
    ElementAssemblyBackend backend,
    double area_tolerance)
{
    if (backend == ElementAssemblyBackend::EfgmiConsistent)
    {
        const auto ctx = efgmi2d::make_context(mesh);
        dpsi_dx.assign(mesh.nodes.size(), 0.0);
        dpsi_dy.assign(mesh.nodes.size(), 0.0);
        for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
        {
            const Node2D &node = mesh.nodes[node_id];
            const auto sample =
                efgmi2d::evaluate_scalar_field(ctx, psi_nodal, node.x, node.y, false);
            dpsi_dx[node_id] = sample.dx;
            dpsi_dy[node_id] = sample.dy;
        }
        return;
    }

    dpsi_dx.assign(mesh.nodes.size(), 0.0);
    dpsi_dy.assign(mesh.nodes.size(), 0.0);
    std::vector<double> area_weight_sum(mesh.nodes.size(), 0.0);

    for (size_t tri_id = 0; tri_id < mesh.tris.size(); ++tri_id)
    {
        const Tri &tri = mesh.tris[tri_id];
        const TriGeom geom = tri_geom(mesh, tri);

        if (geom.A <= area_tolerance)
        {
            std::ostringstream oss;
            oss << "Triangulo degenerado encontrado na reconstrucao de campos: "
                << "tri_id=" << tri_id << " area=" << geom.A;
            throw std::runtime_error(oss.str());
        }

        // Eq. (36) e Eq. (37): o gradiente e constante em cada triangulo P1.
        double tri_dpsi_dx = 0.0;
        double tri_dpsi_dy = 0.0;
        for (int local_i = 0; local_i < 3; ++local_i)
        {
            const double psi_i = psi_nodal[(size_t)tri.v[local_i]];
            tri_dpsi_dx += psi_i * (geom.b[local_i] / (2.0 * geom.A));
            tri_dpsi_dy += psi_i * (geom.c[local_i] / (2.0 * geom.A));
        }

        // A saida final e nodal, entao fazemos uma media ponderada pela area
        // dos triangulos adjacentes. Isso nao altera a fisica do elemento,
        // apenas torna a visualizacao e a exportacao mais legiveis.
        for (int local_i = 0; local_i < 3; ++local_i)
        {
            const int node_id = tri.v[local_i];
            dpsi_dx[(size_t)node_id] += tri_dpsi_dx * geom.A;
            dpsi_dy[(size_t)node_id] += tri_dpsi_dy * geom.A;
            area_weight_sum[(size_t)node_id] += geom.A;
        }
    }

    for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
    {
        if (area_weight_sum[node_id] > 0.0)
        {
            dpsi_dx[node_id] /= area_weight_sum[node_id];
            dpsi_dy[node_id] /= area_weight_sum[node_id];
        }
    }
}

} // namespace

double vacuum_permittivity()
{
    return 8.854187817e-12;
}

double vacuum_permeability()
{
    return 4.0 * std::numbers::pi_v<double> * 1.0e-7;
}

HomogeneousMedium make_homogeneous_medium_from_relative(
    double frequency_hz,
    double eps_r,
    double mu_r)
{
    if (frequency_hz <= 0.0)
        throw std::runtime_error("frequency_hz deve ser > 0");
    if (eps_r <= 0.0)
        throw std::runtime_error("eps_r deve ser > 0");
    if (mu_r <= 0.0)
        throw std::runtime_error("mu_r deve ser > 0");

    HomogeneousMedium medium;
    medium.frequency_hz = frequency_hz;
    medium.omega_rad_s = 2.0 * std::numbers::pi_v<double> * frequency_hz;
    medium.eps_r = eps_r;
    medium.mu_r = mu_r;
    medium.epsilon = eps_r * vacuum_permittivity();
    medium.mu = mu_r * vacuum_permeability();
    medium.k = medium.omega_rad_s * std::sqrt(medium.mu * medium.epsilon);
    return medium;
}

HomogeneousMedium make_homogeneous_medium_from_k(
    double k_medium,
    double eps_r,
    double mu_r)
{
    if (k_medium <= 0.0)
        throw std::runtime_error("k_medium deve ser > 0");
    if (eps_r <= 0.0)
        throw std::runtime_error("eps_r deve ser > 0");
    if (mu_r <= 0.0)
        throw std::runtime_error("mu_r deve ser > 0");

    const double epsilon = eps_r * vacuum_permittivity();
    const double mu = mu_r * vacuum_permeability();
    const double frequency_hz =
        k_medium / (2.0 * std::numbers::pi_v<double> * std::sqrt(mu * epsilon));
    return make_homogeneous_medium_from_relative(frequency_hz, eps_r, mu_r);
}

HomogeneousMedium make_homogeneous_medium_above_kc(
    double kc_max,
    double eps_r,
    double mu_r,
    double safety_factor)
{
    if (eps_r <= 0.0)
        throw std::runtime_error("eps_r deve ser > 0");
    if (mu_r <= 0.0)
        throw std::runtime_error("mu_r deve ser > 0");
    if (safety_factor <= 1.0)
        throw std::runtime_error("safety_factor deve ser > 1");

    const double clamped_kc = std::max(0.0, kc_max);
    const double target_k = std::max(1.0, safety_factor * clamped_kc);
    return make_homogeneous_medium_from_k(target_k, eps_r, mu_r);
}

HomogeneousMedium make_homogeneous_medium_for_tm_reference_ztm(
    double kc_ref,
    double ztm_target,
    double eps_r,
    double mu_r)
{
    if (kc_ref <= 0.0)
        throw std::runtime_error("kc_ref deve ser > 0");
    if (ztm_target <= 0.0)
        throw std::runtime_error("ztm_target deve ser > 0");
    if (eps_r <= 0.0)
        throw std::runtime_error("eps_r deve ser > 0");
    if (mu_r <= 0.0)
        throw std::runtime_error("mu_r deve ser > 0");

    const double epsilon = eps_r * vacuum_permittivity();
    const double mu = mu_r * vacuum_permeability();
    const double eta = std::sqrt(mu / epsilon);
    if (ztm_target >= eta)
    {
        throw std::runtime_error("ztm_target deve ser menor que a impedancia intrinseca do meio");
    }

    // Para TM:
    //   Ztm = eta * sqrt(1 - (kc/k)^2)
    // Fixando Ztm e kc_ref, obtemos o k do meio que realiza essa condicao.
    const double ratio = ztm_target / eta;
    const double denominator = std::sqrt(std::max(1.0e-18, 1.0 - ratio * ratio));
    const double k_medium = kc_ref / denominator;
    return make_homogeneous_medium_from_k(k_medium, eps_r, mu_r);
}

const char *longitudinal_label(LongitudinalScalarKind kind)
{
    switch (kind)
    {
    case LongitudinalScalarKind::TE_Hz:
        return "Hz";
    case LongitudinalScalarKind::TM_Ez:
        return "Ez";
    }
    return "unknown";
}

const char *field_equation_label(LongitudinalScalarKind kind)
{
    switch (kind)
    {
    case LongitudinalScalarKind::TE_Hz:
        return "Eq.(38-39)";
    case LongitudinalScalarKind::TM_Ez:
        return "Eq.(40-41) sem Ztm";
    }
    return "Eq.(?)";
}

ReconstructedField2D reconstruct_transverse_fields(
    const Mesh2D &mesh,
    const std::vector<double> &psi_nodal,
    LongitudinalScalarKind kind,
    double kc,
    const HomogeneousMedium &medium,
    bool normalize_transverse_field,
    double area_tolerance)
{
    return reconstruct_transverse_fields(
        mesh,
        psi_nodal,
        kind,
        kc,
        medium,
        ElementAssemblyBackend::ClosedForm,
        normalize_transverse_field,
        area_tolerance);
}

ReconstructedField2D reconstruct_transverse_fields(
    const Mesh2D &mesh,
    const std::vector<double> &psi_nodal,
    LongitudinalScalarKind kind,
    double kc,
    const HomogeneousMedium &medium,
    ElementAssemblyBackend backend,
    bool normalize_transverse_field,
    double area_tolerance)
{
    if (psi_nodal.size() != mesh.nodes.size())
        throw std::runtime_error("psi_nodal.size() deve coincidir com mesh.nodes.size()");
    if (kc < 0.0)
        throw std::runtime_error("kc deve ser >= 0");
    if (medium.omega_rad_s <= 0.0)
        throw std::runtime_error("omega_rad_s deve ser > 0");
    if (medium.epsilon <= 0.0)
        throw std::runtime_error("epsilon deve ser > 0");
    if (medium.mu <= 0.0)
        throw std::runtime_error("mu deve ser > 0");

    ReconstructedField2D result;
    result.kind = kind;
    result.psi = psi_nodal;
    result.kc = kc;
    result.k = medium.k;
    result.longitudinal_label = longitudinal_label(kind);
    result.field_equations = field_equation_label(kind);

    compute_checked_smoothed_gradient(
        mesh,
        psi_nodal,
        result.dpsi_dx,
        result.dpsi_dy,
        backend,
        area_tolerance);

    result.ex.assign(mesh.nodes.size(), 0.0);
    result.ey.assign(mesh.nodes.size(), 0.0);
    result.ex_without_ztm.assign(mesh.nodes.size(), 0.0);
    result.ey_without_ztm.assign(mesh.nodes.size(), 0.0);
    result.ex_with_ztm.assign(mesh.nodes.size(), 0.0);
    result.ey_with_ztm.assign(mesh.nodes.size(), 0.0);

    const double k_squared = medium.k * medium.k;
    const double kc_squared = kc * kc;
    const double cutoff_tolerance =
        1.0e-12 * std::max({1.0, std::abs(k_squared), std::abs(kc_squared)});
    const double beta_squared = k_squared - kc_squared;

    if (beta_squared <= cutoff_tolerance)
    {
        result.below_cutoff = true;
        result.beta = 0.0;
    }
    else
    {
        result.beta = std::sqrt(beta_squared);
    }

    if (kind == LongitudinalScalarKind::TE_Hz)
    {
        // TE: a variavel escalar representa o potencial longitudinal associado
        // a Hz. Pela convencao escolhida no projeto, os campos transversais
        // sao reconstruidos diretamente pelas Eq. (38) e (39).
        for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
        {
            result.ex[node_id] = -result.dpsi_dy[node_id];
            result.ey[node_id] = result.dpsi_dx[node_id];
            result.ex_without_ztm[node_id] = result.ex[node_id];
            result.ey_without_ztm[node_id] = result.ey[node_id];
        }

        if (result.below_cutoff)
        {
            result.status_label = "te_sign_only_below_cutoff";
            std::ostringstream oss;
            oss << "Modo TE abaixo do corte: k^2-kc^2=" << beta_squared
                << " <= " << cutoff_tolerance
                << ". Os campos exportados preservam apenas o padrao "
                   "transversal das Eq. (38)-(39).";
            result.status_message = oss.str();
        }
        else
        {
            result.status_label = "te_sign_only_above_cutoff";
            std::ostringstream oss;
            oss << "Modo TE acima do corte: beta=" << result.beta
                << ". Os campos exportados seguem as Eq. (38)-(39), sem usar Zte.";
            result.status_message = oss.str();
        }
    }
    else
    {
        // TM: a variavel escalar representa Ez. No modo didatico atual,
        // exportamos Ex e Ey apenas com o gradiente e o sinal das Eq. (40) e
        // (41), sem multiplicar por Ztm. Beta e Ztm continuam sendo
        // calculados como informacao fisica auxiliar quando a frequencia e
        // conhecida.
        for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
        {
            result.ex_without_ztm[node_id] = -result.dpsi_dx[node_id];
            result.ey_without_ztm[node_id] = -result.dpsi_dy[node_id];
            result.ex[node_id] = result.ex_without_ztm[node_id];
            result.ey[node_id] = result.ey_without_ztm[node_id];
        }

        if (!result.below_cutoff)
        {
            result.ztm = result.beta / (medium.omega_rad_s * medium.epsilon);
            result.tm_has_scaled_fields = true;
            for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
            {
                result.ex_with_ztm[node_id] = result.ztm * result.ex_without_ztm[node_id];
                result.ey_with_ztm[node_id] = result.ztm * result.ey_without_ztm[node_id];
            }
        }

        if (result.below_cutoff)
        {
            result.ztm = 0.0;
            result.status_label = "tm_sign_only_below_cutoff";
            std::ostringstream oss;
            oss << "Modo TM abaixo do corte: k^2-kc^2=" << beta_squared
                << " <= " << cutoff_tolerance
                << ". Os campos exportados usam apenas gradiente e sinal; "
                   "Ztm nao foi aplicado.";
            result.status_message = oss.str();
        }
        else
        {
            result.status_label = "tm_sign_only_above_cutoff";
            std::ostringstream oss;
            oss << "Modo TM acima do corte: beta=" << result.beta
                << ", Ztm=" << result.ztm
                << ". Os campos exportados usam apenas gradiente e sinal; "
                   "Ztm nao foi aplicado em Ex e Ey.";
            result.status_message = oss.str();
        }
    }

    if (normalize_transverse_field)
        normalize_vector_field(result.ex, result.ey);

    return result;
}

} // namespace helm10::field_reconstruction
