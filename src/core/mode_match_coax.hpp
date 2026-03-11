/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mode_match_coax.hpp                                      */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo escalar 2D (malha, montagem e identificacao modal).      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.1.          */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/mesh2d.hpp"
#include "core/helm10_scalar_system.hpp"
#include "core/dense.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <limits>

// ---------- utilitarios de Bessel ----------
/******************************************************************************/
/* FUNCAO: Jm                                                                 */
/* DESCRICAO: Avalia a funcao de Bessel de primeira especie J_m(x), usada nas solucoes analiticas em coordenadas cilindricas. */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Jm(int m, double x) { return std::cyl_bessel_j((double)m, x); }
/******************************************************************************/
/* FUNCAO: Ym                                                                 */
/* DESCRICAO: Avalia a funcao de Bessel de segunda especie Y_m(x), necessaria */
/* na combinacao radial da solucao analitica coaxial.                         */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Ym(int m, double x) { return std::cyl_neumann((double)m, x); }

// Jm'(x) = 0.5*(J_{m-1}(x) - J_{m+1}(x)),   Y idem
/******************************************************************************/
/* FUNCAO: Jm_prime                                                           */
/* DESCRICAO: Calcula derivada de J_m(x), usada nas condicoes de contorno     */
/* Neumann (modos TE) do guia coaxial.                                        */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Jm_prime(int m, double x)
{
    if (m == 0)
        return -Jm(1, x); // J0' = -J1
    return 0.5 * (Jm(m - 1, x) - Jm(m + 1, x));
}
/******************************************************************************/
/* FUNCAO: Ym_prime                                                           */
/* DESCRICAO: Avalia a derivada da funcao de Bessel de segunda especie Y_m(x), usada no problema coaxial. */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Ym_prime(int m, double x)
{
    if (m == 0)
        return -Ym(1, x); // Y0' = -Y1
    return 0.5 * (Ym(m - 1, x) - Ym(m + 1, x));
}

// ---------- produtos internos com massa ----------
/******************************************************************************/
/* FUNCAO: mass_inner                                                         */
/* DESCRICAO: Calcula produto interno ponderado pela matriz de massa, base    */
/* para correlacao e normalizacao modal.                                      */
/* ENTRADA: T: const DenseMat &; x: const std::vector<double> &; y: const     */
/* std::vector<double> &.                                                     */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double mass_inner(const DenseMat &T, const std::vector<double> &x, const std::vector<double> &y)
{
    const int n = T.n;
    double s = 0.0;
    for (int i = 0; i < n; i++)
    {
        double row = 0.0;
        const double xi = x[(size_t)i];
        for (int j = 0; j < n; j++)
        {
            row += T(i, j) * y[(size_t)j];
        }
        s += xi * row;
    }
    return s;
}
/******************************************************************************/
/* FUNCAO: mass_norm                                                          */
/* DESCRICAO: Calcula a norma energetica induzida pela matriz de massa.       */
/* ENTRADA: T: const DenseMat &; x: const std::vector<double> &.              */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double mass_norm(const DenseMat &T, const std::vector<double> &x)
{
    double v = mass_inner(T, x, x);
    return (v > 0) ? std::sqrt(v) : 0.0;
}
/******************************************************************************/
/* FUNCAO: mass_correlation_abs                                               */
/* DESCRICAO: Calcula correlacao modal absoluta entre autovetor FEM e modo    */
/* analitico usando produto interno de massa.                                 */
/* ENTRADA: T: const DenseMat &; x: const std::vector<double> &; y: const     */
/* std::vector<double> &.                                                     */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double mass_correlation_abs(const DenseMat &T, const std::vector<double> &x, const std::vector<double> &y)
{
    double xy = mass_inner(T, x, y);
    double nx = mass_norm(T, x);
    double ny = mass_norm(T, y);
    if (nx <= 0 || ny <= 0)
        return 0.0;
    return std::abs(xy / (nx * ny));
}

/******************************************************************************/
/* FUNCAO: extract_mode_from_Zcol                                             */
/* DESCRICAO: Extrai um autovetor (coluna) do formato column-major para vetor */
/* de trabalho.                                                               */
/* ENTRADA: Zcol: const std::vector<double> &; ndof: int; mode_idx: int; x:   */
/* std::vector<double> &.                                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void extract_mode_from_Zcol(const std::vector<double> &Zcol, int ndof, int mode_idx, std::vector<double> &x)
{
    x.assign((size_t)ndof, 0.0);
    auto idx_col = [&](int i, int j)
    { return (size_t)j * ndof + i; };
    for (int i = 0; i < ndof; i++)
        x[(size_t)i] = Zcol[idx_col(i, mode_idx)];
}

/******************************************************************************/
/* FUNCAO: restrict_nodes_to_dofs                                             */
/* DESCRICAO: Realiza transformacao de representacao de campos/modos para     */
/* analise e pos-processamento.                                               */
/* ENTRADA: phi_nodes: const std::vector<double> &; dof_map: const            */
/* std::vector<int> &; ndof: int; phi_dof: std::vector<double> &.             */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void restrict_nodes_to_dofs(const std::vector<double> &phi_nodes,
                                   const std::vector<int> &dof_map,
                                   int ndof,
                                   std::vector<double> &phi_dof)
{
    phi_dof.assign((size_t)ndof, 0.0);
    for (size_t i = 0; i < dof_map.size(); i++)
    {
        int d = dof_map[i];
        if (d >= 0)
            phi_dof[d] = phi_nodes[i];
    }
}

// ---------- busca de raizes ----------
/******************************************************************************/
/* FUNCAO: find_root_bisection                                                */
/* DESCRICAO: Encontra raiz em intervalo com troca de sinal via bisseccao, com*/
/* robustez numerica.                                                         */
/* ENTRADA: f: std::function<double(double)>; a: double; b: double; it: int.  */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double find_root_bisection(std::function<double(double)> f, double a, double b, int it = 100)
{
    double fa = f(a), fb = f(b);
    if (fa * fb > 0)
        throw std::runtime_error("Intervalo sem troca de sinal.");
    for (int k = 0; k < it; k++)
    {
        double c = 0.5 * (a + b);
        double fc = f(c);
        if (fa * fc <= 0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }
    return 0.5 * (a + b);
}

// ---------- determinantes do problema coaxial ----------
// TM (Dirichlet-Dirichlet): Jm(k r1) Ym(k r2) - Jm(k r2) Ym(k r1) = 0
/******************************************************************************/
/* FUNCAO: det_TM                                                             */
/* DESCRICAO: Monta o determinante da equacao caracteristica de modos TM em guia coaxial para busca de raizes. */
/* ENTRADA: m: int; k: double; r1: double; r2: double.                        */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double det_TM(int m, double k, double r1, double r2)
{
    double a = Jm(m, k * r1) * Ym(m, k * r2);
    double b = Jm(m, k * r2) * Ym(m, k * r1);
    return a - b;
}

// TE (Neumann-Neumann): Jm'(k r1) Ym'(k r2) - Jm'(k r2) Ym'(k r1) = 0
/******************************************************************************/
/* FUNCAO: det_TE                                                             */
/* DESCRICAO: Monta o determinante da equacao caracteristica de modos TE em guia coaxial para busca de raizes. */
/* ENTRADA: m: int; k: double; r1: double; r2: double.                        */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double det_TE(int m, double k, double r1, double r2)
{
    double a = Jm_prime(m, k * r1) * Ym_prime(m, k * r2);
    double b = Jm_prime(m, k * r2) * Ym_prime(m, k * r1);
    return a - b;
}

// Encontra as primeiras pmax raízes positivas do determinante para um m
/******************************************************************************/
/* FUNCAO: coax_roots                                                         */
/* DESCRICAO: Varre o intervalo espectral e localiza raizes da equacao caracteristica do guia coaxial. */
/* ENTRADA: m: int; pmax: int; r1: double; r2: double; is_TE: bool.           */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> coax_roots(int m, int pmax, double r1, double r2, bool is_TE)
{
    std::vector<double> roots;
    roots.reserve(pmax);

    auto F = [&](double k)
    {
        return is_TE ? det_TE(m, k, r1, r2) : det_TM(m, k, r1, r2);
    };

    double k = 1e-6;
    double step = 0.02; // coax costuma ter raízes pequenas; começa fino
    double fprev = F(k);

    double kmax = (m + 2 * pmax + 20) * M_PI / (r2 - r1); // limite razoável

    while ((int)roots.size() < pmax && k < kmax)
    {
        double k2 = k + step;
        double f2 = F(k2);

        if (std::isfinite(fprev) && std::isfinite(f2) && fprev * f2 < 0.0)
        {
            double r = find_root_bisection(F, k, k2, 120);
            if (roots.empty() || std::abs(r - roots.back()) > 1e-4)
                roots.push_back(r);
        }

        k = k2;
        fprev = f2;

        step = 0.02 + 0.002 * k;
        if (step > 0.2)
            step = 0.2;
    }

    if ((int)roots.size() < pmax)
        throw std::runtime_error("Nao encontrei raizes coax suficientes (aumente kmax/step).");

    return roots;
}

// ---------- forma modal analitica nos nos ----------
// phi(r,theta) = [Jm(k r) + B Ym(k r)] * cos(m theta)  (ou sin)
// B escolhido para satisfazer BC em r1:
// TM:  Jm(k r1)+B Ym(k r1)=0 => B = -Jm/ Ym
// TE:  Jm'(k r1)+B Ym'(k r1)=0 => B = -Jm'/Ym'
/******************************************************************************/
/* FUNCAO: analytic_phi_coax_on_nodes                                         */
/* DESCRICAO: Avalia o potencial modal analitico coaxial em todos os nos da   */
/* malha.                                                                     */
/* ENTRADA: mesh: const Mesh2D &; r1: double; r2: double; bc: ScalarBC; m:    */
/* int; k: double; use_sin: bool; phi_nodes: std::vector<double> &.           */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void analytic_phi_coax_on_nodes(
    const Mesh2D &mesh,
    double r1, double r2,
    ScalarBC bc,
    int m, double k,
    bool use_sin,
    std::vector<double> &phi_nodes)
{
    (void)r2; // usado só no determinante/raízes
    bool is_TE = (bc == ScalarBC::TE_Neumann);

    double B = 0.0;
    if (!is_TE)
    {
        double denom = Ym(m, k * r1);
        B = -Jm(m, k * r1) / denom;
    }
    else
    {
        double denom = Ym_prime(m, k * r1);
        B = -Jm_prime(m, k * r1) / denom;
    }

    phi_nodes.assign(mesh.nodes.size(), 0.0);
    for (size_t i = 0; i < mesh.nodes.size(); i++)
    {
        double x = mesh.nodes[i].x;
        double y = mesh.nodes[i].y;
        double r = std::sqrt(x * x + y * y);
        double th = std::atan2(y, x);

        double ang = (m == 0) ? 1.0 : (use_sin ? std::sin(m * th) : std::cos(m * th));
        double radial = Jm(m, k * r) + B * Ym(m, k * r);
        phi_nodes[i] = radial * ang;
    }
}

struct CoaxModeID
{
    int m = 0;
    int p = 0; // índice radial (1..)
    double kc_ana = 0.0;
    double rho = 0.0;
};

/******************************************************************************/
/* FUNCAO: match_coax_mode_by_mass_correlation                                */
/* DESCRICAO: Identifica o modo analitico coaxial mais correlacionado ao modo */
/* FEM.                                                                       */
/* ENTRADA: mesh: const Mesh2D &; r1: double; r2: double; sys: const          */
/* ScalarSystem &; Zcol: const std::vector<double> &; mode_idx: int; bc:      */
/* ScalarBC; mmax: int; pmax: int.                                            */
/* SAIDA: CoaxModeID.                                                         */
/******************************************************************************/
inline CoaxModeID match_coax_mode_by_mass_correlation(
    const Mesh2D &mesh,
    double r1, double r2,
    const ScalarSystem &sys,
    const std::vector<double> &Zcol,
    int mode_idx,
    ScalarBC bc,
    int mmax = 6,
    int pmax = 6)
{
    std::vector<double> x;
    extract_mode_from_Zcol(Zcol, sys.ndof, mode_idx, x);

    CoaxModeID best;
    best.rho = -1.0;

    std::vector<double> phi_nodes, phi_dof;

    bool is_TE = (bc == ScalarBC::TE_Neumann);

    for (int m = 0; m <= mmax; m++)
    {
        auto roots = coax_roots(m, pmax, r1, r2, is_TE);

        for (int p = 1; p <= pmax; p++)
        {
            double k = roots[(size_t)(p - 1)];

            // ramo cosseno
            analytic_phi_coax_on_nodes(mesh, r1, r2, bc, m, k, false, phi_nodes);
            restrict_nodes_to_dofs(phi_nodes, sys.dof_map, sys.ndof, phi_dof);
            double rho_cos = mass_correlation_abs(sys.T, x, phi_dof);

            double rho = rho_cos;

            // ramo seno (para m>0)
            if (m > 0)
            {
                analytic_phi_coax_on_nodes(mesh, r1, r2, bc, m, k, true, phi_nodes);
                restrict_nodes_to_dofs(phi_nodes, sys.dof_map, sys.ndof, phi_dof);
                double rho_sin = mass_correlation_abs(sys.T, x, phi_dof);
                rho = std::max(rho_cos, rho_sin);
            }

            if (rho > best.rho)
            {
                best = {m, p, k, rho};
            }
        }
    }

    return best;
}
