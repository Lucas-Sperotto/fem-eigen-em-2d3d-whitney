/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mode_match_circle.hpp                                    */
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
#include "mesh2d.hpp"
#include "helm10_scalar_system.hpp"
#include "dense.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>
#include <functional>
#include <stdexcept>

// ---------- utilitarios de Bessel ----------
/******************************************************************************/
/* FUNCAO: Jm                                                                 */
/* DESCRICAO: Avalia a funcao de Bessel de primeira especie J_m(x), usada nas solucoes analiticas em coordenadas cilindricas. */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Jm(int m, double x)
{
    // C++17/20: std::cyl_bessel_j
    return std::cyl_bessel_j((double)m, x);
}

// Jm'(x) = 0.5*(J_{m-1}(x) - J_{m+1}(x))
/******************************************************************************/
/* FUNCAO: Jm_prime                                                           */
/* DESCRICAO: Avalia a derivada de J_m(x), necessaria nas equacoes de contorno de modos TE/TM cilindricos. */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Jm_prime(int m, double x)
{
    if (m == 0)
    {
        // J0' = -J1
        return -Jm(1, x);
    }
    return 0.5 * (Jm(m - 1, x) - Jm(m + 1, x));
}

/******************************************************************************/
/* FUNCAO: find_root_bisection                                                */
/* DESCRICAO: Encontra raiz em intervalo com troca de sinal via bisseccao, com*/
/* robustez numerica.                                                         */
/* ENTRADA: f: std::function<double(double)>; a: double; b: double; it: int.  */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double find_root_bisection(std::function<double(double)> f, double a, double b, int it = 80)
{
    double fa = f(a), fb = f(b);
    if (fa * fb > 0)
        throw std::runtime_error("Intervalo sem troca de sinal para raiz.");
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

// Encontra as primeiras nroots raízes positivas de f(x)=0 via varredura+ bissecção
/******************************************************************************/
/* FUNCAO: bessel_roots                                                       */
/* DESCRICAO: Calcula zeros de J_m ou J_m' para obter valores analiticos de corte (k_c) em guias circulares. */
/* ENTRADA: m: int; nroots: int; derivative: bool.                            */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> bessel_roots(int m, int nroots, bool derivative)
{
    std::vector<double> roots;
    roots.reserve(nroots);

    auto f = [&](double x)
    {
        return derivative ? Jm_prime(m, x) : Jm(m, x);
    };

    // varredura por intervalos; passo adaptável
    double x = 1e-6;
    double step = 0.05; // suficiente para começar
    double fprev = f(x);

    // limite superior seguro: cresce com nroots e m
    double xmax = (m + 2 * nroots + 10) * M_PI;

    while ((int)roots.size() < nroots && x < xmax)
    {
        double x2 = x + step;
        double f2 = f(x2);

        if (fprev == 0.0)
        {
            roots.push_back(x);
        }
        else if (fprev * f2 < 0.0)
        {
            // bracket achado
            double r = find_root_bisection(f, x, x2, 90);
            // evita duplicata
            if (roots.empty() || std::abs(r - roots.back()) > 1e-3)
                roots.push_back(r);
        }

        x = x2;
        fprev = f2;

        // passo um pouco maior conforme x cresce
        step = 0.05 + 0.002 * x;
        if (step > 0.5)
            step = 0.5;
    }

    if ((int)roots.size() < nroots)
    {
        throw std::runtime_error("Nao encontrei raizes suficientes de Bessel.");
    }
    return roots;
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
/* DESCRICAO: Calcula a correlacao modal absoluta ponderada pela massa entre  */
/* dois vetores de modo.                                                      */
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
inline void extract_mode_from_Zcol(
    const std::vector<double> &Zcol, int ndof, int mode_idx,
    std::vector<double> &x)
{
    x.assign((size_t)ndof, 0.0);
    auto idx_col = [&](int i, int j)
    { return (size_t)j * ndof + i; };
    for (int i = 0; i < ndof; i++)
    {
        x[(size_t)i] = Zcol[idx_col(i, mode_idx)];
    }
}

/******************************************************************************/
/* FUNCAO: restrict_nodes_to_dofs                                             */
/* DESCRICAO: Realiza transformacao de representacao de campos/modos para     */
/* analise e pos-processamento.                                               */
/* ENTRADA: phi_nodes: const std::vector<double> &; dof_map: const            */
/* std::vector<int> &; ndof: int; phi_dof: std::vector<double> &.             */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void restrict_nodes_to_dofs(
    const std::vector<double> &phi_nodes,
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

struct CircModeID
{
    int m = 0;
    int p = 0;         // índice radial (1..)
    double root = 0.0; // x_{m,p} ou x'_{m,p}
    double kc_ana = 0.0;
    double rho = 0.0;
};

/******************************************************************************/
/* FUNCAO: analytic_phi_circle_on_nodes                                       */
/* DESCRICAO: Avalia o potencial modal analitico circular em todos os nos da  */
/* malha.                                                                     */
/* ENTRADA: mesh: const Mesh2D &; R: double; bc: ScalarBC; m: int; p: int;    */
/* roots_m: const std::vector<double> &; use_sin: bool; phi_nodes: // <-- NOVO*/
/* std::vector<double> &.                                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void analytic_phi_circle_on_nodes(
    const Mesh2D &mesh,
    double R,
    ScalarBC bc,
    int m, int p,
    const std::vector<double> &roots_m,
    bool use_sin, // <-- NOVO
    std::vector<double> &phi_nodes)
{
    double root = roots_m[(size_t)(p - 1)];
    double kc = root / R;

    phi_nodes.assign(mesh.nodes.size(), 0.0);
    for (size_t i = 0; i < mesh.nodes.size(); i++)
    {
        double x = mesh.nodes[i].x;
        double y = mesh.nodes[i].y;
        double r = std::sqrt(x * x + y * y);
        double th = std::atan2(y, x);

        double ang = use_sin ? std::sin(m * th) : std::cos(m * th);
        phi_nodes[i] = Jm(m, kc * r) * ang;
    }
}

/******************************************************************************/
/* FUNCAO: match_circle_mode_by_mass_correlation                              */
/* DESCRICAO: Identifica o modo analitico circular mais correlacionado ao modo*/
/* FEM.                                                                       */
/* ENTRADA: mesh: const Mesh2D &; R: double; sys: const ScalarSystem &; Zcol: */
/* const std::vector<double> &; mode_idx: int; bc: ScalarBC; mmax: int; pmax: */
/* int.                                                                       */
/* SAIDA: CircModeID.                                                         */
/******************************************************************************/
inline CircModeID match_circle_mode_by_mass_correlation(
    const Mesh2D &mesh,
    double R,
    const ScalarSystem &sys,
    const std::vector<double> &Zcol,
    int mode_idx,
    ScalarBC bc,
    int mmax = 6,
    int pmax = 6)
{
    std::vector<double> x;
    extract_mode_from_Zcol(Zcol, sys.ndof, mode_idx, x);

    CircModeID best;
    best.rho = -1.0;

    std::vector<double> phi_nodes, phi_dof;

    for (int m = 0; m <= mmax; m++)
    {
        bool deriv = (bc == ScalarBC::TE_Neumann);
        auto roots = bessel_roots(m, pmax, deriv);

        for (int p = 1; p <= pmax; p++)
        {
            // ramo cosseno
            analytic_phi_circle_on_nodes(mesh, R, bc, m, p, roots, false, phi_nodes);
            restrict_nodes_to_dofs(phi_nodes, sys.dof_map, sys.ndof, phi_dof);
            double rho_cos = mass_correlation_abs(sys.T, x, phi_dof);

            double rho = rho_cos;

            // ramo seno (apenas para m>0)
            if (m > 0)
            {
                analytic_phi_circle_on_nodes(mesh, R, bc, m, p, roots, true, phi_nodes);
                restrict_nodes_to_dofs(phi_nodes, sys.dof_map, sys.ndof, phi_dof);
                double rho_sin = mass_correlation_abs(sys.T, x, phi_dof);
                rho = std::max(rho_cos, rho_sin);
            }

            if (rho > best.rho)
            {
                double root = roots[(size_t)(p - 1)];
                best = {m, p, root, root / R, rho};
            }
        }
    }

    return best;
}
