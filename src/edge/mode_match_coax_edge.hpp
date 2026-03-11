/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge/mode_match_coax_edge.hpp                                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo 2D de elementos de aresta (DOFs, base de Whitney e       */
/* montagem).                                                                 */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.1.        */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/mesh2d.hpp"
#include "core/dense.hpp"
#include "edge/edge_dofs.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <utility>

// ======================================================
//  Utilitarios de Bessel (mesma ideia do caso coaxial escalar)
// ======================================================
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
/* na formulacao analitica do guia coaxial.                                   */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Ym(int m, double x) { return std::cyl_neumann((double)m, x); }

// Jm'(x) = 0.5*(J_{m-1}-J_{m+1}), J0'=-J1
/******************************************************************************/
/* FUNCAO: Jm_prime                                                           */
/* DESCRICAO: Calcula a derivada de J_m(x), usada nas condicoes de contorno   */
/* de modos TE e na equacao transcendental do coaxial.                        */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Jm_prime(int m, double x)
{
    if (m == 0) return -Jm(1, x);
    return 0.5 * (Jm(m - 1, x) - Jm(m + 1, x));
}

// Ym'(x) = 0.5*(Y_{m-1}-Y_{m+1}), Y0'=-Y1
/******************************************************************************/
/* FUNCAO: Ym_prime                                                           */
/* DESCRICAO: Avalia a derivada da funcao de Bessel de segunda especie Y_m(x), usada no problema coaxial. */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Ym_prime(int m, double x)
{
    if (m == 0) return -Ym(1, x);
    return 0.5 * (Ym(m - 1, x) - Ym(m + 1, x));
}

// ======================================================
//  Correlacao de massa (mesmo estilo do modulo escalar)
// ======================================================
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
        for (int j = 0; j < n; j++) row += T(i, j) * y[(size_t)j];
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
    if (nx <= 0 || ny <= 0) return 0.0;
    return std::abs(xy / (nx * ny));
}

/******************************************************************************/
/* FUNCAO: extract_edge_mode_from_Zcol                                        */
/* DESCRICAO: Extrai um modo no espaco de arestas a partir da matriz de       */
/* autovetores column-major.                                                  */
/* ENTRADA: Zcol: const std::vector<double> &; ndof: int; mode_idx: int; x:   */
/* std::vector<double> &.                                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void extract_edge_mode_from_Zcol(
    const std::vector<double> &Zcol, int ndof, int mode_idx,
    std::vector<double> &x)
{
    x.assign((size_t)ndof, 0.0);
    auto idx_col = [&](int i, int j) { return (size_t)j * ndof + i; };
    for (int i = 0; i < ndof; i++) x[(size_t)i] = Zcol[idx_col(i, mode_idx)];
}

// ======================================================
//  Busca de raizes + determinantes (igual ao coaxial escalar)
// ======================================================
/******************************************************************************/
/* FUNCAO: find_root_bisection                                                */
/* DESCRICAO: Encontra raiz em intervalo com troca de sinal via bisseccao, com*/
/* robustez numerica.                                                         */
/* ENTRADA: f: std::function<double(double)>; a: double; b: double; it: int.  */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double find_root_bisection(std::function<double(double)> f, double a, double b, int it = 120)
{
    double fa = f(a), fb = f(b);
    if (fa * fb > 0) throw std::runtime_error("Intervalo sem troca de sinal.");
    for (int k = 0; k < it; k++)
    {
        double c = 0.5 * (a + b);
        double fc = f(c);
        if (fa * fc <= 0) { b = c; fb = fc; }
        else { a = c; fa = fc; }
    }
    return 0.5 * (a + b);
}

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

// Encontra as primeiras pmax raízes positivas para um m
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

    auto F = [&](double k) { return is_TE ? det_TE(m, k, r1, r2) : det_TM(m, k, r1, r2); };

    double k = 1e-6;
    double step = 0.02;
    double fprev = F(k);

    double kmax = (m + 2 * pmax + 20) * M_PI / (r2 - r1);

    while ((int)roots.size() < pmax && k < kmax)
    {
        double k2 = k + step;
        double f2 = F(k2);

        if (std::isfinite(fprev) && std::isfinite(f2) && fprev * f2 < 0.0)
        {
            double r = find_root_bisection(F, k, k2, 160);
            if (roots.empty() || std::abs(r - roots.back()) > 1e-4)
                roots.push_back(r);
        }

        k = k2;
        fprev = f2;

        step = 0.02 + 0.002 * k;
        if (step > 0.2) step = 0.2;
    }

    if ((int)roots.size() < pmax)
        throw std::runtime_error("Nao encontrei raizes coax suficientes (aumente kmax/step).");

    return roots;
}

// ======================================================
//  Potencial coaxial U(r,theta) e suas derivadas
// ======================================================
// U(r,theta) = R(r) * ang(m theta)
// R(r) = Jm(k r) + B Ym(k r)
// B escolhido para satisfazer a BC em r1:
//   TM:  R(r1)=0      => B = -Jm(k r1)/Ym(k r1)
//   TE:  R'(r1)=0     => B = -Jm'(k r1)/Ym'(k r1)  (derivada no argumento)
/******************************************************************************/
/* FUNCAO: coax_B                                                             */
/* DESCRICAO: Calcula o fator B da combinacao radial J_m + B Y_m no guia coaxial.     */
/* ENTRADA: m: int; k: double; r1: double; is_TE: bool.                       */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double coax_B(int m, double k, double r1, bool is_TE)
{
    if (!is_TE)
    {
        double denom = Ym(m, k * r1);
        return -Jm(m, k * r1) / denom;
    }
    else
    {
        double denom = Ym_prime(m, k * r1);
        return -Jm_prime(m, k * r1) / denom;
    }
}

/******************************************************************************/
/* FUNCAO: coax_ang                                                           */
/* DESCRICAO: Avalia o fator angular cos(m phi) ou sin(m phi) associado ao indice p do modo coaxial. */
/* ENTRADA: m: int; th: double; use_sin: bool; ang: double &; dang_dth: double*/
/* &.                                                                         */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void coax_ang(int m, double th, bool use_sin, double &ang, double &dang_dth)
{
    if (m == 0)
    {
        ang = 1.0;
        dang_dth = 0.0;
        return;
    }
    if (!use_sin)
    {
        ang = std::cos(m * th);
        dang_dth = -m * std::sin(m * th);
    }
    else
    {
        ang = std::sin(m * th);
        dang_dth =  m * std::cos(m * th);
    }
}

// Calcula dU/dr e (1/r) dU/dtheta em (r,theta)
// R'(r) = d/dr [Jm(k r) + B Ym(k r)] = k*(Jm'(k r) + B Ym'(k r))
/******************************************************************************/
/* FUNCAO: coax_polar_derivatives_U                                           */
/* DESCRICAO: Calcula derivadas de U em coordenadas polares para modos analiticos do coaxial. */
/* ENTRADA: m: int; k: double; B: double; r: double; th: double; use_sin:     */
/* bool; dUdr: double &; dUdth_over_r: double &.                              */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void coax_polar_derivatives_U(
    int m, double k, double B,
    double r, double th, bool use_sin,
    double &dUdr, double &dUdth_over_r)
{
    double ang, dang;
    coax_ang(m, th, use_sin, ang, dang);

    double kr = k * r;
    double R = Jm(m, kr) + B * Ym(m, kr);
    double Rp = k * (Jm_prime(m, kr) + B * Ym_prime(m, kr));

    dUdr = Rp * ang;

    if (r > 1e-14)
        dUdth_over_r = (R * dang) / r;
    else
        dUdth_over_r = 0.0;
}

/******************************************************************************/
/* FUNCAO: gradU_cart                                                         */
/* DESCRICAO: Converte gradiente de U de coordenadas polares para cartesianas (x,y).  */
/* ENTRADA: dUdr: double; dUdth_over_r: double; th: double; dUdx: double &;   */
/* dUdy: double &.                                                            */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void gradU_cart(double dUdr, double dUdth_over_r, double th, double &dUdx, double &dUdy)
{
    double c = std::cos(th), s = std::sin(th);
    dUdx = c * dUdr - s * dUdth_over_r;
    dUdy = s * dUdr + c * dUdth_over_r;
}

// Ft = z x grad U = (-dU/dy, dU/dx)
/******************************************************************************/
/* FUNCAO: Ft_from_gradU                                                      */
/* DESCRICAO: Monta o campo transversal F_t a partir do gradiente de U no caso coaxial. */
/* ENTRADA: dUdx: double; dUdy: double; Fx: double &; Fy: double &.           */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void Ft_from_gradU(double dUdx, double dUdy, double &Fx, double &Fy)
{
    Fx = -dUdy;
    Fy =  dUdx;
}

// ======================================================
//  Projecao em aresta com Gauss de 2 pontos (template C++17)
// ======================================================
template <class FieldXY>
/******************************************************************************/
/* FUNCAO: edge_dof_from_field_gauss2                                         */
/* DESCRICAO: Integra o campo tangencial ao longo da aresta via Gauss-2 para obter grau de liberdade de aresta. */
/* ENTRADA: x0: double; y0: double; x1: double; y1: double; field_xy: FieldXY */
/* &&.                                                                        */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double edge_dof_from_field_gauss2(
    double x0, double y0, double x1, double y1,
    FieldXY &&field_xy)
{
    double dx = x1 - x0, dy = y1 - y0;
    double L = std::sqrt(dx * dx + dy * dy);
    if (L <= 0) return 0.0;

    double tx = dx / L, ty = dy / L;

    double s1 = 0.5 * L * (1.0 - 1.0 / std::sqrt(3.0));
    double s2 = 0.5 * L * (1.0 + 1.0 / std::sqrt(3.0));
    double w  = 0.5 * L;

    auto sample = [&](double s)
    {
        double x = x0 + tx * s;
        double y = y0 + ty * s;
        auto F = field_xy(x, y);
        return F.first * tx + F.second * ty;
    };

    double I = w * sample(s1) + w * sample(s2);
    return I / L; // média ao longo da aresta
}

// Constroi vetor analitico de DOFs de aresta y para um (m,p), usando k=raiz
// y[dof] = (1/L) ∫ Ft·t ds
/******************************************************************************/
/* FUNCAO: analytic_coax_edges_dofs                                           */
/* DESCRICAO: Calcula expressao analitica de referencia para comparacao e     */
/* validacao dos resultados numericos.                                        */
/* ENTRADA: mesh: const Mesh2D &; ed: const EdgeDofs &; r1: double; r2:       */
/* double; m: int; p: int; roots_m: const std::vector<double> &; is_TE: bool; */
/* use_sin: bool; y: std::vector<double> &.                                   */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void analytic_coax_edges_dofs(
    const Mesh2D &mesh,
    const EdgeDofs &ed,
    double r1, double r2,
    int m, int p,
    const std::vector<double> &roots_m,
    bool is_TE,
    bool use_sin,
    std::vector<double> &y)
{
    (void)r2; // r2 entra via roots_m (determinante)

    double k = roots_m[(size_t)(p - 1)];
    double B = coax_B(m, k, r1, is_TE);

    y.assign((size_t)ed.ndof, 0.0);

    for (int eid = 0; eid < (int)ed.edges.size(); eid++)
    {
        int dof = ed.edge_to_dof[eid];
        if (dof < 0) continue; // eliminado (TE/Et=0)

        const auto &E = ed.edges[eid];
        const auto &n0 = mesh.nodes[E.n0];
        const auto &n1 = mesh.nodes[E.n1];

        auto field_xy = [&](double x, double y_)
        {
            double r = std::sqrt(x * x + y_ * y_);
            double th = std::atan2(y_, x);

            // dominio coaxial: r deveria ser >= r1, mas a malha numerica pode ter pequeno desvio
            // se r<r1 ainda funciona, mas Ym perto de 0 e perigoso. Em coax, r1>0.
            double dUdr, dUdth_over_r;
            coax_polar_derivatives_U(m, k, B, r, th, use_sin, dUdr, dUdth_over_r);

            double dUdx, dUdy;
            gradU_cart(dUdr, dUdth_over_r, th, dUdx, dUdy);

            double Fx, Fy;
            Ft_from_gradU(dUdx, dUdy, Fx, Fy);

            return std::pair<double,double>{Fx, Fy};
        };

        y[(size_t)dof] = edge_dof_from_field_gauss2(n0.x, n0.y, n1.x, n1.y, field_xy);
    }
}

// ======================================================
//  API de casamento modal (TE e TM) - mesmo estilo do caso circular
// ======================================================
struct CoaxEdgeModeID
{
    int m = 0;
    int p = 0;
    double kc_ana = 0.0; // = k
    double rho = 0.0;
    bool use_sin = false;
};

/******************************************************************************/
/* FUNCAO: match_coax_edge_mode_by_mass_correlation_TE                        */
/* DESCRICAO: Classifica modo de aresta coaxial TE por correlacao de massa.   */
/* ENTRADA: mesh: const Mesh2D &; r1: double; r2: double; T: const DenseMat &;*/
/* ed: const EdgeDofs &; Zcol: const std::vector<double> &; mode_idx: int;    */
/* mmax: int; pmax: int.                                                      */
/* SAIDA: CoaxEdgeModeID.                                                     */
/******************************************************************************/
inline CoaxEdgeModeID match_coax_edge_mode_by_mass_correlation_TE(
    const Mesh2D &mesh,
    double r1, double r2,
    const DenseMat &T,
    const EdgeDofs &ed,
    const std::vector<double> &Zcol,
    int mode_idx,
    int mmax = 6,
    int pmax = 6)
{
    std::vector<double> x;
    extract_edge_mode_from_Zcol(Zcol, ed.ndof, mode_idx, x);

    CoaxEdgeModeID best;
    best.rho = -1.0;

    for (int m = 0; m <= mmax; m++)
    {
        bool is_TE = true;
        auto roots = coax_roots(m, pmax, r1, r2, /*is_TE*/ true);

        for (int p = 1; p <= pmax; p++)
        {
            std::vector<double> y_cos, y_sin;

            analytic_coax_edges_dofs(mesh, ed, r1, r2, m, p, roots, is_TE, false, y_cos);
            double rho_cos = mass_correlation_abs(T, x, y_cos);

            double rho = rho_cos;
            bool use_sin = false;

            if (m > 0)
            {
                analytic_coax_edges_dofs(mesh, ed, r1, r2, m, p, roots, is_TE, true, y_sin);
                double rho_sin = mass_correlation_abs(T, x, y_sin);
                if (rho_sin > rho_cos) { rho = rho_sin; use_sin = true; }
            }

            if (rho > best.rho)
            {
                double k = roots[(size_t)(p - 1)];
                best = {m, p, k, rho, use_sin};
            }
        }
    }

    return best;
}

/******************************************************************************/
/* FUNCAO: match_coax_edge_mode_by_mass_correlation_TM                        */
/* DESCRICAO: Classifica modo de aresta coaxial TM por correlacao de massa.   */
/* ENTRADA: mesh: const Mesh2D &; r1: double; r2: double; T: const DenseMat &;*/
/* ed: const EdgeDofs &; Zcol: const std::vector<double> &; mode_idx: int;    */
/* mmax: int; pmax: int.                                                      */
/* SAIDA: CoaxEdgeModeID.                                                     */
/******************************************************************************/
inline CoaxEdgeModeID match_coax_edge_mode_by_mass_correlation_TM(
    const Mesh2D &mesh,
    double r1, double r2,
    const DenseMat &T,
    const EdgeDofs &ed,
    const std::vector<double> &Zcol,
    int mode_idx,
    int mmax = 6,
    int pmax = 6)
{
    std::vector<double> x;
    extract_edge_mode_from_Zcol(Zcol, ed.ndof, mode_idx, x);

    CoaxEdgeModeID best;
    best.rho = -1.0;

    for (int m = 0; m <= mmax; m++)
    {
        bool is_TE = false;
        auto roots = coax_roots(m, pmax, r1, r2, /*is_TE*/ false);

        for (int p = 1; p <= pmax; p++)
        {
            std::vector<double> y_cos, y_sin;

            analytic_coax_edges_dofs(mesh, ed, r1, r2, m, p, roots, is_TE, false, y_cos);
            double rho_cos = mass_correlation_abs(T, x, y_cos);

            double rho = rho_cos;
            bool use_sin = false;

            if (m > 0)
            {
                analytic_coax_edges_dofs(mesh, ed, r1, r2, m, p, roots, is_TE, true, y_sin);
                double rho_sin = mass_correlation_abs(T, x, y_sin);
                if (rho_sin > rho_cos) { rho = rho_sin; use_sin = true; }
            }

            if (rho > best.rho)
            {
                double k = roots[(size_t)(p - 1)];
                best = {m, p, k, rho, use_sin};
            }
        }
    }

    return best;
}
