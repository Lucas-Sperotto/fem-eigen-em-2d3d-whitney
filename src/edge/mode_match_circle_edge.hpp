/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge/mode_match_circle_edge.hpp                               */
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

// ---------------- utilitarios de Bessel (igual ao escalar) ----------------
/******************************************************************************/
/* FUNCAO: Jm                                                                 */
/* DESCRICAO: Avalia a funcao de Bessel de primeira especie J_m(x), usada nas solucoes analiticas em coordenadas cilindricas. */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Jm(int m, double x) { return std::cyl_bessel_j((double)m, x); }

// Jm'(x) = 0.5*(J_{m-1}(x) - J_{m+1}(x))
/******************************************************************************/
/* FUNCAO: Jm_prime                                                           */
/* DESCRICAO: Calcula derivada de J_m(x), necessaria para modos TE e para a   */
/* busca de raizes de J_m' no guia circular.                                  */
/* ENTRADA: m: int; x: double.                                                */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double Jm_prime(int m, double x)
{
    if (m == 0) return -Jm(1, x); // J0' = -J1
    return 0.5 * (Jm(m - 1, x) - Jm(m + 1, x));
}

/******************************************************************************/
/* FUNCAO: find_root_bisection                                                */
/* DESCRICAO: Encontra raiz em intervalo com troca de sinal via bisseccao, com*/
/* robustez numerica.                                                         */
/* ENTRADA: f: std::function<double(double)>; a: double; b: double; it: int.  */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double find_root_bisection(std::function<double(double)> f, double a, double b, int it = 90)
{
    double fa = f(a), fb = f(b);
    if (fa * fb > 0) throw std::runtime_error("Intervalo sem troca de sinal para raiz.");
    for (int k = 0; k < it; k++)
    {
        double c = 0.5 * (a + b);
        double fc = f(c);
        if (fa * fc <= 0) { b = c; fb = fc; }
        else { a = c; fa = fc; }
    }
    return 0.5 * (a + b);
}

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

    auto f = [&](double x) { return derivative ? Jm_prime(m, x) : Jm(m, x); };

    double x = 1e-6;
    double step = 0.05;
    double fprev = f(x);

    double xmax = (m + 2 * nroots + 10) * M_PI;

    while ((int)roots.size() < nroots && x < xmax)
    {
        double x2 = x + step;
        double f2 = f(x2);

        if (fprev == 0.0) roots.push_back(x);
        else if (fprev * f2 < 0.0)
        {
            double r = find_root_bisection(f, x, x2, 90);
            if (roots.empty() || std::abs(r - roots.back()) > 1e-3)
                roots.push_back(r);
        }

        x = x2;
        fprev = f2;

        step = 0.05 + 0.002 * x;
        if (step > 0.5) step = 0.5;
    }

    if ((int)roots.size() < nroots)
        throw std::runtime_error("Nao encontrei raizes suficientes de Bessel.");

    return roots;
}

// ---------------- correlacao de massa nos DOFs de aresta ----------------
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

// ---------------- utilitarios geometricos ----------------
struct CircEdgeModeID
{
    int m = 0;
    int p = 0;         // radial index (1..)
    double root = 0.0; // chi_{m,p} or chi'_{m,p}
    double kc_ana = 0.0;
    double rho = 0.0;
    bool use_sin = false; // qual base angular forneceu melhor correspondencia
};

// dU/dr e (1/r)dU/dtheta para U = Jm(k r) * ang(m theta)
/******************************************************************************/
/* FUNCAO: polar_derivatives_U                                                */
/* DESCRICAO: Avalia derivadas radial e angular do potencial modal U em coordenadas polares. */
/* ENTRADA: R: double; m: int; kc: double; r: double; th: double; use_sin:    */
/* bool; dUdr: double &; dUdth_over_r: double &.                              */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void polar_derivatives_U(double R, int m, double kc, double r, double th,
                                bool use_sin,
                                double &dUdr, double &dUdth_over_r)
{
    // parte angular e sua derivada
    double ang, dang;
    if (!use_sin)
    {
        // cos(m th)
        ang  = std::cos(m * th);
        dang = -m * std::sin(m * th);
    }
    else
    {
        // sin(m th)
        ang  = std::sin(m * th);
        dang =  m * std::cos(m * th);
    }

    double kr = kc * r;

    // derivada radial: d/dr [Jm(k r)] = k Jm'(k r)
    dUdr = kc * Jm_prime(m, kr) * ang;

    // derivada angular: (1/r) d/dtheta [Jm(k r) ang] = (1/r) Jm(k r) dang
    if (r > 1e-14)
        dUdth_over_r = (Jm(m, kr) * dang) / r;
    else
        dUdth_over_r = 0.0;
}

// converte derivadas polares para gradiente cartesiano de U
/******************************************************************************/
/* FUNCAO: gradU_cart                                                         */
/* DESCRICAO: Converte gradiente de U de coordenadas polares para cartesianas (x,y).  */
/* ENTRADA: dUdr: double; dUdth_over_r: double; th: double; dUdx: double &;   */
/* dUdy: double &.                                                            */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void gradU_cart(double dUdr, double dUdth_over_r, double th,
                       double &dUdx, double &dUdy)
{
    // dUdx = cos th * dUdr - sin th * (1/r dUdth)
    // dUdy = sin th * dUdr + cos th * (1/r dUdth)
    double c = std::cos(th), s = std::sin(th);
    dUdx = c * dUdr - s * dUdth_over_r;
    dUdy = s * dUdr + c * dUdth_over_r;
}

// campo Ft = z x grad U = (-dU/dy, dU/dx)
/******************************************************************************/
/* FUNCAO: field_Ft_from_U                                                    */
/* DESCRICAO: Reconstrui campo transversal F_t a partir do potencial escalar U no guia circular. */
/* ENTRADA: dUdx: double; dUdy: double; Fx: double &; Fy: double &.           */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void field_Ft_from_U(double dUdx, double dUdy, double &Fx, double &Fy)
{
    Fx = -dUdy;
    Fy =  dUdx;
}

// Gauss de 2 pontos na aresta, retornando a componente tangencial media:
//   (1/L) * integral Ft.t ds
/******************************************************************************/
/* FUNCAO: edge_dof_from_field_gauss2                                         */
/* DESCRICAO: Integra o campo tangencial ao longo da aresta via Gauss-2 para obter grau de liberdade de aresta. */
/* ENTRADA: x0: double; y0: double; x1: double; y1: double; field_xy: auto.   */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double edge_dof_from_field_gauss2(
    double x0, double y0, double x1, double y1,
    auto field_xy) // funcao: field_xy(x,y)->pair(Fx,Fy)
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
    return I / L;
}

// DOFs analiticos de aresta para o caso circular:
//   y[dof] = (1/L) * integral Ft.t ds
/******************************************************************************/
/* FUNCAO: analytic_circle_edges_dofs                                         */
/* DESCRICAO: Calcula expressao analitica de referencia para comparacao e     */
/* validacao dos resultados numericos.                                        */
/* ENTRADA: mesh: const Mesh2D &; ed: const EdgeDofs &; R: double; m: int; p: */
/* int; roots_m: const std::vector<double> &; derivative_roots: bool; param:  */
/* // TE:true (J'); use_sin: TM:false (J) bool; y: std::vector<double> &.     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void analytic_circle_edges_dofs(
    const Mesh2D &mesh,
    const EdgeDofs &ed,
    double R,
    int m, int p,
    const std::vector<double> &roots_m,
    bool derivative_roots, // TE:true (J'), TM:false (J)
    bool use_sin,
    std::vector<double> &y)
{
    (void)derivative_roots; // usado apenas para documentar a origem das raizes

    double root = roots_m[(size_t)(p - 1)];
    double kc = root / R;

    y.assign((size_t)ed.ndof, 0.0);

    for (int eid = 0; eid < (int)ed.edges.size(); eid++)
    {
        int dof = ed.edge_to_dof[eid];
        if (dof < 0) continue; // eliminado (BC de TE)

        const auto &E = ed.edges[eid];
        const auto &n0 = mesh.nodes[E.n0];
        const auto &n1 = mesh.nodes[E.n1];

        auto field_xy = [&](double x, double y_)
        {
            double r = std::sqrt(x * x + y_ * y_);
            double th = std::atan2(y_, x);

            double dUdr, dUdth_over_r;
            polar_derivatives_U(R, m, kc, r, th, use_sin, dUdr, dUdth_over_r);

            double dUdx, dUdy;
            gradU_cart(dUdr, dUdth_over_r, th, dUdx, dUdy);

            double Fx, Fy;
            field_Ft_from_U(dUdx, dUdy, Fx, Fy);
            return std::pair<double,double>{Fx, Fy};
        };

        y[(size_t)dof] = edge_dof_from_field_gauss2(n0.x, n0.y, n1.x, n1.y, field_xy);
    }
}

// ---------- Casamento modal: TE e TM ----------
/******************************************************************************/
/* FUNCAO: match_circle_edge_mode_by_mass_correlation_TE                      */
/* DESCRICAO: Classifica modo de aresta circular TE por correlacao de massa.  */
/* ENTRADA: mesh: const Mesh2D &; R: double; T: const DenseMat &; ed: const   */
/* EdgeDofs &; Zcol: const std::vector<double> &; mode_idx: int; mmax: int;   */
/* pmax: int.                                                                 */
/* SAIDA: CircEdgeModeID.                                                     */
/******************************************************************************/
inline CircEdgeModeID match_circle_edge_mode_by_mass_correlation_TE(
    const Mesh2D &mesh,
    double R,
    const DenseMat &T,
    const EdgeDofs &ed,
    const std::vector<double> &Zcol,
    int mode_idx,
    int mmax = 6,
    int pmax = 6)
{
    std::vector<double> x;
    extract_edge_mode_from_Zcol(Zcol, ed.ndof, mode_idx, x);

    CircEdgeModeID best;
    best.rho = -1.0;

    for (int m = 0; m <= mmax; m++)
    {
        bool deriv = true; // TE usa J'(raiz)=0
        auto roots = bessel_roots(m, pmax, deriv);

        for (int p = 1; p <= pmax; p++)
        {
            std::vector<double> y_cos, y_sin;

            analytic_circle_edges_dofs(mesh, ed, R, m, p, roots, true, false, y_cos);
            double rho_cos = mass_correlation_abs(T, x, y_cos);

            double rho = rho_cos;
            bool use_sin = false;

            if (m > 0)
            {
                analytic_circle_edges_dofs(mesh, ed, R, m, p, roots, true, true, y_sin);
                double rho_sin = mass_correlation_abs(T, x, y_sin);
                if (rho_sin > rho_cos) { rho = rho_sin; use_sin = true; }
            }

            if (rho > best.rho)
            {
                double root = roots[(size_t)(p - 1)];
                best = {m, p, root, root / R, rho, use_sin};
            }
        }
    }

    return best;
}

/******************************************************************************/
/* FUNCAO: match_circle_edge_mode_by_mass_correlation_TM                      */
/* DESCRICAO: Classifica modo de aresta circular TM por correlacao de massa.  */
/* ENTRADA: mesh: const Mesh2D &; R: double; T: const DenseMat &; ed: const   */
/* EdgeDofs &; Zcol: const std::vector<double> &; mode_idx: int; mmax: int;   */
/* pmax: int.                                                                 */
/* SAIDA: CircEdgeModeID.                                                     */
/******************************************************************************/
inline CircEdgeModeID match_circle_edge_mode_by_mass_correlation_TM(
    const Mesh2D &mesh,
    double R,
    const DenseMat &T,
    const EdgeDofs &ed,
    const std::vector<double> &Zcol,
    int mode_idx,
    int mmax = 6,
    int pmax = 6)
{
    std::vector<double> x;
    extract_edge_mode_from_Zcol(Zcol, ed.ndof, mode_idx, x);

    CircEdgeModeID best;
    best.rho = -1.0;

    for (int m = 0; m <= mmax; m++)
    {
        bool deriv = false; // TM usa J(raiz)=0
        auto roots = bessel_roots(m, pmax, deriv);

        for (int p = 1; p <= pmax; p++)
        {
            std::vector<double> y_cos, y_sin;

            analytic_circle_edges_dofs(mesh, ed, R, m, p, roots, false, false, y_cos);
            double rho_cos = mass_correlation_abs(T, x, y_cos);

            double rho = rho_cos;
            bool use_sin = false;

            if (m > 0)
            {
                analytic_circle_edges_dofs(mesh, ed, R, m, p, roots, false, true, y_sin);
                double rho_sin = mass_correlation_abs(T, x, y_sin);
                if (rho_sin > rho_cos) { rho = rho_sin; use_sin = true; }
            }

            if (rho > best.rho)
            {
                double root = roots[(size_t)(p - 1)];
                best = {m, p, root, root / R, rho, use_sin};
            }
        }
    }

    return best;
}
