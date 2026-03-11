/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge/mode_match_rect_edge_tm.hpp                              */
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

struct RectEdgeModeID_TM
{
    int m = 0, n = 0;
    double kc_ana = 0;
    double rho = 0;
};

/******************************************************************************/
/* FUNCAO: kc_rect_analytic                                                   */
/* DESCRICAO: Retorna o valor analitico de k_c para o modo retangular (m,n), conforme separacao de variaveis. */
/* ENTRADA: a: double; b: double; m: int; n: int.                             */
/* SAIDA: double.                                                             */
/******************************************************************************/
static inline double kc_rect_analytic(double a, double b, int m, int n)
{
    double kx = (m * M_PI) / a;
    double ky = (n * M_PI) / b;
    return std::sqrt(kx * kx + ky * ky);
}

// TM: psi = sin(kx x) sin(ky y)
// Ht = z x grad(psi) = (-dpsi/dy, dpsi/dx)
/******************************************************************************/
/* FUNCAO: tm_field_rect                                                      */
/* DESCRICAO: Avalia o campo TM analitico no guia retangular para comparacao com resultados FEM. */
/* ENTRADA: a: double; b: double; m: int; n: int; x: double; y: double; Hx:   */
/* double &; Hy: double &.                                                    */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static inline void tm_field_rect(double a, double b, int m, int n,
                                 double x, double y,
                                 double &Hx, double &Hy)
{
    double kx = (m * M_PI) / a;
    double ky = (n * M_PI) / b;

    // dpsi/dx = kx cos(kx x) sin(ky y)
    // dpsi/dy = ky sin(kx x) cos(ky y)
    // Hx = -dpsi/dy
    // Hy =  dpsi/dx
    Hx = -ky * std::sin(kx * x) * std::cos(ky * y);
    Hy = kx * std::cos(kx * x) * std::sin(ky * y);
}

/******************************************************************************/
/* FUNCAO: edge_dof_tm_gauss2                                                 */
/* DESCRICAO: Projeta campo TM analitico nos graus de liberdade de aresta via quadratura Gauss-2. */
/* ENTRADA: a: double; b: double; m: int; n: int; x0: double; y0: double; x1: */
/* double; y1: double.                                                        */
/* SAIDA: double.                                                             */
/******************************************************************************/
static inline double edge_dof_tm_gauss2(double a, double b, int m, int n,
                                        double x0, double y0, double x1, double y1)
{
    double dx = x1 - x0, dy = y1 - y0;
    double L = std::sqrt(dx * dx + dy * dy);
    if (L <= 0)
        return 0.0;
    double tx = dx / L, ty = dy / L;

    double s1 = 0.5 * L * (1.0 - 1.0 / std::sqrt(3.0));
    double s2 = 0.5 * L * (1.0 + 1.0 / std::sqrt(3.0));
    double w = 0.5 * L;

    auto sample = [&](double s)
    {
        double x = x0 + tx * s;
        double y = y0 + ty * s;
        double Hx, Hy;
        tm_field_rect(a, b, m, n, x, y, Hx, Hy);
        return Hx * tx + Hy * ty;
    };

    double I = w * sample(s1) + w * sample(s2);
    return I / L;
}

/******************************************************************************/
/* FUNCAO: mass_inner                                                         */
/* DESCRICAO: Calcula produto interno ponderado pela matriz de massa, base    */
/* para correlacao e normalizacao modal.                                      */
/* ENTRADA: T: const DenseMat &; x: const std::vector<double> &; y: const     */
/* std::vector<double> &.                                                     */
/* SAIDA: double.                                                             */
/******************************************************************************/
static inline double mass_inner(const DenseMat &T, const std::vector<double> &x, const std::vector<double> &y)
{
    int n = T.n;
    double s = 0.0;
    for (int i = 0; i < n; i++)
    {
        double row = 0.0;
        for (int j = 0; j < n; j++)
            row += T(i, j) * y[j];
        s += x[i] * row;
    }
    return s;
}

/******************************************************************************/
/* FUNCAO: mass_norm                                                          */
/* DESCRICAO: Calcula a norma energetica induzida pela matriz de massa.       */
/* ENTRADA: T: const DenseMat &; x: const std::vector<double> &.              */
/* SAIDA: double.                                                             */
/******************************************************************************/
static inline double mass_norm(const DenseMat &T, const std::vector<double> &x)
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
static inline double mass_correlation_abs(const DenseMat &T,
                                          const std::vector<double> &x,
                                          const std::vector<double> &y)
{
    double num = std::abs(mass_inner(T, x, y));
    double den = mass_norm(T, x) * mass_norm(T, y);
    if (den <= 0)
        return 0.0;
    return num / den;
}

/******************************************************************************/
/* FUNCAO: extract_edge_mode_from_Zcol                                        */
/* DESCRICAO: Extrai um modo no espaco de arestas a partir da matriz de       */
/* autovetores column-major.                                                  */
/* ENTRADA: Zcol: const std::vector<double> &; ndof: int; mode_idx: int; x:   */
/* std::vector<double> &.                                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static inline void extract_edge_mode_from_Zcol(const std::vector<double> &Zcol,
                                               int ndof, int mode_idx,
                                               std::vector<double> &x)
{
    x.assign(ndof, 0.0);
    auto idx = [&](int i, int j)
    { return (size_t)j * ndof + i; };
    for (int i = 0; i < ndof; i++)
        x[i] = Zcol[idx(i, mode_idx)];
}

/******************************************************************************/
/* FUNCAO: analytic_tm_edges_rect                                             */
/* DESCRICAO: Calcula expressao analitica de referencia para comparacao e     */
/* validacao dos resultados numericos.                                        */
/* ENTRADA: mesh: const Mesh2D &; ed: const EdgeDofs &; a: double; b: double; */
/* m: int; n: int; y: std::vector<double> &.                                  */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
static inline void analytic_tm_edges_rect(const Mesh2D &mesh,
                                          const EdgeDofs &ed,
                                          double a, double b,
                                          int m, int n,
                                          std::vector<double> &y)
{
    y.assign(ed.ndof, 0.0);
    for (int eid = 0; eid < (int)ed.edges.size(); eid++)
    {
        int dof = ed.edge_to_dof[eid];
        if (dof < 0)
            continue; // (para TM normalmente ninguém é -1)
        const auto &E = ed.edges[eid];
        const auto &n0 = mesh.nodes[E.n0];
        const auto &n1 = mesh.nodes[E.n1];
        y[dof] = edge_dof_tm_gauss2(a, b, m, n, n0.x, n0.y, n1.x, n1.y);
    }
}

/******************************************************************************/
/* FUNCAO: match_rect_edge_mode_by_mass_correlation_tm                        */
/* DESCRICAO: Classifica modo de aresta retangular TM por correlacao de massa.*/
/* ENTRADA: mesh: const Mesh2D &; a: double; b: double; T: const DenseMat &;  */
/* ed: const EdgeDofs &; Zcol: const std::vector<double> &; mode_idx: int;    */
/* mmax: int; nmax: int.                                                      */
/* SAIDA: RectEdgeModeID_TM.                                                  */
/******************************************************************************/
static inline RectEdgeModeID_TM match_rect_edge_mode_by_mass_correlation_tm(
    const Mesh2D &mesh,
    double a, double b,
    const DenseMat &T,
    const EdgeDofs &ed,
    const std::vector<double> &Zcol,
    int mode_idx,
    int mmax = 8, int nmax = 8)
{
    std::vector<double> x;
    extract_edge_mode_from_Zcol(Zcol, ed.ndof, mode_idx, x);

    RectEdgeModeID_TM best;
    best.rho = -1.0;

    for (int m = 1; m <= mmax; m++)
    {
        for (int n = 1; n <= nmax; n++)
        {
            std::vector<double> y;
            analytic_tm_edges_rect(mesh, ed, a, b, m, n, y);
            double rho = mass_correlation_abs(T, x, y);
            if (rho > best.rho)
            {
                best = {m, n, kc_rect_analytic(a, b, m, n), rho};
            }
        }
    }
    return best;
}
