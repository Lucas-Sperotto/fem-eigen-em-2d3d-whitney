/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_rect_reference.hpp                             */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema misto vetorial+escalar para kc, separando blocos        */
/* transverso/longitudinal.                                                   */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.2, Eq.    */
/* (92).                                                                      */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

struct RectAnaMode
{
    int m = 0;
    int n = 0;
    double kc = 0.0;
};

// Familia analitica TE em guia retangular:
//   kc(m,n) = sqrt((m*pi/a)^2 + (n*pi/b)^2), except (0,0).
/******************************************************************************/
/* FUNCAO: analytic_rect_te                                                   */
/* DESCRICAO: Calcula expressao analitica de referencia para comparacao e     */
/* validacao dos resultados numericos.                                        */
/* ENTRADA: a: double; b: double; n_modes: int.                               */
/* SAIDA: std::vector<RectAnaMode>.                                           */
/******************************************************************************/
inline std::vector<RectAnaMode> analytic_rect_te(double a, double b, int n_modes)
{
    const double pi = 3.14159265358979323846;
    std::vector<RectAnaMode> out;

    constexpr int mmax = 30;
    constexpr int nmax = 30;
    for (int m = 0; m <= mmax; ++m)
    {
        for (int n = 0; n <= nmax; ++n)
        {
            if (m == 0 && n == 0)
                continue;
            const double kc = std::sqrt(std::pow(m * pi / a, 2) + std::pow(n * pi / b, 2));
            out.push_back({m, n, kc});
        }
    }

    std::sort(out.begin(), out.end(), [](const RectAnaMode &x, const RectAnaMode &y)
              { return x.kc < y.kc; });
    if ((int)out.size() > n_modes)
        out.resize((size_t)n_modes);
    return out;
}

// Familia analitica TM em guia retangular:
//   kc(m,n) = sqrt((m*pi/a)^2 + (n*pi/b)^2), com m>=1 e n>=1.
/******************************************************************************/
/* FUNCAO: analytic_rect_tm                                                   */
/* DESCRICAO: Calcula expressao analitica de referencia para comparacao e     */
/* validacao dos resultados numericos.                                        */
/* ENTRADA: a: double; b: double; n_modes: int.                               */
/* SAIDA: std::vector<RectAnaMode>.                                           */
/******************************************************************************/
inline std::vector<RectAnaMode> analytic_rect_tm(double a, double b, int n_modes)
{
    const double pi = 3.14159265358979323846;
    std::vector<RectAnaMode> out;

    constexpr int mmax = 30;
    constexpr int nmax = 30;
    for (int m = 1; m <= mmax; ++m)
    {
        for (int n = 1; n <= nmax; ++n)
        {
            const double kc = std::sqrt(std::pow(m * pi / a, 2) + std::pow(n * pi / b, 2));
            out.push_back({m, n, kc});
        }
    }

    std::sort(out.begin(), out.end(), [](const RectAnaMode &x, const RectAnaMode &y)
              { return x.kc < y.kc; });
    if ((int)out.size() > n_modes)
        out.resize((size_t)n_modes);
    return out;
}

// Mantenha este formato de saida estavel, pois scripts/validate_2d_22.py o interpreta.
/******************************************************************************/
/* FUNCAO: print_rect_compare_table                                           */
/* DESCRICAO: Gera saida didatica dos resultados para inspecao no terminal ou */
/* em arquivo.                                                                */
/* ENTRADA: title: const std::string &; ana: const std::vector<RectAnaMode> &;*/
/* calc: const std::vector<double> &; rho: double.                            */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_rect_compare_table(
    const std::string &title,
    const std::vector<RectAnaMode> &ana,
    const std::vector<double> &calc,
    double rho)
{
    std::cout << "\n" << title << "\n";
    std::cout << "rho(filter kc^2) = " << rho << "\n";
    std::cout << " #   mode     kc(ana)        kc(calc)       err(%)\n";

    const int n = std::min((int)ana.size(), (int)calc.size());
    for (int i = 0; i < n; ++i)
    {
        const double kc_a = ana[(size_t)i].kc;
        const double kc_c = calc[(size_t)i];
        const double errp = 100.0 * (kc_c - kc_a) / kc_a;

        std::cout << " " << (i + 1)
                  << "   (" << ana[(size_t)i].m << "," << ana[(size_t)i].n << ")"
                  << "  " << kc_a
                  << "   " << kc_c
                  << "   " << errp
                  << "\n";
    }
}
