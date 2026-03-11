/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/fem3d/fem3d_compare.hpp                                       */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Infraestrutura de comparacao dos casos 3D com tabelas do artigo.*/
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 3.1.          */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fem3d
{
// Linha usada nas tabelas de comparacao com o artigo:
// - campo 'analytical': valor analitico em forma fechada
// - ref_paper: valor reportado na implementacao de referencia do artigo (FEM3D1/ref.17)
struct RefRow
{
  const char *mode = "";
  double analytical = 0.0;
  double ref_paper = 0.0;
};

/******************************************************************************/
/* FUNCAO: first_positive_k0                                                  */
/* DESCRICAO: Extrai o primeiro autovalor positivo (k0) da lista numerica, ignorando residuos nao fisicos. */
/* Essa conversao e aplicada apos resolver Sx=lambdaTx, com lambda=k0^2.      */
/* ENTRADA: lambda: const std::vector<double> &; nmax: int; tol: double.      */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> first_positive_k0(
    const std::vector<double> &lambda,
    int nmax,
    double tol = 1e-10)
{
  std::vector<double> out;
  out.reserve((size_t)nmax);
  for (double l : lambda)
  {
    if (l <= tol)
      continue;
    out.push_back(std::sqrt(l));
    if ((int)out.size() >= nmax)
      break;
  }
  return out;
}

// Casamento guloso com agrupamento explicito para valores analiticos iguais.
// Isso mantem familias degeneradas (mesmo k0 analitico) tratadas de forma consistente.
/******************************************************************************/
/* FUNCAO: match_by_reference_with_degeneracy                                 */
/* DESCRICAO: Associa resultados FEM a referencias analiticas usando criterio */
/* de proximidade, preservando grupos degenerados reportados nas Tabelas      */
/* 12-15.                                                                     */
/* ENTRADA: ref: const std::vector<RefRow> &; k0_all: const                   */
/* std::vector<double> &; scan_limit: int.                                    */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> match_by_reference_with_degeneracy(
    const std::vector<RefRow> &ref,
    const std::vector<double> &k0_all,
    int scan_limit)
{
  struct Group
  {
    double target = 0.0;
    std::vector<int> iref;
  };

  const int ns = std::min((int)k0_all.size(), scan_limit);
  std::vector<double> kscan;
  kscan.reserve((size_t)ns);
  for (int i = 0; i < ns; ++i)
    kscan.push_back(k0_all[(size_t)i]);

  std::vector<Group> groups;
  constexpr double kDegTol = 1e-9;
  for (int i = 0; i < (int)ref.size(); ++i)
  {
    bool placed = false;
    for (auto &g : groups)
    {
      if (std::abs(g.target - ref[(size_t)i].analytical) < kDegTol)
      {
        g.iref.push_back(i);
        placed = true;
        break;
      }
    }
    if (!placed)
      groups.push_back({ref[(size_t)i].analytical, {i}});
  }

  std::vector<int> used((size_t)ns, 0);
  std::vector<double> out(ref.size(), std::nan(""));

  for (const auto &g : groups)
  {
    struct Cand
    {
      double err = 0.0;
      int ik = -1;
    };
    std::vector<Cand> cand;
    cand.reserve((size_t)ns);
    for (int ik = 0; ik < ns; ++ik)
    {
      if (used[(size_t)ik])
        continue;
      const double e = std::abs(kscan[(size_t)ik] - g.target) / g.target;
      cand.push_back({e, ik});
    }
    std::sort(cand.begin(), cand.end(), [](const Cand &a, const Cand &b)
              { return a.err < b.err; });

    const int take = std::min((int)g.iref.size(), (int)cand.size());
    std::vector<int> picked;
    picked.reserve((size_t)take);
    for (int i = 0; i < take; ++i)
      picked.push_back(cand[(size_t)i].ik);
    std::sort(picked.begin(), picked.end(), [&](int a, int b)
              { return kscan[(size_t)a] < kscan[(size_t)b]; });

    for (int i = 0; i < take; ++i)
    {
      const int ir = g.iref[(size_t)i];
      const int ik = picked[(size_t)i];
      used[(size_t)ik] = 1;
      out[(size_t)ir] = kscan[(size_t)ik];
    }
  }

  return out;
}

/******************************************************************************/
/* FUNCAO: print_table_compare                                                */
/* DESCRICAO: Gera saida didatica dos resultados para inspecao no terminal ou */
/* em arquivo.                                                                */
/* ENTRADA: title: const std::string &; ref: const std::vector<RefRow> &;     */
/* k0_fem: const std::vector<double> &.                                       */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_table_compare(
    const std::string &title,
    const std::vector<RefRow> &ref,
    const std::vector<double> &k0_fem)
{
  std::cout << "\n" << title << "\n";
  std::cout << "idx  mode      k0_ana     k0_fem     err_ana(%)   ref_paper  err_ref(%)\n";
  for (int i = 0; i < (int)ref.size(); ++i)
  {
    const double kf = (i < (int)k0_fem.size()) ? k0_fem[(size_t)i] : std::nan("");
    const double ea = std::isfinite(kf) ? 100.0 * std::abs(kf - ref[(size_t)i].analytical) / ref[(size_t)i].analytical : std::nan("");
    const double er = std::isfinite(kf) ? 100.0 * std::abs(kf - ref[(size_t)i].ref_paper) / ref[(size_t)i].ref_paper : std::nan("");
    std::cout << std::setw(3) << (i + 1) << "  "
              << std::setw(6) << ref[(size_t)i].mode << "  "
              << std::setw(8) << std::fixed << std::setprecision(3) << ref[(size_t)i].analytical << "  "
              << std::setw(8) << std::setprecision(3) << kf << "  "
              << std::setw(10) << std::setprecision(3) << ea << "  "
              << std::setw(8) << std::setprecision(3) << ref[(size_t)i].ref_paper << "  "
              << std::setw(10) << std::setprecision(3) << er << "\n";
  }
}

} // namespace fem3d
