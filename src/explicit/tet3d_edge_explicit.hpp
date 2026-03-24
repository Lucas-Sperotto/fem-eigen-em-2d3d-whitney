/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/explicit/tet3d_edge_explicit.hpp                              */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Formas fechadas (closed-form) do elemento de aresta tetraedrico */
/* 3D da Secao 3.1.                                                           */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Eq. (162)-(182).    */
/*****************************************************************************/
/* Observacao: Este cabecalho materializa os coeficientes geometricos do      */
/* tetraedro e os termos I1..I10 usados na massa vetorial explicita.          */
/*****************************************************************************/

#pragma once

#include "edge3d/edge3d_basis.hpp"

#include <array>

namespace explicit_tet3d
{

struct Tet3DEdgeClosedFormData
{
  double volume = 0.0;
  double xtet = 0.0;
  double ytet = 0.0;
  double ztet = 0.0;
  double sum_x2 = 0.0;
  double sum_y2 = 0.0;
  double sum_z2 = 0.0;
  double sum_xy = 0.0;
  double sum_xz = 0.0;
  double sum_yz = 0.0;
  std::array<double, 6> edge_len{};
  std::array<double, 6> A_x{};
  std::array<double, 6> B_x{};
  std::array<double, 6> C_x{};
  std::array<double, 6> A_y{};
  std::array<double, 6> B_y{};
  std::array<double, 6> C_y{};
  std::array<double, 6> A_z{};
  std::array<double, 6> B_z{};
  std::array<double, 6> C_z{};
};

/******************************************************************************/
/* FUNCAO: tet3d_edge_closed_form_data                                        */
/* DESCRICAO: Construi todos os coeficientes geometricos usados nas formas    */
/* fechadas 3D. Em particular, reaproveita a Eq. (162) para formar as Eq.     */
/* (164)-(172) e prepara os momentos geometricos do tetraedro. Esses dados    */
/* alimentam diretamente as matrizes locais Eq. (181)-(182) e, por            */
/* consequencia, a montagem global Eq. (178).                                 */
/* ENTRADA: tg: const TetGeomEdge &.                                          */
/* SAIDA: Tet3DEdgeClosedFormData.                                            */
/******************************************************************************/
inline Tet3DEdgeClosedFormData tet3d_edge_closed_form_data(const TetGeomEdge &tg)
{
  Tet3DEdgeClosedFormData data;
  data.volume = tg.V;
  data.edge_len = tg.L;

  for (int i = 0; i < 4; ++i)
  {
    data.xtet += tg.X[i].x;
    data.ytet += tg.X[i].y;
    data.ztet += tg.X[i].z;
    data.sum_x2 += tg.X[i].x * tg.X[i].x;
    data.sum_y2 += tg.X[i].y * tg.X[i].y;
    data.sum_z2 += tg.X[i].z * tg.X[i].z;
    data.sum_xy += tg.X[i].x * tg.X[i].y;
    data.sum_xz += tg.X[i].x * tg.X[i].z;
    data.sum_yz += tg.X[i].y * tg.X[i].z;
  }

  data.xtet /= 4.0;
  data.ytet /= 4.0;
  data.ztet /= 4.0;

  constexpr std::array<std::array<int, 2>, 6> kLocalEdges = {{
      {{0, 1}},
      {{0, 2}},
      {{0, 3}},
      {{1, 2}},
      {{1, 3}},
      {{2, 3}},
  }};

  for (int m = 0; m < 6; ++m)
  {
    const int i = kLocalEdges[m][0];
    const int j = kLocalEdges[m][1];
    const TetSimplexCoeff3D &li = tg.lambda_coeff[i];
    const TetSimplexCoeff3D &lj = tg.lambda_coeff[j];

    // Eq. (164)-(172): coeficientes da base de Whitney 3D apos substituir
    // a Eq. (162) na Eq. (163).
    data.A_x[m] = li.a * lj.b - lj.a * li.b; // Eq. (164)
    data.B_x[m] = li.c * lj.b - lj.c * li.b; // Eq. (165)
    data.C_x[m] = li.d * lj.b - lj.d * li.b; // Eq. (166)

    data.A_y[m] = li.a * lj.c - lj.a * li.c; // Eq. (167)
    data.B_y[m] = li.b * lj.c - lj.b * li.c; // Eq. (168)
    data.C_y[m] = li.d * lj.c - lj.d * li.c; // Eq. (169)

    data.A_z[m] = li.a * lj.d - lj.a * li.d; // Eq. (170)
    data.B_z[m] = li.b * lj.d - lj.b * li.d; // Eq. (171)
    data.C_z[m] = li.c * lj.d - lj.c * li.d; // Eq. (172)
  }

  return data;
}

/******************************************************************************/
/* FUNCAO: tet3d_edge_i_terms_eq_182                                          */
/* DESCRICAO: Calcula os dez termos auxiliares I1..I10 usados na Eq. (182)   */
/* para a massa vetorial do tetraedro. Esses termos entram no bloco local que */
/* depois e acumulado no sistema global da Eq. (178).                         */
/* ENTRADA: data: const Tet3DEdgeClosedFormData &; m: int; n: int.            */
/* SAIDA: std::array<double, 10>.                                             */
/******************************************************************************/
inline std::array<double, 10> tet3d_edge_i_terms_eq_182(
    const Tet3DEdgeClosedFormData &data,
    int m,
    int n)
{
  const double xy_moment = (data.sum_xy + 16.0 * data.xtet * data.ytet) / 20.0;
  const double yz_moment = (data.sum_yz + 16.0 * data.ytet * data.ztet) / 20.0;
  const double xz_moment = (data.sum_xz + 16.0 * data.xtet * data.ztet) / 20.0;
  const double x2_moment = (data.sum_x2 + 16.0 * data.xtet * data.xtet) / 20.0;
  const double y2_moment = (data.sum_y2 + 16.0 * data.ytet * data.ytet) / 20.0;
  const double z2_moment = (data.sum_z2 + 16.0 * data.ztet * data.ztet) / 20.0;

  return {{
      // Eq. (182): I1
      data.A_x[m] * data.A_x[n] +
          data.A_y[m] * data.A_y[n] +
          data.A_z[m] * data.A_z[n],
      // Eq. (182): I2
      (data.A_y[m] * data.B_y[n] + data.A_y[n] * data.B_y[m] +
       data.A_z[m] * data.B_z[n] + data.A_z[n] * data.B_z[m]) *
          data.xtet,
      // Eq. (182): I3
      (data.A_x[m] * data.B_x[n] + data.A_x[n] * data.B_x[m] +
       data.A_z[m] * data.C_z[n] + data.A_z[n] * data.C_z[m]) *
          data.ytet,
      // Eq. (182): I4
      (data.A_x[m] * data.C_x[n] + data.A_x[n] * data.C_x[m] +
       data.A_y[m] * data.C_y[n] + data.A_y[n] * data.C_y[m]) *
          data.ztet,
      // Eq. (182): I5
      (data.B_z[m] * data.C_z[n] + data.B_z[n] * data.C_z[m]) * xy_moment,
      // Eq. (182): I6
      (data.B_x[m] * data.C_x[n] + data.B_x[n] * data.C_x[m]) * yz_moment,
      // Eq. (182): I7
      (data.B_y[m] * data.C_y[n] + data.B_y[n] * data.C_y[m]) * xz_moment,
      // Eq. (182): I8
      (data.B_y[m] * data.B_y[n] + data.B_z[m] * data.B_z[n]) * x2_moment,
      // Eq. (182): I9
      (data.B_x[m] * data.B_x[n] + data.C_z[m] * data.C_z[n]) * y2_moment,
      // Eq. (182): I10
      (data.C_x[m] * data.C_x[n] + data.C_y[m] * data.C_y[n]) * z2_moment,
  }};
}

/******************************************************************************/
/* FUNCAO: tet3d_edge_curlcurl_entry_eq_181                                   */
/* DESCRICAO: Calcula uma entrada local do bloco curl-curl 3D usando a forma  */
/* fechada obtida apos substituir as Eq. (163)-(172) na Eq. (181). Este bloco */
/* local e uma das parcelas que, apos assembleia global, formam a Eq. (178).  */
/* ENTRADA: data: const Tet3DEdgeClosedFormData &; m: int; n: int.            */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double tet3d_edge_curlcurl_entry_eq_181(
    const Tet3DEdgeClosedFormData &data,
    int m,
    int n)
{
  const double scale =
      (data.edge_len[m] * data.edge_len[n]) /
      (324.0 * data.volume * data.volume * data.volume);

  // Eq. (181): forma fechada do bloco curl-curl. Nesta implementacao usamos a
  // convencao algebrica diretamente induzida pelas Eq. (164)-(172), o que
  // leva ao somatorio C_zm*C_zn + C_xm*C_xn + B_ym*B_yn.
  return scale *
         (data.C_z[m] * data.C_z[n] +
          data.C_x[m] * data.C_x[n] +
          data.B_y[m] * data.B_y[n]);
}

/******************************************************************************/
/* FUNCAO: tet3d_edge_mass_entry_eq_182                                       */
/* DESCRICAO: Calcula uma entrada local do bloco de massa vetorial 3D pela    */
/* forma fechada da Eq. (182), ja agregando os termos I1..I10. Este bloco     */
/* local tambem e acumulado depois na Eq. (178).                              */
/* ENTRADA: data: const Tet3DEdgeClosedFormData &; m: int; n: int.            */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double tet3d_edge_mass_entry_eq_182(
    const Tet3DEdgeClosedFormData &data,
    int m,
    int n)
{
  const auto I = tet3d_edge_i_terms_eq_182(data, m, n);
  const double scale =
      (data.edge_len[m] * data.edge_len[n]) /
      (1296.0 * data.volume * data.volume * data.volume);

  double sum = 0.0;
  for (double term : I)
    sum += term;
  return scale * sum;
}

/******************************************************************************/
/* FUNCAO: tet3d_edge_closed_form_eq_181_182                                  */
/* DESCRICAO: Monta os blocos locais Sel/Tel do tetraedro 3D usando as formas */
/* fechadas Eq. (181) e Eq. (182). Esses blocos sao a versao closed-form do   */
/* elemento local que alimenta a montagem global da Eq. (178).                */
/* ENTRADA: tg: const TetGeomEdge &; inv_mu_r: double; eps_r: double; Sel:    */
/* double[6][6]; Tel: double[6][6].                                           */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void tet3d_edge_closed_form_eq_181_182(
    const TetGeomEdge &tg,
    double inv_mu_r,
    double eps_r,
    double Sel[6][6],
    double Tel[6][6])
{
  const Tet3DEdgeClosedFormData data = tet3d_edge_closed_form_data(tg);
  for (int m = 0; m < 6; ++m)
  {
    for (int n = 0; n < 6; ++n)
    {
      Sel[m][n] = inv_mu_r * tet3d_edge_curlcurl_entry_eq_181(data, m, n);
      Tel[m][n] = eps_r * tet3d_edge_mass_entry_eq_182(data, m, n);
    }
  }
}

} // namespace explicit_tet3d
