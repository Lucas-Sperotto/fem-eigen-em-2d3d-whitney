/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/lapack_eig.hpp                                           */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios numericos para autovalor generalizado e operacoes   */
/* matriciais.                                                                */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Apoio numerico as   */
/* Secoes 2.x e 3.1.                                                          */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "dense.hpp"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <lapacke.h>

struct GenEigResult
{
  std::vector<double> w;    // autovalores
  std::vector<double> Zcol; // autovetores em column-major (cada coluna = vetor)
  int n = 0;
};

/******************************************************************************/
/* FUNCAO: generalized_eigs_sym_vec                                           */
/* DESCRICAO: Resolve problema de autovalor generalizado simetrico e retorna  */
/* autovalores/autovetores.                                                   */
/* ENTRADA: S: DenseMat; T: DenseMat.                                         */
/* SAIDA: GenEigResult.                                                       */
/******************************************************************************/
inline GenEigResult generalized_eigs_sym_vec(DenseMat S, DenseMat T)
{
  const int n = S.n;
  if (T.n != n)
    throw std::runtime_error("Dimensoes diferentes em S e T.");

  std::vector<double> A((size_t)n * n), B((size_t)n * n), w((size_t)n);

  auto idx_col = [&](int i, int j)
  { return (size_t)j * n + i; };
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      A[idx_col(i, j)] = S(i, j);
      B[idx_col(i, j)] = T(i, j);
    }
  }

  // jobz='V' => autovalores + autovetores (A vira Z)
  int info = LAPACKE_dsygv(LAPACK_COL_MAJOR, 1, 'V', 'U', n, A.data(), n, B.data(), n, w.data());
  if (info != 0)
  {
    throw std::runtime_error("LAPACKE_dsygv falhou. info=" + std::to_string(info));
  }

  // A agora contem Z (autovetores)
  return {w, A, n};
}

struct GenEigGeneralResult
{
  std::vector<double> lambda_re; // real(lambda)
  std::vector<double> lambda_im; // imag(lambda)
  std::vector<double> beta;      // generalized denominator
  std::vector<double> VRcol;     // right eigenvectors (column-major)
  int n = 0;
};

/******************************************************************************/
/* FUNCAO: generalized_eigs_real_vec                                          */
/* DESCRICAO: Resolve problema de autovalor generalizado real nao simetrico e */
/* retorna espectro completo.                                                 */
/* ENTRADA: Arow: DenseMat; Brow: DenseMat.                                   */
/* SAIDA: GenEigGeneralResult.                                                */
/******************************************************************************/
inline GenEigGeneralResult generalized_eigs_real_vec(DenseMat Arow, DenseMat Brow)
{
  const int n = Arow.n;
  if (Brow.n != n)
    throw std::runtime_error("Dimensoes diferentes em A e B.");

  std::vector<double> A((size_t)n * n), B((size_t)n * n);
  std::vector<double> alphar((size_t)n), alphai((size_t)n), beta((size_t)n);
  std::vector<double> VR((size_t)n * n);

  // row-major DenseMat -> col-major LAPACK
  auto idx_col = [&](int i, int j)
  { return (size_t)j * n + i; };
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      A[idx_col(i, j)] = Arow(i, j);
      B[idx_col(i, j)] = Brow(i, j);
    }
  }

  int info = LAPACKE_dggev(
      LAPACK_COL_MAJOR,
      'N', // no VL
      'V', // compute VR
      n,
      A.data(),
      n,
      B.data(),
      n,
      alphar.data(),
      alphai.data(),
      beta.data(),
      nullptr,
      n,
      VR.data(),
      n);

  if (info != 0)
  {
    throw std::runtime_error("LAPACKE_dggev falhou. info=" + std::to_string(info));
  }

  std::vector<double> lambda_re((size_t)n), lambda_im((size_t)n);
  for (int i = 0; i < n; ++i)
  {
    if (std::abs(beta[i]) < 1e-14)
    {
      lambda_re[i] = std::numeric_limits<double>::infinity();
      lambda_im[i] = std::numeric_limits<double>::infinity();
    }
    else
    {
      lambda_re[i] = alphar[i] / beta[i];
      lambda_im[i] = alphai[i] / beta[i];
    }
  }

  return {lambda_re, lambda_im, beta, VR, n};
}
