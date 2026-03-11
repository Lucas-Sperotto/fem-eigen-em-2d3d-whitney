/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/dense.hpp                                                */
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
#include <vector>
#include <cassert>

struct DenseMat
{
  int n = 0;
  std::vector<double> a; // row-major n x n

  DenseMat() = default;
/******************************************************************************/
/* FUNCAO: DenseMat                                                           */
/* DESCRICAO: Construtor da matriz densa quadrada, inicializando dimensao e armazenamento com zeros. */
/* ENTRADA: n_: int.                                                          */
/* SAIDA: sem retorno explicito (construtor/destrutor).                       */
/******************************************************************************/
  explicit DenseMat(int n_) : n(n_), a((size_t)n_ * n_, 0.0) {}

/******************************************************************************/
/* FUNCAO: operator                                                           */
/* DESCRICAO: Acessa elemento (i,j) com referencia mutavel para escrita na    */
/* matriz densa (row-major).                                                  */
/* ENTRADA: i: int; j: int.                                                   */
/* SAIDA: double &.                                                           */
/******************************************************************************/
  double &operator()(int i, int j)
  {
    assert(i >= 0 && i < n && j >= 0 && j < n);
    return a[(size_t)i * n + j];
  }
/******************************************************************************/
/* FUNCAO: operator                                                           */
/* DESCRICAO: Acessa elemento (i,j) apenas para leitura na matriz densa       */
/* (sobrecarga const).                                                        */
/* ENTRADA: i: int; j: int.                                                   */
/* SAIDA: double.                                                             */
/******************************************************************************/
  double operator()(int i, int j) const
  {
    assert(i >= 0 && i < n && j >= 0 && j < n);
    return a[(size_t)i * n + j];
  }
};
