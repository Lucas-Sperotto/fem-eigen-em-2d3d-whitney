/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/sparse_sym.hpp                                           */
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
#include <unordered_map>
#include <vector>
#include <stdexcept>

// Armazenamento esparso simetrico da parte triangular inferior.
struct SparseSymMat
{
  int n = 0;
  std::vector<std::unordered_map<int, double>> row_to_lower;

  SparseSymMat() = default;
/******************************************************************************/
/* FUNCAO: SparseSymMat                                                       */
/* DESCRICAO: Construtor da matriz simetrica esparsa em formato de dicionario da parte triangular inferior. */
/* ENTRADA: n_: int.                                                          */
/* SAIDA: sem retorno explicito (construtor/destrutor).                       */
/******************************************************************************/
  explicit SparseSymMat(int n_) : n(n_), row_to_lower((size_t)n_) {}

/******************************************************************************/
/* FUNCAO: add                                                                */
/* DESCRICAO: Acumula valor em (i,j) explorando simetria, gravando somente na */
/* metade inferior para reduzir memoria.                                      */
/* ENTRADA: i: int; j: int; v: double.                                        */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
  void add(int i, int j, double v)
  {
    if (i < 0 || i >= n || j < 0 || j >= n)
      throw std::runtime_error("Indice fora do intervalo em SparseSymMat::add.");
    if (i < j)
      std::swap(i, j);
    row_to_lower[(size_t)i][j] += v;
  }

/******************************************************************************/
/* FUNCAO: nnz_lower                                                          */
/* DESCRICAO: Conta quantos termos nao nulos existem na parte triangular inferior da matriz esparsa. */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: std::size_t.                                                        */
/******************************************************************************/
  std::size_t nnz_lower() const
  {
    std::size_t nnz = 0;
    for (const auto &row : row_to_lower)
      nnz += row.size();
    return nnz;
  }

/******************************************************************************/
/* FUNCAO: to_dense                                                           */
/* DESCRICAO: Expande a representacao esparsa simetrica para matriz densa, util para depuracao e validacao. */
/* ENTRADA: sem parametros.                                                   */
/* SAIDA: DenseMat.                                                           */
/******************************************************************************/
  DenseMat to_dense() const
  {
    DenseMat out(n);
    for (int i = 0; i < n; ++i)
    {
      for (const auto &kv : row_to_lower[(size_t)i])
      {
        const int j = kv.first;
        const double v = kv.second;
        out(i, j) += v;
        if (i != j)
          out(j, i) += v;
      }
    }
    return out;
  }
};
