/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/block_ops.hpp                                            */
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

// bloco-diagonal: [A 0; 0 B]
/******************************************************************************/
/* FUNCAO: block_diag                                                         */
/* DESCRICAO: Monta uma matriz bloco-diagonal combinando dois blocos quadrados*/
/* independentes.                                                             */
/* ENTRADA: A: const DenseMat&; B: const DenseMat&.                           */
/* SAIDA: DenseMat.                                                           */
/******************************************************************************/
inline DenseMat block_diag(const DenseMat& A, const DenseMat& B)
{
    const int nA = A.n;
    const int nB = B.n;
    DenseMat C(nA + nB);

    // A
    for (int i = 0; i < nA; ++i)
        for (int j = 0; j < nA; ++j)
            C(i, j) = A(i, j);

    // B
    for (int i = 0; i < nB; ++i)
        for (int j = 0; j < nB; ++j)
            C(nA + i, nA + j) = B(i, j);

    return C;
}
