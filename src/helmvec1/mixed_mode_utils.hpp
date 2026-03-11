/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_mode_utils.hpp                                 */
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
#include "core/lapack_eig.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// Separa autovetores generalizados pelo bloco dominante de energia.
//
// Em sistemas mistos no formato da Eq. (92), o autovetor e:
//   x = [block0 ; block1]
//
// Para cada autopar (lambda, x), esta rotina calcula ||block0||^2 e
// ||block1||^2, classifica pelo bloco dominante e armazena
// k = sqrt(lambda) na lista correspondente.
/******************************************************************************/
/* FUNCAO: split_modes_by_block_energy                                        */
/* DESCRICAO: Separa modos por predominancia energetica entre blocos de um    */
/* sistema misto.                                                             */
/* ENTRADA: res: const GenEigResult &; n_block0: int; n_block1: int;          */
/* lambda_min: double; k_block0: std::vector<double> &; k_block1:             */
/* std::vector<double> &.                                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void split_modes_by_block_energy(
    const GenEigResult &res,
    int n_block0,
    int n_block1,
    double lambda_min,
    std::vector<double> &k_block0,
    std::vector<double> &k_block1)
{
    k_block0.clear();
    k_block1.clear();

    const int n = res.n;
    auto vec_val = [&](int i, int k) -> double
    {
        return res.Zcol[(size_t)k * n + i];
    };

    for (int k = 0; k < n; ++k)
    {
        const double lam = res.w[(size_t)k];
        if (lam < lambda_min)
            continue;

        double e0 = 0.0, e1 = 0.0;
        for (int i = 0; i < n_block0; ++i)
        {
            const double v = vec_val(i, k);
            e0 += v * v;
        }
        for (int i = 0; i < n_block1; ++i)
        {
            const double v = vec_val(n_block0 + i, k);
            e1 += v * v;
        }

        const double kval = std::sqrt(lam);
        if (e0 >= e1)
            k_block0.push_back(kval);
        else
            k_block1.push_back(kval);
    }

    std::sort(k_block0.begin(), k_block0.end());
    std::sort(k_block1.begin(), k_block1.end());
}

// Impressao compacta de espectro usada nos casos circulares/coaxiais.
/******************************************************************************/
/* FUNCAO: print_first_modes                                                  */
/* DESCRICAO: Gera saida didatica dos resultados para inspecao no terminal ou */
/* em arquivo.                                                                */
/* ENTRADA: title: const char *; kvals: const std::vector<double> &; nprint:  */
/* int.                                                                       */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_first_modes(
    const char *title,
    const std::vector<double> &kvals,
    int nprint = 8)
{
    std::cout << title << "\n";
    for (int i = 0; i < (int)kvals.size() && i < nprint; ++i)
    {
        std::cout << " " << (i + 1) << "  " << kvals[(size_t)i] << "\n";
    }
}
