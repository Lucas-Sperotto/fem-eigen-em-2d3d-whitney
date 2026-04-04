/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_mode_utils.hpp                                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
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
#include <string>
#include <vector>

struct BlockEnergyMode
{
    int eig_index = -1;
    double lambda = 0.0;
    double kc = 0.0;
    double edge_energy = 0.0;
    double scalar_energy = 0.0;
    double dominant_energy_ratio = 0.0;
    std::string dominant_block;
};

/******************************************************************************/
/* FUNCAO: collect_modes_by_block_energy                                      */
/* DESCRICAO: Coleta modos fisicos com suas energias por bloco, preservando o */
/* indice espectral e a classificacao por bloco dominante.                    */
/* ENTRADA: res: const GenEigResult &; n_edge: int; n_scalar: int;            */
/* lambda_min: double; edge_modes: std::vector<BlockEnergyMode> &;            */
/* scalar_modes: std::vector<BlockEnergyMode> &.                              */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void collect_modes_by_block_energy(
    const GenEigResult &res,
    int n_edge,
    int n_scalar,
    double lambda_min,
    std::vector<BlockEnergyMode> &edge_modes,
    std::vector<BlockEnergyMode> &scalar_modes)
{
    edge_modes.clear();
    scalar_modes.clear();

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

        double e_edge = 0.0;
        double e_scalar = 0.0;
        for (int i = 0; i < n_edge; ++i)
        {
            const double v = vec_val(i, k);
            e_edge += v * v;
        }
        for (int i = 0; i < n_scalar; ++i)
        {
            const double v = vec_val(n_edge + i, k);
            e_scalar += v * v;
        }

        BlockEnergyMode rec;
        rec.eig_index = k;
        rec.lambda = lam;
        rec.kc = std::sqrt(lam);
        rec.edge_energy = e_edge;
        rec.scalar_energy = e_scalar;
        const double etot = e_edge + e_scalar;
        rec.dominant_energy_ratio =
            (etot > 0.0) ? std::max(e_edge, e_scalar) / etot : 0.0;
        rec.dominant_block = (e_edge >= e_scalar) ? "edge" : "scalar";

        if (rec.dominant_block == "edge")
            edge_modes.push_back(rec);
        else
            scalar_modes.push_back(rec);
    }

    auto by_kc = [](const BlockEnergyMode &a, const BlockEnergyMode &b)
    {
        if (a.kc != b.kc)
            return a.kc < b.kc;
        return a.eig_index < b.eig_index;
    };
    std::sort(edge_modes.begin(), edge_modes.end(), by_kc);
    std::sort(scalar_modes.begin(), scalar_modes.end(), by_kc);
}

/******************************************************************************/
/* FUNCAO: extract_kc_values                                                  */
/* DESCRICAO: Extrai apenas os valores kc de um conjunto de modos ja          */
/* classificados por energia dominante.                                       */
/* ENTRADA: modes: const std::vector<BlockEnergyMode> &.                      */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
inline std::vector<double> extract_kc_values(const std::vector<BlockEnergyMode> &modes)
{
    std::vector<double> out;
    out.reserve(modes.size());
    for (const BlockEnergyMode &mode : modes)
        out.push_back(mode.kc);
    return out;
}

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
    std::vector<BlockEnergyMode> block0_modes;
    std::vector<BlockEnergyMode> block1_modes;
    collect_modes_by_block_energy(
        res,
        n_block0,
        n_block1,
        lambda_min,
        block0_modes,
        block1_modes);
    k_block0 = extract_kc_values(block0_modes);
    k_block1 = extract_kc_values(block1_modes);
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
