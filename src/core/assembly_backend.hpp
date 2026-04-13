/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/assembly_backend.hpp                                     */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Seletores de backend de montagem elementar.                     */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secoes 2 e 3.       */
/*****************************************************************************/
/* Observacao: Este cabecalho permite escolher entre formas fechadas          */
/* (closed-form) e integracao por quadratura de Gauss/cubatura.              */
/*****************************************************************************/

#pragma once

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

enum class ElementAssemblyBackend
{
    ClosedForm,
    GaussianQuadrature,
    EfgmiConsistent
};

/******************************************************************************/
/* FUNCAO: element_assembly_backend_name                                      */
/* DESCRICAO: Retorna um nome curto e estavel para o backend de montagem      */
/* selecionado, util para logs, CSVs e mensagens de linha de comando.         */
/* ENTRADA: backend: ElementAssemblyBackend.                                  */
/* SAIDA: const char *.                                                       */
/******************************************************************************/
inline const char *element_assembly_backend_name(ElementAssemblyBackend backend)
{
    switch (backend)
    {
        case ElementAssemblyBackend::ClosedForm:
            return "closed-form";
        case ElementAssemblyBackend::GaussianQuadrature:
            return "gauss";
        case ElementAssemblyBackend::EfgmiConsistent:
            return "efgmi";
    }
    return "desconhecido";
}

/******************************************************************************/
/* FUNCAO: normalize_backend_token                                            */
/* DESCRICAO: Normaliza uma string de linha de comando removendo separadores e*/
/* convertendo para minusculas, o que simplifica o parse de aliases.          */
/* ENTRADA: token: std::string.                                               */
/* SAIDA: std::string.                                                        */
/******************************************************************************/
inline std::string normalize_backend_token(std::string token)
{
    std::transform(
        token.begin(),
        token.end(),
        token.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });

    std::string out;
    out.reserve(token.size());
    for (unsigned char ch : token)
    {
        if (std::isalnum(ch))
        {
            out.push_back(static_cast<char>(ch));
        }
    }
    return out;
}

/******************************************************************************/
/* FUNCAO: parse_element_assembly_backend                                     */
/* DESCRICAO: Interpreta um token textual e converte para o backend de        */
/* montagem correspondente. Aceita aliases em ingles para facilitar scripts.  */
/* ENTRADA: token: const std::string &.                                       */
/* SAIDA: ElementAssemblyBackend.                                             */
/******************************************************************************/
inline ElementAssemblyBackend parse_element_assembly_backend(const std::string &token)
{
    const std::string normalized = normalize_backend_token(token);

    if (normalized == "closedform" ||
        normalized == "closed" ||
        normalized == "explicit" ||
        normalized == "cf")
    {
        return ElementAssemblyBackend::ClosedForm;
    }

    if (normalized == "gauss" ||
        normalized == "gaussian" ||
        normalized == "gaussianquadrature" ||
        normalized == "quadrature" ||
        normalized == "gq")
    {
        return ElementAssemblyBackend::GaussianQuadrature;
    }

    if (normalized == "efgmi" ||
        normalized == "meshfree" ||
        normalized == "consistency" ||
        normalized == "consistent" ||
        normalized == "efgmiconsistent")
    {
        return ElementAssemblyBackend::EfgmiConsistent;
    }

    throw std::runtime_error(
        "backend invalido: use 'closed-form', 'gauss' ou 'efgmi'.");
}
