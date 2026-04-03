/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/scalar_cli_options.hpp                                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios de linha de comando para os executaveis escalares   */
/* da Secao 2.1.                                                              */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.1.          */
/*****************************************************************************/
/* Observacao: O objetivo aqui e padronizar a escolha de backend sem quebrar  */
/* os argumentos posicionais legados usados pelos scripts atuais.             */
/*****************************************************************************/

#pragma once

#include "core/assembly_backend.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace helm10
{

struct ScalarCliOptions
{
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm;
    bool debug_local_blocks = false;
    bool debug_candidates = false;
    std::vector<std::string> positionals;
};

/******************************************************************************/
/* FUNCAO: parse_scalar_cli_options                                           */
/* DESCRICAO: Separa argumentos posicionais legados de opcoes nomeadas e      */
/* interpreta backend e flags de depuracao para os executaveis escalares da   */
/* Secao 2.1.                                                                 */
/* ENTRADA: argc: int; argv: char **.                                         */
/* SAIDA: ScalarCliOptions.                                                   */
/******************************************************************************/
inline ScalarCliOptions parse_scalar_cli_options(int argc, char **argv)
{
    ScalarCliOptions opts;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--backend")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --backend");
            }
            opts.backend = parse_element_assembly_backend(argv[++i]);
            continue;
        }

        const std::string prefix = "--backend=";
        if (arg.rfind(prefix, 0) == 0)
        {
            opts.backend = parse_element_assembly_backend(arg.substr(prefix.size()));
            continue;
        }

        if (arg == "--debug" || arg == "--debug-all")
        {
            opts.debug_local_blocks = true;
            opts.debug_candidates = true;
            continue;
        }
        if (arg == "--debug-local-blocks")
        {
            opts.debug_local_blocks = true;
            continue;
        }
        if (arg == "--debug-candidates")
        {
            opts.debug_candidates = true;
            continue;
        }

        if (!arg.empty() && arg[0] == '-')
        {
            throw std::runtime_error(
                "opcao desconhecida na formulacao escalar 2D: " + arg);
        }

        opts.positionals.push_back(arg);
    }

    return opts;
}

} // namespace helm10
