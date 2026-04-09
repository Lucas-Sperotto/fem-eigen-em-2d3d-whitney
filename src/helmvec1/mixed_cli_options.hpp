/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_cli_options.hpp                                */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios de linha de comando para os executaveis mistos da   */
/* Secao 2.2.2.                                                               */
/*****************************************************************************/

#pragma once

#include "core/assembly_backend.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace helmvec1
{

struct MixedCliOptions
{
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm;
    bool debug_local_blocks = false;
    bool debug_candidates = false;
    bool nx_was_provided = false;
    int nx = 0;
    bool ny_was_provided = false;
    int ny = 0;
    bool nr_was_provided = false;
    int nr = 0;
    bool nt_was_provided = false;
    int nt = 0;
    std::vector<std::string> positionals;
};

inline int parse_positive_mixed_cli_int(
    const std::string &text,
    const std::string &name)
{
    try
    {
        size_t idx = 0;
        const int value = std::stoi(text, &idx);
        if (idx != text.size())
            throw std::runtime_error(name + " invalido: " + text);
        if (value <= 0)
        {
            throw std::runtime_error(name + " deve ser > 0");
        }
        return value;
    }
    catch (const std::runtime_error &)
    {
        throw;
    }
    catch (const std::exception &)
    {
        throw std::runtime_error(name + " invalido: " + text);
    }
}

inline bool mixed_cli_requests_help(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--help")
            return true;
    }
    return false;
}

/******************************************************************************/
/* FUNCAO: parse_mixed_cli_options                                            */
/* DESCRICAO: Separa argumentos posicionais e processa a chave --backend para */
/* os executaveis mistos da Secao 2.2.2, incluindo flags de depuracao.        */
/* ENTRADA: argc: int; argv: char **.                                         */
/* SAIDA: MixedCliOptions.                                                    */
/******************************************************************************/
inline MixedCliOptions parse_mixed_cli_options(int argc, char **argv)
{
    MixedCliOptions opts;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--backend")
        {
            if (i + 1 >= argc)
                throw std::runtime_error("faltou valor apos --backend");
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

        if (arg == "--nx")
        {
            if (i + 1 >= argc)
                throw std::runtime_error("faltou valor apos --nx");
            opts.nx_was_provided = true;
            opts.nx = parse_positive_mixed_cli_int(argv[++i], "nx");
            continue;
        }
        if (arg.rfind("--nx=", 0) == 0)
        {
            opts.nx_was_provided = true;
            opts.nx = parse_positive_mixed_cli_int(
                arg.substr(std::string("--nx=").size()),
                "nx");
            continue;
        }

        if (arg == "--ny")
        {
            if (i + 1 >= argc)
                throw std::runtime_error("faltou valor apos --ny");
            opts.ny_was_provided = true;
            opts.ny = parse_positive_mixed_cli_int(argv[++i], "ny");
            continue;
        }
        if (arg.rfind("--ny=", 0) == 0)
        {
            opts.ny_was_provided = true;
            opts.ny = parse_positive_mixed_cli_int(
                arg.substr(std::string("--ny=").size()),
                "ny");
            continue;
        }

        if (arg == "--nr")
        {
            if (i + 1 >= argc)
                throw std::runtime_error("faltou valor apos --nr");
            opts.nr_was_provided = true;
            opts.nr = parse_positive_mixed_cli_int(argv[++i], "nr");
            continue;
        }
        if (arg.rfind("--nr=", 0) == 0)
        {
            opts.nr_was_provided = true;
            opts.nr = parse_positive_mixed_cli_int(
                arg.substr(std::string("--nr=").size()),
                "nr");
            continue;
        }

        if (arg == "--nt")
        {
            if (i + 1 >= argc)
                throw std::runtime_error("faltou valor apos --nt");
            opts.nt_was_provided = true;
            opts.nt = parse_positive_mixed_cli_int(argv[++i], "nt");
            continue;
        }
        if (arg.rfind("--nt=", 0) == 0)
        {
            opts.nt_was_provided = true;
            opts.nt = parse_positive_mixed_cli_int(
                arg.substr(std::string("--nt=").size()),
                "nt");
            continue;
        }

        if (!arg.empty() && arg[0] == '-')
            throw std::runtime_error("opcao desconhecida na formulacao mista: " + arg);

        opts.positionals.push_back(arg);
    }

    return opts;
}

} // namespace helmvec1
