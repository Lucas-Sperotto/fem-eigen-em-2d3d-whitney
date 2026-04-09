/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec2/coupled_cli_options.hpp                              */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios de linha de comando para os executaveis acoplados   */
/* das secoes 2.2.3 e 2.2.4.                                                  */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secoes 2.2.3-2.2.4. */
/*****************************************************************************/

#pragma once

#include "core/assembly_backend.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace helmvec2
{

struct CoupledCliOptions
{
    ElementAssemblyBackend backend = ElementAssemblyBackend::ClosedForm;
    bool debug_local_blocks = false;
    bool debug_candidates = false;
    bool beta_was_provided = false;
    double beta = 0.0;
    bool d_over_a_preview_was_provided = false;
    double d_over_a_preview = 0.0;
    bool nx_was_provided = false;
    int nx = 0;
    bool ny_was_provided = false;
    int ny = 0;
    std::vector<std::string> positionals;
};

inline double parse_cli_real(
    const std::string &text,
    const std::string &name)
{
    try
    {
        size_t idx = 0;
        const double value = std::stod(text, &idx);
        if (idx != text.size())
            throw std::runtime_error("");
        return value;
    }
    catch (const std::exception &)
    {
        throw std::runtime_error(name + " invalido: " + text);
    }
}

inline int parse_positive_coupled_cli_int(
    const std::string &text,
    const std::string &name)
{
    const int value = std::stoi(text);
    if (value <= 0)
    {
        throw std::runtime_error(name + " deve ser > 0");
    }
    return value;
}

/******************************************************************************/
/* FUNCAO: parse_coupled_cli_options                                          */
/* DESCRICAO: Separa argumentos posicionais legados de opcoes nomeadas e      */
/* interpreta as chaves de backend e depuracao para os drivers acoplados.     */
/* ENTRADA: argc: int; argv: char **.                                         */
/* SAIDA: CoupledCliOptions.                                                  */
/******************************************************************************/
inline CoupledCliOptions parse_coupled_cli_options(int argc, char **argv)
{
    CoupledCliOptions opts;

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

        if (arg == "--beta")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --beta");
            }
            opts.beta_was_provided = true;
            opts.beta = parse_cli_real(argv[++i], "beta");
            continue;
        }
        if (arg.rfind("--beta=", 0) == 0)
        {
            opts.beta_was_provided = true;
            opts.beta = parse_cli_real(
                arg.substr(std::string("--beta=").size()),
                "beta");
            continue;
        }

        if (arg == "--d-over-a-preview")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --d-over-a-preview");
            }
            opts.d_over_a_preview_was_provided = true;
            opts.d_over_a_preview = parse_cli_real(argv[++i], "d_over_a_preview");
            continue;
        }
        if (arg.rfind("--d-over-a-preview=", 0) == 0)
        {
            opts.d_over_a_preview_was_provided = true;
            opts.d_over_a_preview = parse_cli_real(
                arg.substr(std::string("--d-over-a-preview=").size()),
                "d_over_a_preview");
            continue;
        }

        if (arg == "--nx")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --nx");
            }
            opts.nx_was_provided = true;
            opts.nx = parse_positive_coupled_cli_int(argv[++i], "nx");
            continue;
        }
        if (arg.rfind("--nx=", 0) == 0)
        {
            opts.nx_was_provided = true;
            opts.nx = parse_positive_coupled_cli_int(
                arg.substr(std::string("--nx=").size()),
                "nx");
            continue;
        }

        if (arg == "--ny")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --ny");
            }
            opts.ny_was_provided = true;
            opts.ny = parse_positive_coupled_cli_int(argv[++i], "ny");
            continue;
        }
        if (arg.rfind("--ny=", 0) == 0)
        {
            opts.ny_was_provided = true;
            opts.ny = parse_positive_coupled_cli_int(
                arg.substr(std::string("--ny=").size()),
                "ny");
            continue;
        }

        if (!arg.empty() && arg[0] == '-')
        {
            throw std::runtime_error("opcao desconhecida na formulacao acoplada: " + arg);
        }

        opts.positionals.push_back(arg);
    }

    return opts;
}

} // namespace helmvec2
