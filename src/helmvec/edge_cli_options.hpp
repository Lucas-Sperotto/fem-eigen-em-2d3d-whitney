/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec/edge_cli_options.hpp                                  */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios de linha de comando para os executaveis vetoriais   */
/* transversais da Secao 2.2.1.                                               */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.1.        */
/*****************************************************************************/
/* Observacao: O parser preserva os argumentos posicionais atuais e adiciona  */
/* apenas a selecao do backend local por --backend.                           */
/*****************************************************************************/

#pragma once

#include "core/assembly_backend.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace helmvec
{

struct EdgeCliOptions
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
    bool nmodos_was_provided = false;
    int nmodos = 0;
    std::vector<std::string> positionals;
};

inline int parse_positive_edge_cli_int(
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

inline int parse_nonnegative_edge_cli_int(
    const std::string &text,
    const std::string &name)
{
    const int value = std::stoi(text);
    if (value < 0)
    {
        throw std::runtime_error(name + " deve ser >= 0");
    }
    return value;
}

/******************************************************************************/
/* FUNCAO: parse_edge_cli_options                                             */
/* DESCRICAO: Separa argumentos posicionais e processa a chave --backend para */
/* os executaveis de aresta da Secao 2.2.1, incluindo flags de depuracao.     */
/* ENTRADA: argc: int; argv: char **.                                         */
/* SAIDA: EdgeCliOptions.                                                     */
/******************************************************************************/
inline EdgeCliOptions parse_edge_cli_options(int argc, char **argv)
{
    EdgeCliOptions opts;

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

        if (arg == "--nx")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --nx");
            }
            opts.nx_was_provided = true;
            opts.nx = parse_positive_edge_cli_int(argv[++i], "nx");
            continue;
        }
        if (arg.rfind("--nx=", 0) == 0)
        {
            opts.nx_was_provided = true;
            opts.nx = parse_positive_edge_cli_int(
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
            opts.ny = parse_positive_edge_cli_int(argv[++i], "ny");
            continue;
        }
        if (arg.rfind("--ny=", 0) == 0)
        {
            opts.ny_was_provided = true;
            opts.ny = parse_positive_edge_cli_int(
                arg.substr(std::string("--ny=").size()),
                "ny");
            continue;
        }

        if (arg == "--nr")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --nr");
            }
            opts.nr_was_provided = true;
            opts.nr = parse_positive_edge_cli_int(argv[++i], "nr");
            continue;
        }
        if (arg.rfind("--nr=", 0) == 0)
        {
            opts.nr_was_provided = true;
            opts.nr = parse_positive_edge_cli_int(
                arg.substr(std::string("--nr=").size()),
                "nr");
            continue;
        }

        if (arg == "--nt")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --nt");
            }
            opts.nt_was_provided = true;
            opts.nt = parse_positive_edge_cli_int(argv[++i], "nt");
            continue;
        }
        if (arg.rfind("--nt=", 0) == 0)
        {
            opts.nt_was_provided = true;
            opts.nt = parse_positive_edge_cli_int(
                arg.substr(std::string("--nt=").size()),
                "nt");
            continue;
        }

        if (arg == "--nmodos")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --nmodos");
            }
            opts.nmodos_was_provided = true;
            opts.nmodos = parse_nonnegative_edge_cli_int(argv[++i], "nmodos");
            continue;
        }
        if (arg.rfind("--nmodos=", 0) == 0)
        {
            opts.nmodos_was_provided = true;
            opts.nmodos = parse_nonnegative_edge_cli_int(
                arg.substr(std::string("--nmodos=").size()),
                "nmodos");
            continue;
        }

        if (!arg.empty() && arg[0] == '-')
        {
            throw std::runtime_error(
                "opcao desconhecida na formulacao vetorial 2D: " + arg);
        }

        opts.positionals.push_back(arg);
    }

    return opts;
}

} // namespace helmvec
