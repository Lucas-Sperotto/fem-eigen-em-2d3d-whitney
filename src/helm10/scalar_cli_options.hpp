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
    bool ar_m_was_provided = false;
    double ar_m = 1.0;
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
    bool frequency_was_provided = false;
    double frequency_hz = 0.0;
    double eps_r = 1.0;
    double mu_r = 1.0;
    std::vector<std::string> positionals;
};

/******************************************************************************/
/* FUNCAO: parse_positive_cli_real                                            */
/* DESCRICAO: Converte um argumento textual em numero real positivo para uso  */
/* nos parametros fisicos opcionais do HELM10.                                */
/* ENTRADA: text: const std::string &; name: const std::string &.             */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double parse_positive_cli_real(
    const std::string &text,
    const std::string &name)
{
    const double value = std::stod(text);
    if (value <= 0.0)
    {
        throw std::runtime_error(name + " deve ser > 0");
    }
    return value;
}

/******************************************************************************/
/* FUNCAO: parse_positive_cli_int                                             */
/* DESCRICAO: Converte um argumento textual em inteiro positivo para uso nos  */
/* parametros discretos dos executaveis HELM10.                               */
/* ENTRADA: text: const std::string &; name: const std::string &.             */
/* SAIDA: int.                                                                */
/******************************************************************************/
inline int parse_positive_cli_int(
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
/* FUNCAO: parse_nonnegative_cli_int                                          */
/* DESCRICAO: Converte um argumento textual em inteiro nao negativo para      */
/* controlar quantidades opcionais de modos exportados.                       */
/* ENTRADA: text: const std::string &; name: const std::string &.             */
/* SAIDA: int.                                                                */
/******************************************************************************/
inline int parse_nonnegative_cli_int(
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

        if (arg == "--freq-hz")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --freq-hz");
            }
            opts.frequency_was_provided = true;
            opts.frequency_hz = parse_positive_cli_real(argv[++i], "freq_hz");
            continue;
        }
        if (arg.rfind("--freq-hz=", 0) == 0)
        {
            opts.frequency_was_provided = true;
            opts.frequency_hz = parse_positive_cli_real(
                arg.substr(std::string("--freq-hz=").size()),
                "freq_hz");
            continue;
        }

        if (arg == "--ar-m")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --ar-m");
            }
            opts.ar_m_was_provided = true;
            opts.ar_m = parse_positive_cli_real(argv[++i], "ar_m");
            continue;
        }
        if (arg.rfind("--ar-m=", 0) == 0)
        {
            opts.ar_m_was_provided = true;
            opts.ar_m = parse_positive_cli_real(
                arg.substr(std::string("--ar-m=").size()),
                "ar_m");
            continue;
        }

        if (arg == "--nx")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --nx");
            }
            opts.nx_was_provided = true;
            opts.nx = parse_positive_cli_int(argv[++i], "nx");
            continue;
        }
        if (arg.rfind("--nx=", 0) == 0)
        {
            opts.nx_was_provided = true;
            opts.nx = parse_positive_cli_int(
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
            opts.ny = parse_positive_cli_int(argv[++i], "ny");
            continue;
        }
        if (arg.rfind("--ny=", 0) == 0)
        {
            opts.ny_was_provided = true;
            opts.ny = parse_positive_cli_int(
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
            opts.nr = parse_positive_cli_int(argv[++i], "nr");
            continue;
        }
        if (arg.rfind("--nr=", 0) == 0)
        {
            opts.nr_was_provided = true;
            opts.nr = parse_positive_cli_int(
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
            opts.nt = parse_positive_cli_int(argv[++i], "nt");
            continue;
        }
        if (arg.rfind("--nt=", 0) == 0)
        {
            opts.nt_was_provided = true;
            opts.nt = parse_positive_cli_int(
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
            opts.nmodos = parse_nonnegative_cli_int(argv[++i], "nmodos");
            continue;
        }
        if (arg.rfind("--nmodos=", 0) == 0)
        {
            opts.nmodos_was_provided = true;
            opts.nmodos = parse_nonnegative_cli_int(
                arg.substr(std::string("--nmodos=").size()),
                "nmodos");
            continue;
        }

        if (arg == "--eps-r")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --eps-r");
            }
            opts.eps_r = parse_positive_cli_real(argv[++i], "eps_r");
            continue;
        }
        if (arg.rfind("--eps-r=", 0) == 0)
        {
            opts.eps_r = parse_positive_cli_real(
                arg.substr(std::string("--eps-r=").size()),
                "eps_r");
            continue;
        }

        if (arg == "--mu-r")
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("faltou valor apos --mu-r");
            }
            opts.mu_r = parse_positive_cli_real(argv[++i], "mu_r");
            continue;
        }
        if (arg.rfind("--mu-r=", 0) == 0)
        {
            opts.mu_r = parse_positive_cli_real(
                arg.substr(std::string("--mu-r=").size()),
                "mu_r");
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
