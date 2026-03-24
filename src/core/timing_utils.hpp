/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/timing_utils.hpp                                         */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Utilitarios simples de medicao de tempo para benchmark local.   */
/*****************************************************************************/
/* Observacao: A intencao e padronizar a comparacao entre backends            */
/* closed-form e quadratura de Gauss sem depender de cronometro externo.      */
/*****************************************************************************/

#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

namespace timing
{

struct Breakdown
{
  double assembly_ms = 0.0;
  double solve_ms = 0.0;
  double post_ms = 0.0;
  double total_ms = 0.0;
};

class Stopwatch
{
public:
  /****************************************************************************/
  /* FUNCAO: Stopwatch                                                        */
  /* DESCRICAO: Inicializa o cronometro no instante atual.                    */
  /* ENTRADA: sem parametros.                                                 */
  /* SAIDA: objeto Stopwatch.                                                 */
  /****************************************************************************/
  Stopwatch()
      : t0_(clock::now())
  {
  }

  /****************************************************************************/
  /* FUNCAO: reset                                                            */
  /* DESCRICAO: Reinicia o cronometro a partir do instante atual.             */
  /* ENTRADA: sem parametros.                                                 */
  /* SAIDA: sem retorno explicito (void).                                     */
  /****************************************************************************/
  void reset()
  {
    t0_ = clock::now();
  }

  /****************************************************************************/
  /* FUNCAO: elapsed_ms                                                       */
  /* DESCRICAO: Retorna o tempo decorrido, em milissegundos, desde o ultimo   */
  /* reset (ou desde a construcao do objeto).                                 */
  /* ENTRADA: sem parametros.                                                 */
  /* SAIDA: double.                                                           */
  /****************************************************************************/
  double elapsed_ms() const
  {
    return std::chrono::duration<double, std::milli>(clock::now() - t0_).count();
  }

private:
  using clock = std::chrono::steady_clock;
  clock::time_point t0_;
};

/******************************************************************************/
/* FUNCAO: print_breakdown                                                    */
/* DESCRICAO: Imprime uma linha padronizada de benchmark para facilitar       */
/* parsing por scripts e comparacao entre executaveis/backends.               */
/* ENTRADA: label: const std::string &; b: const Breakdown &.                 */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
inline void print_breakdown(const std::string &label, const Breakdown &b)
{
  auto old_flags = std::cout.flags();
  auto old_prec = std::cout.precision();
  std::cout << std::fixed << std::setprecision(3)
            << "[timing] label=" << label
            << " assembly_ms=" << b.assembly_ms
            << " solve_ms=" << b.solve_ms
            << " post_ms=" << b.post_ms
            << " total_ms=" << b.total_ms << "\n";
  std::cout.flags(old_flags);
  std::cout.precision(old_prec);
}

} // namespace timing
