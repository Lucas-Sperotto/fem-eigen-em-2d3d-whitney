/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec1/mixed_coax_edge_match_wrapper.hpp                    */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Wrapper local para evitar colisao de simbolos entre os matchers */
/* escalares e de aresta do caso coaxial.                                     */
/*****************************************************************************/

#pragma once

#define Jm helmvec1_edge_coax_Jm
#define Ym helmvec1_edge_coax_Ym
#define Jm_prime helmvec1_edge_coax_Jm_prime
#define Ym_prime helmvec1_edge_coax_Ym_prime
#define mass_inner helmvec1_edge_coax_mass_inner
#define mass_norm helmvec1_edge_coax_mass_norm
#define mass_correlation_abs helmvec1_edge_coax_mass_correlation_abs
#define find_root_bisection helmvec1_edge_coax_find_root_bisection
#define det_TM helmvec1_edge_coax_det_TM
#define det_TE helmvec1_edge_coax_det_TE
#define coax_roots helmvec1_edge_coax_coax_roots
#include "edge/mode_match_coax_edge.hpp"
#undef coax_roots
#undef det_TE
#undef det_TM
#undef find_root_bisection
#undef mass_correlation_abs
#undef mass_norm
#undef mass_inner
#undef Ym_prime
#undef Jm_prime
#undef Ym
#undef Jm
