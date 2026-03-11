/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/helm10_scalar_system.hpp                                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo escalar 2D (malha, montagem e identificacao modal).      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.1.          */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "mesh2d.hpp"
#include "dense.hpp"
#include <vector>

enum class ScalarBC
{
    TE_Neumann,
    TM_Dirichlet
};

struct ScalarSystem
{
    DenseMat S;
    DenseMat T;
    int ndof = 0;
    std::vector<int> dof_map; // node -> dof (>=0) ou -1 eliminado
};

/******************************************************************************/
/* FUNCAO: build_helm10_scalar_system                                         */
/* DESCRICAO: Monta o sistema escalar generalizado da secao 2.1 com materiais */
/* e BCs informados. Implementa a formulacao escalar da Secao 2.1.            */
/* ENTRADA: mesh: const Mesh2D &; bc: ScalarBC.                               */
/* SAIDA: ScalarSystem.                                                       */
/******************************************************************************/
ScalarSystem build_helm10_scalar_system(const Mesh2D &mesh, ScalarBC bc);

// Variante nao homogenea (parametros por triangulo).
ScalarSystem build_helm10_scalar_system(
    const Mesh2D &mesh,
    ScalarBC bc,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri
);
