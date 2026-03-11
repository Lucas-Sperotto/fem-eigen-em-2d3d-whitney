/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge/edge_dofs.hpp                                            */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo 2D de elementos de aresta (DOFs, base de Whitney e       */
/* montagem).                                                                 */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.1.        */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#pragma once
#include "core/mesh2d.hpp"
#include <vector>
#include <array>
#include <unordered_map>
#include <cstdint>

struct Edge
{
    int n0 = -1, n1 = -1;     // orientação global: n0 < n1
    int triL = -1, triR = -1; // triângulos incidentes (opcional)
    bool is_boundary = false; // PEC
};

enum class EdgeBC
{
    TE_PEC_TangentialZero, // elimina arestas de contorno
    TM_PEC_NormalZero      // mantém arestas de contorno (DOF tangencial livre)
};

struct TriEdges
{
    std::array<int, 3> e;   // ids globais das 3 arestas do triângulo
    std::array<int, 3> sgn; // +1 se orientação local coincide com global, -1 caso contrário
};

struct EdgeDofs
{
    std::vector<Edge> edges;
    std::vector<TriEdges> tri_edges;
    std::vector<int> edge_to_dof; // -1 se contorno (PEC), senao 0..ndof-1
    int ndof = 0;
};

/******************************************************************************/
/* FUNCAO: build_edge_dofs                                                    */
/* DESCRICAO: Monta e retorna a estrutura principal do modulo, combinando     */
/* malha, materiais e condicoes de contorno.                                  */
/* ENTRADA: mesh: const Mesh2D &.                                             */
/* SAIDA: EdgeDofs.                                                           */
/******************************************************************************/
EdgeDofs build_edge_dofs(const Mesh2D &mesh);

EdgeDofs build_edge_dofs(const Mesh2D &mesh, EdgeBC bc);
