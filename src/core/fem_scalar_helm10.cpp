/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/fem_scalar_helm10.cpp                                    */
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

#include "fem_scalar.hpp"
#include "lapack_eig.hpp"
#include <vector>
#include <algorithm>
#include <stdexcept>

enum class ScalarBC
{
    TE_Neumann,
    TM_Dirichlet
};

struct ScalarResult
{
    std::vector<double> lambdas;   // kc^2
    std::vector<double> modes_col; // autovetores em col-major (ndof x ndof)
    int ndof = 0;
    std::vector<int> dof_map; // node -> dof (-1 se eliminado)
};

/******************************************************************************/
/* FUNCAO: build_dof_map                                                      */
/* DESCRICAO: Gera o mapeamento de no para grau de liberdade conforme a       */
/* condicao de contorno.                                                      */
/* ENTRADA: m: const Mesh2D &; bc: ScalarBC.                                  */
/* SAIDA: std::vector<int>.                                                   */
/******************************************************************************/
static std::vector<int> build_dof_map(const Mesh2D &m, ScalarBC bc)
{
    std::vector<int> map(m.nodes.size(), -1);
    int cnt = 0;
    for (size_t i = 0; i < m.nodes.size(); i++)
    {
        bool bd = m.nodes[i].is_boundary;
        if (bc == ScalarBC::TM_Dirichlet && bd)
        {
            map[i] = -1;
        }
        else
        {
            map[i] = cnt++;
        }
    }
    return map;
}

/******************************************************************************/
/* FUNCAO: solve_helm10_scalar                                                */
/* DESCRICAO: Resolve o problema numerico associado e retorna as grandezas de */
/* interesse.                                                                 */
/* ENTRADA: mesh: const Mesh2D &; bc: ScalarBC.                               */
/* SAIDA: ScalarResult.                                                       */
/******************************************************************************/
ScalarResult solve_helm10_scalar(const Mesh2D &mesh, ScalarBC bc)
{
    auto map = build_dof_map(mesh, bc);
    int ndof = 0;
    for (int v : map)
        if (v >= 0)
            ndof = std::max(ndof, v + 1);
    if (ndof <= 0)
        throw std::runtime_error("ndof invalido.");

    DenseMat S(ndof), T(ndof);

    for (const auto &tri : mesh.tris)
    {
        TriGeom g = tri_geom(mesh, tri);

        double Se[3][3], Te[3][3];
        element_mats_scalar(g, Se, Te);

        for (int a = 0; a < 3; a++)
        {
            int ia = map[tri.v[a]];
            if (ia < 0)
                continue;
            for (int b = 0; b < 3; b++)
            {
                int ib = map[tri.v[b]];
                if (ib < 0)
                    continue;
                S(ia, ib) += Se[a][b];
                T(ia, ib) += Te[a][b];
            }
        }
    }

    auto res = generalized_eigs_sym_vec(S, T);
    return {res.w, res.Zcol, ndof, map};
}
