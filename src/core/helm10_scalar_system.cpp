/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/helm10_scalar_system.cpp                                 */
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

#include "helm10_scalar_system.hpp"
#include "fem_scalar.hpp"
#include <algorithm>
#include <stdexcept>

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
            map[i] = -1;
        else
            map[i] = cnt++;
    }
    return map;
}

/******************************************************************************/
/* FUNCAO: build_helm10_scalar_system                                         */
/* DESCRICAO: Monta o sistema escalar generalizado da secao 2.1 com materiais */
/* e BCs informados. Implementa a formulacao escalar da Secao 2.1.            */
/* ENTRADA: mesh: const Mesh2D &; bc: ScalarBC.                               */
/* SAIDA: ScalarSystem.                                                       */
/******************************************************************************/
ScalarSystem build_helm10_scalar_system(const Mesh2D &mesh, ScalarBC bc)
{
    // Sistema escalar da Secao 2.1:
    //   S u = lambda T u
    // com S associado a grad-grad e T a massa consistente.
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

    return {S, T, ndof, map};
}


/******************************************************************************/
/* FUNCAO: build_helm10_scalar_system                                         */
/* DESCRICAO: Monta o sistema escalar generalizado da secao 2.1 com materiais */
/* e BCs informados. Implementa a formulacao escalar da Secao 2.1.            */
/* ENTRADA: mesh: const Mesh2D &; bc: ScalarBC; eps_r_tri: const              */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &.              */
/* SAIDA: ScalarSystem.                                                       */
/******************************************************************************/
ScalarSystem build_helm10_scalar_system(
    const Mesh2D &mesh,
    ScalarBC bc,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri
){
    if ((int)eps_r_tri.size() != (int)mesh.tris.size())
        throw std::runtime_error("eps_r_tri.size() != mesh.tris.size()");
    if ((int)mu_r_tri.size() != (int)mesh.tris.size())
        throw std::runtime_error("mu_r_tri.size() != mesh.tris.size()");

    auto map = build_dof_map(mesh, bc);
    int ndof = 0;
    for (int v : map)
        if (v >= 0)
            ndof = std::max(ndof, v + 1);
    if (ndof <= 0)
        throw std::runtime_error("ndof invalido.");

    DenseMat S(ndof), T(ndof);

    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const auto &tri = mesh.tris[tid];
        TriGeom g = tri_geom(mesh, tri);

        const double eps_r = eps_r_tri[tid];
        const double mu_r  = mu_r_tri[tid];

        double Se[3][3], Te[3][3];
        element_mats_scalar(g, Se, Te);

        // Aplica ponderacoes inhomogeneas conforme formulacao do artigo:
        // Sz ~ (1/mu) * grad-grad
        // Tz ~ eps * mass
        for (int a = 0; a < 3; ++a)
        {
            for (int b2 = 0; b2 < 3; ++b2)
            {
                Se[a][b2] *= (1.0 / mu_r);
                Te[a][b2] *= eps_r;
            }
        }

        for (int a = 0; a < 3; a++)
        {
            int ia = map[tri.v[a]];
            if (ia < 0)
                continue;
            for (int b2 = 0; b2 < 3; b2++)
            {
                int ib = map[tri.v[b2]];
                if (ib < 0)
                    continue;
                S(ia, ib) += Se[a][b2];
                T(ia, ib) += Te[a][b2];
            }
        }
    }

    return {S, T, ndof, map};
}
