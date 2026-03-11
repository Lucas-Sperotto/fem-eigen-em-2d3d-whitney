/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge/edge_assembly.cpp                                        */
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

#include "edge_assembly.hpp"
#include "edge_basis.hpp"
#include <array>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{
constexpr std::array<std::array<double, 3>, 3> kTriQuadP2 = {{
    {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
    {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
}};

/******************************************************************************/
/* FUNCAO: assemble_edge_system_with_tri_material                             */
/* DESCRICAO: Monta as matrizes globais de aresta (rigidez e massa) somando   */
/* contribuicoes por triangulo. Implementa Sx=lambdaTx da formulacao          */
/* transversal (Secao 2.2.1).                                                 */
/* ENTRADA: mesh: const Mesh2D &; edge_dofs: EdgeDofs; eps_r_tri: const       */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &.              */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
EdgeSystem assemble_edge_system_with_tri_material(
    const Mesh2D &mesh,
    EdgeDofs edge_dofs,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri)
{
    EdgeSystem sys;
    sys.ed = std::move(edge_dofs);

    const int ndof = sys.ed.ndof;
    sys.S = DenseMat(ndof);
    sys.T = DenseMat(ndof);

    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const Tri &t = mesh.tris[tid];
        TriGeomEdge tg = tri_geom_edge(mesh, t);
        const double eps_r = eps_r_tri[tid];
        const double mu_r = mu_r_tri[tid];

        double Sel[3][3] = {{0}}, Tel[3][3] = {{0}};

        // Bloco local de rigidez (curl-curl):
        //   Sel_ij = int_T (1/mu_r) curl(W_i) curl(W_j) dA
        // Para triangulo linear, curl(W_i) e constante no elemento.
        double curl[3];
        for (int m = 0; m < 3; ++m)
            curl[m] = whitney_curl_local(m, tg);

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                Sel[i][j] = (tg.g.A / mu_r) * (curl[i] * curl[j]);
            }
        }

        // Bloco local de massa vetorial:
        //   Tel_ij = int_T eps_r (W_i . W_j) dA
        // Integrando com quadratura exata para polinomios de grau 2.
        for (const auto &lam : kTriQuadP2)
        {
            double w = tg.g.A / 3.0;
            Vec2 W[3];
            for (int m = 0; m < 3; ++m)
                W[m] = whitney_W_local(m, tg, lam);

            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    Tel[i][j] += eps_r * w * (W[i].x * W[j].x + W[i].y * W[j].y);
                }
            }
        }

        // Espalhamento local->global com correcao de orientacao:
        // o sinal local da aresta precisa ser alinhado com a orientacao global.
        const TriEdges &te = sys.ed.tri_edges[tid];

        for (int li = 0; li < 3; ++li)
        {
            int eid_i = te.e[li];
            int si = te.sgn[li];
            int I = sys.ed.edge_to_dof[eid_i];
            if (I < 0)
                continue;

            for (int lj = 0; lj < 3; ++lj)
            {
                int eid_j = te.e[lj];
                int sj = te.sgn[lj];
                int J = sys.ed.edge_to_dof[eid_j];
                if (J < 0)
                    continue;

                sys.S(I, J) += (double)(si * sj) * Sel[li][lj];
                sys.T(I, J) += (double)(si * sj) * Tel[li][lj];
            }
        }
    }

    return sys;
}

/******************************************************************************/
/* FUNCAO: make_uniform_tri_data                                              */
/* DESCRICAO: Gera e retorna dados geometricos/discretizacao usados nas       */
/* simulacoes.                                                                */
/* ENTRADA: ntri: int; value: double.                                         */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<double> make_uniform_tri_data(int ntri, double value)
{
    return std::vector<double>((size_t)ntri, value);
}
} // namespace

/******************************************************************************/
/* FUNCAO: build_helm10_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial transversal de aresta da secao 2.2.1.  */
/* Implementa Sx=lambdaTx da formulacao transversal (Secao 2.2.1).            */
/* ENTRADA: mesh: const Mesh2D &; eps_r: double; mu_r: double.                */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
EdgeSystem build_helm10_edge_system(const Mesh2D &mesh, double eps_r, double mu_r)
{
    auto eps_r_tri = make_uniform_tri_data((int)mesh.tris.size(), eps_r);
    auto mu_r_tri = make_uniform_tri_data((int)mesh.tris.size(), mu_r);
    return assemble_edge_system_with_tri_material(mesh, build_edge_dofs(mesh), eps_r_tri, mu_r_tri);
}

/******************************************************************************/
/* FUNCAO: build_helm10_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial transversal de aresta da secao 2.2.1.  */
/* Implementa Sx=lambdaTx da formulacao transversal (Secao 2.2.1).            */
/* ENTRADA: mesh: const Mesh2D &; bc: EdgeBC; eps_r: double; mu_r: double.    */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
EdgeSystem build_helm10_edge_system(const Mesh2D &mesh, EdgeBC bc, double eps_r, double mu_r)
{
    auto eps_r_tri = make_uniform_tri_data((int)mesh.tris.size(), eps_r);
    auto mu_r_tri = make_uniform_tri_data((int)mesh.tris.size(), mu_r);
    return assemble_edge_system_with_tri_material(mesh, build_edge_dofs(mesh, bc), eps_r_tri, mu_r_tri);
}

/******************************************************************************/
/* FUNCAO: build_helm10_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial transversal de aresta da secao 2.2.1.  */
/* Implementa Sx=lambdaTx da formulacao transversal (Secao 2.2.1).            */
/* ENTRADA: mesh: const Mesh2D &; bc: EdgeBC; eps_r_tri: const                */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &.              */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
EdgeSystem build_helm10_edge_system(
    const Mesh2D &mesh,
    EdgeBC bc,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri)
{
    if ((int)eps_r_tri.size() != (int)mesh.tris.size())
        throw std::runtime_error("eps_r_tri.size() != mesh.tris.size()");
    if ((int)mu_r_tri.size() != (int)mesh.tris.size())
        throw std::runtime_error("mu_r_tri.size() != mesh.tris.size()");

    return assemble_edge_system_with_tri_material(mesh, build_edge_dofs(mesh, bc), eps_r_tri, mu_r_tri);
}
