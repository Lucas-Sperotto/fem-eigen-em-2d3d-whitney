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
#include "explicit/tri2d_edge_explicit.hpp"
#include "meshfree/efgmi_2d.hpp"
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

EdgeSystem assemble_edge_system_efgmi(
    const Mesh2D &mesh,
    EdgeDofs edge_dofs,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri);

/******************************************************************************/
/* FUNCAO: assemble_edge_system_with_tri_material                             */
/* DESCRICAO: Monta as matrizes globais de aresta (rigidez e massa) somando   */
/* contribuicoes por triangulo. Implementa Sx=lambdaTx da formulacao          */
/* transversal (Secao 2.2.1). Apos a assembleia local->global, este e o       */
/* sistema matricial da Eq. (65), no papel exercido pelo HELMVEC do apendice. */
/* ENTRADA: mesh: const Mesh2D &; edge_dofs: EdgeDofs; eps_r_tri: const       */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
EdgeSystem assemble_edge_system_with_tri_material(
    const Mesh2D &mesh,
    EdgeDofs edge_dofs,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend)
{
    if (backend == ElementAssemblyBackend::EfgmiConsistent)
        return assemble_edge_system_efgmi(mesh, std::move(edge_dofs), eps_r_tri, mu_r_tri);

    EdgeSystem sys;
    sys.ed = std::move(edge_dofs);
    sys.backend = backend;

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
        if (backend == ElementAssemblyBackend::ClosedForm)
        {
            explicit_tri2d::tri2d_edge_closed_form_eq_66_67(tg, 1.0 / mu_r, eps_r, Sel, Tel);
        }
        else
        {
            // Backend por quadratura de Gauss, inlined para clareza.
            double curl[3];
            for (int m = 0; m < 3; ++m)
                curl[m] = whitney_curl_local(m, tg);

            for (const auto &lam : kTriQuadP2)
            {
                const double w = tg.g.A / 3.0;
                Vec2 W[3];
                for (int m = 0; m < 3; ++m)
                    W[m] = whitney_W_local(m, tg, lam);

                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        Sel[i][j] += (w / mu_r) * (curl[i] * curl[j]);
                        Tel[i][j] += eps_r * w * (W[i].x * W[j].x + W[i].y * W[j].y);
                    }
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

EdgeSystem assemble_edge_system_efgmi(
    const Mesh2D &mesh,
    EdgeDofs edge_dofs,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri)
{
    EdgeSystem sys;
    sys.ed = std::move(edge_dofs);
    sys.backend = ElementAssemblyBackend::EfgmiConsistent;

    const int ndof = sys.ed.ndof;
    sys.S = DenseMat(ndof);
    sys.T = DenseMat(ndof);

    const auto ctx = efgmi2d::make_edge_context(mesh, sys.ed);
    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const Tri &tri = mesh.tris[(size_t)tid];
        const Node2D &n0 = mesh.nodes[(size_t)tri.v[0]];
        const Node2D &n1 = mesh.nodes[(size_t)tri.v[1]];
        const Node2D &n2 = mesh.nodes[(size_t)tri.v[2]];
        const TriEdges &te = sys.ed.tri_edges[(size_t)tid];
        const TriGeom geom = tri_geom(mesh, tri);
        const double eps_r = eps_r_tri[(size_t)tid];
        const double mu_r = mu_r_tri[(size_t)tid];

        for (const auto &qp : efgmi2d::kTriQuadP5)
        {
            const double x =
                qp.lambda[0] * n0.x +
                qp.lambda[1] * n1.x +
                qp.lambda[2] * n2.x;
            const double y =
                qp.lambda[0] * n0.y +
                qp.lambda[1] * n1.y +
                qp.lambda[2] * n2.y;
            const double w = geom.A * qp.weight;
            const auto edge_basis = efgmi2d::evaluate_triangle_edge_basis(ctx, tid, x, y);

            for (int a = 0; a < 3; ++a)
            {
                const int I = sys.ed.edge_to_dof[(size_t)te.e[a]];
                if (I < 0)
                    continue;

                for (int b = 0; b < 3; ++b)
                {
                    const int J = sys.ed.edge_to_dof[(size_t)te.e[b]];
                    if (J < 0)
                        continue;

                    const double sign = static_cast<double>(te.sgn[a] * te.sgn[b]);
                    sys.S(I, J) += sign * (w / mu_r) *
                                   edge_basis.curl[(size_t)a] * edge_basis.curl[(size_t)b];
                    sys.T(I, J) += sign * eps_r * w *
                                   (edge_basis.vx[(size_t)a] * edge_basis.vx[(size_t)b] +
                                    edge_basis.vy[(size_t)a] * edge_basis.vy[(size_t)b]);
                }
            }
        }
    }

    return sys;
}
} // namespace

/******************************************************************************/
/* FUNCAO: build_helm10_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial transversal de aresta da secao 2.2.1.  */
/* Implementa Sx=lambdaTx da formulacao transversal (Secao 2.2.1).            */
/* Este wrapper continua representando a Eq. (65), apenas com meio homogeneo. */
/* ENTRADA: mesh: const Mesh2D &; eps_r: double; mu_r: double; backend:       */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
EdgeSystem build_helm10_edge_system(
    const Mesh2D &mesh,
    double eps_r,
    double mu_r,
    ElementAssemblyBackend backend)
{
    auto eps_r_tri = make_uniform_tri_data((int)mesh.tris.size(), eps_r);
    auto mu_r_tri = make_uniform_tri_data((int)mesh.tris.size(), mu_r);
    return assemble_edge_system_with_tri_material(mesh, build_edge_dofs(mesh), eps_r_tri, mu_r_tri, backend);
}

/******************************************************************************/
/* FUNCAO: build_helm10_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial transversal de aresta da secao 2.2.1.  */
/* Implementa Sx=lambdaTx da formulacao transversal (Secao 2.2.1).            */
/* Este wrapper continua representando a Eq. (65), agora com BC explicita.    */
/* ENTRADA: mesh: const Mesh2D &; bc: EdgeBC; eps_r: double; mu_r: double;    */
/* backend: ElementAssemblyBackend.                                           */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
EdgeSystem build_helm10_edge_system(
    const Mesh2D &mesh,
    EdgeBC bc,
    double eps_r,
    double mu_r,
    ElementAssemblyBackend backend)
{
    auto eps_r_tri = make_uniform_tri_data((int)mesh.tris.size(), eps_r);
    auto mu_r_tri = make_uniform_tri_data((int)mesh.tris.size(), mu_r);
    return assemble_edge_system_with_tri_material(mesh, build_edge_dofs(mesh, bc), eps_r_tri, mu_r_tri, backend);
}

/******************************************************************************/
/* FUNCAO: build_helm10_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial transversal de aresta da secao 2.2.1.  */
/* Implementa Sx=lambdaTx da formulacao transversal (Secao 2.2.1).            */
/* Esta variante heterogenea tambem desemboca no sistema global da Eq. (65).  */
/* ENTRADA: mesh: const Mesh2D &; bc: EdgeBC; eps_r_tri: const                */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &; backend:     */
/* ElementAssemblyBackend.                                                    */
/* SAIDA: EdgeSystem.                                                         */
/******************************************************************************/
EdgeSystem build_helm10_edge_system(
    const Mesh2D &mesh,
    EdgeBC bc,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend)
{
    if ((int)eps_r_tri.size() != (int)mesh.tris.size())
        throw std::runtime_error("eps_r_tri.size() != mesh.tris.size()");
    if ((int)mu_r_tri.size() != (int)mesh.tris.size())
        throw std::runtime_error("mu_r_tri.size() != mesh.tris.size()");

    return assemble_edge_system_with_tri_material(mesh, build_edge_dofs(mesh, bc), eps_r_tri, mu_r_tri, backend);
}
