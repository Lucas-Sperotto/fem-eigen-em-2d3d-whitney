/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/edge3d/edge3d_assembly.cpp                                    */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Nucleo 3D de elementos de aresta tetraedricos (Whitney 1-form) e*/
/* montagem.                                                                  */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 3.1, integrais*/
/* I1..I10.                                                                   */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "edge3d_assembly.hpp"
#include "edge3d_basis.hpp"
#include <array>
#include <stdexcept>
#include <utility>

namespace
{
// Quadratura tetraedrica de 4 pontos (exata para polinomios de grau 2).
// Necessaria para os termos de massa com bases de Whitney 1-forma.
constexpr std::array<std::array<double, 4>, 4> kTetQuadP2 = {{
    {{0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105}},
    {{0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.1381966011250105}},
    {{0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105}},
    {{0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685}},
}};

/******************************************************************************/
/* FUNCAO: dot3                                                               */
/* DESCRICAO: Calcula produto interno entre vetores 3D, usado nas integrais de rigidez e massa tetraedricas. */
/* ENTRADA: a: Vec3d; b: Vec3d.                                               */
/* SAIDA: double.                                                             */
/******************************************************************************/
double dot3(Vec3d a, Vec3d b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

/******************************************************************************/
/* FUNCAO: uniform_data                                                       */
/* DESCRICAO: Cria vetor de material uniforme por tetraedro para facilitar    */
/* chamadas com meios homogeneos.                                             */
/* ENTRADA: n: int; val: double.                                              */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<double> uniform_data(int n, double val)
{
  return std::vector<double>((size_t)n, val);
}

template <typename AddFn>
/******************************************************************************/
/* FUNCAO: assemble_generic                                                   */
/* DESCRICAO: Monta contribuicoes locais 3D e acumula no sistema global       */
/* conforme politica de armazenamento. Implementa a formulacao vetorial 3D da */
/* Secao 3.1, incluindo os termos equivalentes aos coeficientes I1..I10.      */
/* ENTRADA: mesh: const Mesh3D &; ed: const EdgeDofs3D &; eps_r_tet: const    */
/* std::vector<double> &; mu_r_tet: const std::vector<double> &; add_global:  */
/* AddFn.                                                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void assemble_generic(
    const Mesh3D &mesh,
    const EdgeDofs3D &ed,
    const std::vector<double> &eps_r_tet,
    const std::vector<double> &mu_r_tet,
    AddFn add_global)
{
  for (int tid = 0; tid < (int)mesh.tets.size(); ++tid)
  {
    const Tet &t = mesh.tets[(size_t)tid];
    const TetEdges &te = ed.tet_edges[(size_t)tid];
    const TetGeomEdge tg = tet_geom_edge(mesh, t);
    const double eps_r = eps_r_tet[(size_t)tid];
    const double mu_r = mu_r_tet[(size_t)tid];

    double Sel[6][6] = {{0.0}};
    double Tel[6][6] = {{0.0}};

    // Bloco local curl-curl (Eq. 176):
    //   Sel_ij = int (1/mu_r) curl(W_i).curl(W_j) dV
    // Em tetraedros lineares, curl(W_i) e constante em cada elemento:
    // avaliacao exata por multiplicacao do volume.
    Vec3d curlW[6];
    for (int m = 0; m < 6; ++m)
      curlW[m] = whitney_curl_local_3d(m, tg);

    for (int i = 0; i < 6; ++i)
    {
      for (int j = 0; j < 6; ++j)
      {
        Sel[i][j] = (tg.V / mu_r) * dot3(curlW[i], curlW[j]);
      }
    }

    // Bloco local de massa vetorial (Eq. 177):
    //   Tel_ij = int eps_r W_i.W_j dV
    // Para base de Whitney linear, W_i.W_j e polinomio de grau 2:
    // a quadratura de 4 pontos adotada integra exatamente esse termo.
    for (const auto &lam : kTetQuadP2)
    {
      const double w = tg.V / 4.0;
      Vec3d W[6];
      for (int m = 0; m < 6; ++m)
        W[m] = whitney_W_local_3d(m, tg, lam);

      for (int i = 0; i < 6; ++i)
      {
        for (int j = 0; j < 6; ++j)
        {
          Tel[i][j] += eps_r * w * dot3(W[i], W[j]);
        }
      }
    }

    for (int li = 0; li < 6; ++li)
    {
      const int eid_i = te.e[(size_t)li];
      const int sgn_i = te.sgn[(size_t)li];
      const int I = ed.edge_to_dof[(size_t)eid_i];
      if (I < 0)
        continue;

      for (int lj = 0; lj < 6; ++lj)
      {
        const int eid_j = te.e[(size_t)lj];
        const int sgn_j = te.sgn[(size_t)lj];
        const int J = ed.edge_to_dof[(size_t)eid_j];
        if (J < 0)
          continue;

        // Correcao de orientacao da aresta:
        // a direcao local (tetra) deve respeitar a orientacao global do DOF.
        // O fator sgn_i*sgn_j preserva continuidade tangencial entre elementos.
        const double s = (double)(sgn_i * sgn_j);
        add_global(I, J, s * Sel[li][lj], s * Tel[li][lj]);
      }
    }
  }
}

} // namespace

/******************************************************************************/
/* FUNCAO: build_helm3d_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial 3D de aresta para autovalores de       */
/* cavidades. Wrapper homogeneo com BC PEC padrao.                            */
/* ENTRADA: mesh: const Mesh3D &; eps_r: double; mu_r: double.                */
/* SAIDA: EdgeSystem3D.                                                       */
/******************************************************************************/
EdgeSystem3D build_helm3d_edge_system(
    const Mesh3D &mesh,
    double eps_r,
    double mu_r)
{
  // Wrapper para meio homogeneo: expande eps_r/mu_r para vetores por tetra.
  auto eps = uniform_data((int)mesh.tets.size(), eps_r);
  auto mu = uniform_data((int)mesh.tets.size(), mu_r);
  return build_helm3d_edge_system(mesh, Edge3DBC::PEC_TangentialZero, eps, mu);
}

/******************************************************************************/
/* FUNCAO: build_helm3d_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial 3D de aresta para autovalores de       */
/* cavidades. Wrapper homogeneo com BC explicita.                             */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r: double; mu_r: double.  */
/* SAIDA: EdgeSystem3D.                                                       */
/******************************************************************************/
EdgeSystem3D build_helm3d_edge_system(
    const Mesh3D &mesh,
    Edge3DBC bc,
    double eps_r,
    double mu_r)
{
  // Wrapper homogeneo com BC explicita.
  auto eps = uniform_data((int)mesh.tets.size(), eps_r);
  auto mu = uniform_data((int)mesh.tets.size(), mu_r);
  return build_helm3d_edge_system(mesh, bc, eps, mu);
}

/******************************************************************************/
/* FUNCAO: build_helm3d_edge_system                                           */
/* DESCRICAO: Monta o sistema vetorial 3D de aresta para autovalores de       */
/* cavidades em meio heterogeneo, conforme formulacao vetorial da Secao 3.1.  */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r_tet: const              */
/* std::vector<double> &; mu_r_tet: const std::vector<double> &.              */
/* SAIDA: EdgeSystem3D.                                                       */
/******************************************************************************/
EdgeSystem3D build_helm3d_edge_system(
    const Mesh3D &mesh,
    Edge3DBC bc,
    const std::vector<double> &eps_r_tet,
    const std::vector<double> &mu_r_tet)
{
  if ((int)eps_r_tet.size() != (int)mesh.tets.size())
    throw std::runtime_error("eps_r_tet.size() != mesh.tets.size()");
  if ((int)mu_r_tet.size() != (int)mesh.tets.size())
    throw std::runtime_error("mu_r_tet.size() != mesh.tets.size()");

  EdgeSystem3D out;
  out.ed = build_edge_dofs_3d(mesh, bc);
  out.S = DenseMat(out.ed.ndof);
  out.T = DenseMat(out.ed.ndof);

  // Montagem global densa do problema Sx=lambdaTx.
  assemble_generic(
      mesh,
      out.ed,
      eps_r_tet,
      mu_r_tet,
      [&](int I, int J, double s_ij, double t_ij)
      {
        out.S(I, J) += s_ij;
        out.T(I, J) += t_ij;
      });

  return out;
}

/******************************************************************************/
/* FUNCAO: build_helm3d_edge_system_sparse                                    */
/* DESCRICAO: Monta o sistema vetorial 3D em formato esparso simetrico        */
/* (triangulo inferior). Implementa a formulacao vetorial 3D da Secao 3.1.    */
/* ENTRADA: mesh: const Mesh3D &; bc: Edge3DBC; eps_r_tet: const              */
/* std::vector<double> &; mu_r_tet: const std::vector<double> &.              */
/* SAIDA: EdgeSystem3DSparse.                                                 */
/******************************************************************************/
EdgeSystem3DSparse build_helm3d_edge_system_sparse(
    const Mesh3D &mesh,
    Edge3DBC bc,
    const std::vector<double> &eps_r_tet,
    const std::vector<double> &mu_r_tet)
{
  if ((int)eps_r_tet.size() != (int)mesh.tets.size())
    throw std::runtime_error("eps_r_tet.size() != mesh.tets.size()");
  if ((int)mu_r_tet.size() != (int)mesh.tets.size())
    throw std::runtime_error("mu_r_tet.size() != mesh.tets.size()");

  EdgeSystem3DSparse out;
  out.ed = build_edge_dofs_3d(mesh, bc);
  out.S = SparseSymMat(out.ed.ndof);
  out.T = SparseSymMat(out.ed.ndof);

  // Montagem global esparsa simetrica (triangulo inferior).
  assemble_generic(
      mesh,
      out.ed,
      eps_r_tet,
      mu_r_tet,
      [&](int I, int J, double s_ij, double t_ij)
      {
        // assemble_generic visita (I,J) e (J,I); no armazenamento simetrico
        // guardamos apenas a parte inferior para nao dobrar off-diagonais.
        if (I >= J)
        {
          out.S.add(I, J, s_ij);
          out.T.add(I, J, t_ij);
        }
      });

  return out;
}
