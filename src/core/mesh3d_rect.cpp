/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/core/mesh3d_rect.cpp                                          */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Geracao de malha tetraedrica retangular e perfis de material    */
/* para os casos 3D de validacao (Figura 15 e Figura 16).                     */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Apoio as Secoes 2.x */
/* e 3.1.                                                                     */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "mesh3d_rect.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace
{
/******************************************************************************/
/* FUNCAO: signed_tet_volume                                                  */
/* DESCRICAO: Calcula o volume assinado do tetraedro para validar orientacao e detectar degeneracao geometrica. */
/* Sinal positivo define orientacao consistente da conectividade local.        */
/* ENTRADA: a: const Node3D &; b: const Node3D &; c: const Node3D &; d: const */
/* Node3D &.                                                                  */
/* SAIDA: double.                                                             */
/******************************************************************************/
double signed_tet_volume(
    const Node3D &a,
    const Node3D &b,
    const Node3D &c,
    const Node3D &d)
{
  const double bax = b.x - a.x, bay = b.y - a.y, baz = b.z - a.z;
  const double cax = c.x - a.x, cay = c.y - a.y, caz = c.z - a.z;
  const double dax = d.x - a.x, day = d.y - a.y, daz = d.z - a.z;
  const double det =
      bax * (cay * daz - caz * day) -
      bay * (cax * daz - caz * dax) +
      baz * (cax * day - cay * dax);
  return det / 6.0;
}
} // namespace

/******************************************************************************/
/* FUNCAO: make_rect_tet_mesh                                                 */
/* DESCRICAO: Gera malha tetraedrica estruturada para paralelepipedo          */
/* retangular com 6 tetraedros por celula hexaedrica (Freudenthal).           */
/* Esta malha e a base dos casos da Figura 15 (ar) e Figura 16 (meio          */
/* preenchido).                                                               */
/* ENTRADA: lx: double; ly: double; lz: double; nx: int; ny: int; nz: int.    */
/* SAIDA: Mesh3D.                                                             */
/******************************************************************************/
Mesh3D make_rect_tet_mesh(double lx, double ly, double lz, int nx, int ny, int nz)
{
  if (nx <= 0 || ny <= 0 || nz <= 0)
    throw std::runtime_error("nx, ny, nz devem ser > 0.");
  if (!(lx > 0.0 && ly > 0.0 && lz > 0.0))
    throw std::runtime_error("Dimensoes lx, ly, lz devem ser > 0.");

  Mesh3D mesh;
  mesh.nodes.reserve((size_t)(nx + 1) * (ny + 1) * (nz + 1));
  mesh.tets.reserve((size_t)6 * nx * ny * nz);

  auto node_id = [&](int i, int j, int k)
  { return (k * (ny + 1) + j) * (nx + 1) + i; };

  for (int k = 0; k <= nz; ++k)
  {
    const double z = lz * (double)k / (double)nz;
    for (int j = 0; j <= ny; ++j)
    {
      const double y = ly * (double)j / (double)ny;
      for (int i = 0; i <= nx; ++i)
      {
        const double x = lx * (double)i / (double)nx;
        const bool is_bd = (i == 0 || i == nx || j == 0 || j == ny || k == 0 || k == nz);
        mesh.nodes.push_back({x, y, z, is_bd});
      }
    }
  }

  for (int k = 0; k < nz; ++k)
  {
    for (int j = 0; j < ny; ++j)
    {
      for (int i = 0; i < nx; ++i)
      {
        const int v0 = node_id(i, j, k);
        const int v1 = node_id(i + 1, j, k);
        const int v2 = node_id(i, j + 1, k);
        const int v3 = node_id(i + 1, j + 1, k);
        const int v4 = node_id(i, j, k + 1);
        const int v5 = node_id(i + 1, j, k + 1);
        const int v6 = node_id(i, j + 1, k + 1);
        const int v7 = node_id(i + 1, j + 1, k + 1);

        // Particionamento de Freudenthal em 6 tetraedros com diagonal espacial v0 -> v7.
        // Esse padrao evita tetraedros degenerados e garante preenchimento
        // conformante de toda a celula hexaedrica.
        std::array<std::array<int, 4>, 6> local = {{
            {{v0, v1, v3, v7}},
            {{v0, v3, v2, v7}},
            {{v0, v2, v6, v7}},
            {{v0, v6, v4, v7}},
            {{v0, v4, v5, v7}},
            {{v0, v5, v1, v7}},
        }};

        for (auto tetv : local)
        {
          Tet t = {{tetv[0], tetv[1], tetv[2], tetv[3]}};
          const Node3D &a = mesh.nodes[t.v[0]];
          const Node3D &b = mesh.nodes[t.v[1]];
          const Node3D &c = mesh.nodes[t.v[2]];
          const Node3D &d = mesh.nodes[t.v[3]];
          // Ajusta orientacao local para volume positivo, requisito para
          // manter sinais consistentes nas integrais elementares.
          if (signed_tet_volume(a, b, c, d) < 0.0)
          {
            std::swap(t.v[2], t.v[3]);
          }
          mesh.tets.push_back(t);
        }
      }
    }
  }

  return mesh;
}

/******************************************************************************/
/* FUNCAO: make_eps_r_tets_by_z                                               */
/* DESCRICAO: Atribui permissividade por tetraedro com base na posicao do     */
/* centroide em z. Modelo de meio por partes para o caso meio preenchido      */
/* (Figura 16 / Tabela 13).                                                   */
/* ENTRADA: mesh: const Mesh3D &; z_cut: double; eps_r_below: double;         */
/* eps_r_above: double.                                                       */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<double> make_eps_r_tets_by_z(
    const Mesh3D &mesh,
    double z_cut,
    double eps_r_below,
    double eps_r_above)
{
  std::vector<double> eps((size_t)mesh.tets.size(), eps_r_above);
  for (int tid = 0; tid < (int)mesh.tets.size(); ++tid)
  {
    const Tet &t = mesh.tets[tid];
    const Node3D &n0 = mesh.nodes[t.v[0]];
    const Node3D &n1 = mesh.nodes[t.v[1]];
    const Node3D &n2 = mesh.nodes[t.v[2]];
    const Node3D &n3 = mesh.nodes[t.v[3]];
    const double zc = 0.25 * (n0.z + n1.z + n2.z + n3.z);
    eps[tid] = (zc < z_cut) ? eps_r_below : eps_r_above;
  }
  return eps;
}
