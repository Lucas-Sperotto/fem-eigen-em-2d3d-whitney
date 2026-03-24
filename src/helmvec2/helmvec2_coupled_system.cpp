/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helmvec2/helmvec2_coupled_system.cpp                          */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Sistema acoplado vetorial+escalar para obter k0 dado beta.      */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Secao 2.2.3, Eq.    */
/* (108)-(109), Fig. 11, Tabela 8.                                            */
/*****************************************************************************/
/* Observacao: Comentarios priorizam didatica, rastreabilidade e validacao.   */
/*****************************************************************************/

#include "helmvec2_coupled_system.hpp"
#include "edge/edge_basis.hpp"
#include "explicit/tri2d_coupled_explicit.hpp"
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace
{
// Quadratura triangular de 3 pontos (exata para grau 2 em elementos lineares).
// Usada no bloco de acoplamento vetorial-escalar:
//   C_mj = int_T (1/mu_r) W_m . grad(N_j) dA
constexpr std::array<std::array<double, 3>, 3> kTriQuadP2 = {{
    {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
    {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
}};

struct CoupledContextE
{
    EdgeSystem edge;
    ScalarSystem scal;
    DenseMat mt_muinv; // int (1/mu_r) W_i.W_j dA na base de arestas
    int nt = 0;
    int nz = 0;
};

/******************************************************************************/
/* FUNCAO: validate_tri_data                                                  */
/* DESCRICAO: Valida consistencia geometrica/topologica dos triangulos antes da montagem do sistema acoplado. */
/* ENTRADA: mesh: const Mesh2D &; eps_r_tri: const std::vector<double> &;     */
/* mu_r_tri: const std::vector<double> &.                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void validate_tri_data(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri)
{
    if ((int)eps_r_tri.size() != (int)mesh.tris.size())
        throw std::runtime_error("eps_r_tri.size() != mesh.tris.size()");
    if ((int)mu_r_tri.size() != (int)mesh.tris.size())
        throw std::runtime_error("mu_r_tri.size() != mesh.tris.size()");

    for (int i = 0; i < (int)mesh.tris.size(); ++i)
    {
        if (eps_r_tri[(size_t)i] <= 0.0)
            throw std::runtime_error("eps_r_tri deve ser positivo em todo triangulo.");
        if (std::abs(mu_r_tri[(size_t)i]) < 1e-14)
            throw std::runtime_error("mu_r_tri contem valor muito proximo de zero.");
    }
}

/******************************************************************************/
/* FUNCAO: inverse_per_triangle_checked                                       */
/* DESCRICAO: Calcula inversa local por triangulo com verificacoes de singularidade e robustez numerica. */
/* ENTRADA: v: const std::vector<double> &.                                   */
/* SAIDA: std::vector<double>.                                                */
/******************************************************************************/
std::vector<double> inverse_per_triangle_checked(const std::vector<double> &v)
{
    std::vector<double> out(v.size(), 0.0);
    for (size_t i = 0; i < v.size(); ++i)
    {
        if (v[i] == 0.0)
            throw std::runtime_error("mu_r_tri contem zero.");
        out[i] = 1.0 / v[i];
    }
    return out;
}

/******************************************************************************/
/* FUNCAO: add_block_scaled                                                   */
/* DESCRICAO: Acumula bloco local escalonado na matriz global, preservando simetria e sinal fisico. */
/* ENTRADA: dst: DenseMat &; row0: int; col0: int; src: const DenseMat &;     */
/* alpha: double.                                                             */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void add_block_scaled(
    DenseMat &dst,
    int row0,
    int col0,
    const DenseMat &src,
    double alpha)
{
    const int n = src.n;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            dst(row0 + i, col0 + j) += alpha * src(i, j);
}

/******************************************************************************/
/* FUNCAO: build_edge_mass_with_inverse_mu                                    */
/* DESCRICAO: Constroi matriz de massa vetorial ponderada por 1/mu_r, isto e, */
/* M_t^(1/mu)=int_T (1/mu_r) W_i.W_j dA, usada nos termos com beta^2 das      */
/* Eq. (113), (116) e no rearranjo da Eq. (126)-(127).                        */
/* ENTRADA: mesh: const Mesh2D &; bc: EdgeBC; mu_r_tri: const                 */
/* std::vector<double> &.                                                     */
/* SAIDA: DenseMat.                                                           */
/******************************************************************************/
DenseMat build_edge_mass_with_inverse_mu(
    const Mesh2D &mesh,
    EdgeBC bc,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend)
{
    // Reuso do montador de aresta configurando:
    //   eps_proxy = 1/mu_r  -> a massa T vira int (1/mu_r) W_i.W_j dA
    //   mu_proxy  = 1       -> sem escala extra no termo curl-curl
    std::vector<double> one(mesh.tris.size(), 1.0);
    auto inv_mu = inverse_per_triangle_checked(mu_r_tri);
    auto sys = build_helm10_edge_system(mesh, bc, inv_mu, one, backend);
    return sys.T;
}

/******************************************************************************/
/* FUNCAO: build_context_E                                                    */
/* DESCRICAO: Prepara blocos base da formulacao em E (Et,Ez): sistema de      */
/* aresta, sistema escalar e massa vetorial com 1/mu para o acoplamento das   */
/* Secoes 2.2.3 e 2.2.4.                                                       */
/* ENTRADA: mesh: const Mesh2D &; eps_r_tri: const std::vector<double> &;     */
/* mu_r_tri: const std::vector<double> &.                                     */
/* SAIDA: CoupledContextE.                                                    */
/******************************************************************************/
CoupledContextE build_context_E(
    const Mesh2D &mesh,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend)
{
    validate_tri_data(mesh, eps_r_tri, mu_r_tri);

    CoupledContextE ctx;
    ctx.edge = build_helm10_edge_system(
        mesh,
        EdgeBC::TE_PEC_TangentialZero,
        eps_r_tri,
        mu_r_tri,
        backend);
    ctx.scal = build_helm10_scalar_system(
        mesh,
        ScalarBC::TM_Dirichlet,
        eps_r_tri,
        mu_r_tri,
        backend);
    ctx.nt = ctx.edge.ed.ndof;
    ctx.nz = ctx.scal.ndof;
    ctx.mt_muinv = build_edge_mass_with_inverse_mu(
        mesh,
        EdgeBC::TE_PEC_TangentialZero,
        mu_r_tri,
        backend);
    return ctx;
}

/******************************************************************************/
/* FUNCAO: assemble_coupling_block_C                                          */
/* DESCRICAO: Monta o bloco de acoplamento vetorial-escalar e aplica          */
/* orientacao local/global. Este bloco representa int(1/mu) W_m.grad(N_j)dA   */
/* que aparece nos termos mistos das Eq. (114)-(115) e Eq. (126)-(127).       */
/* ENTRADA: full: DenseMat &; nt: int; coeff_top_right: double;               */
/* coeff_bottom_left: double; mesh: const Mesh2D &; ed: const EdgeDofs &;     */
/* node_to_dof: const std::vector<int> &; mu_r_tri: const std::vector<double> */
/* &.                                                                         */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void assemble_coupling_block_C(
    DenseMat &full,
    int nt,
    double coeff_top_right,
    double coeff_bottom_left,
    const Mesh2D &mesh,
    const EdgeDofs &ed,
    const std::vector<int> &node_to_dof,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend)
{
    // Monta o bloco C e seu acoplamento transposto com orientacao:
    //   C_mj = int (1/mu_r) W_m . grad(N_j) dA
    // e escreve:
    //   bloco sup-dir  += coeff_top_right  * C
    //   bloco inf-esq  += coeff_bottom_left* C_orientado^T
    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const Tri &t = mesh.tris[(size_t)tid];
        const TriGeomEdge tg = tri_geom_edge(mesh, t);
        const TriEdges &te = ed.tri_edges[(size_t)tid];

        Vec2 gradN[3];
        for (int j = 0; j < 3; ++j)
        {
            gradN[j].x = tg.g.b[j] / (2.0 * tg.g.A);
            gradN[j].y = tg.g.c[j] / (2.0 * tg.g.A);
        }

        double C_loc[3][3] = {{0.0}};
        const double inv_mu = 1.0 / mu_r_tri[(size_t)tid];
        if (backend == ElementAssemblyBackend::ClosedForm)
        {
            const auto edge_cf = explicit_tri2d::tri2d_edge_closed_form_data(tg);
            for (int m = 0; m < 3; ++m)
            {
                for (int j = 0; j < 3; ++j)
                {
                    // Eq. (121)-(122), sem o fator beta^2: bloco base C.
                    C_loc[m][j] =
                        inv_mu *
                        explicit_tri2d::tri2d_edge_scalar_coupling_entry_base(
                            edge_cf,
                            m,
                            j);
                }
            }
        }
        else
        {
            for (const auto &lam : kTriQuadP2)
            {
                const double w = inv_mu * (tg.g.A / 3.0);
                Vec2 W[3];
                for (int m = 0; m < 3; ++m)
                    W[m] = whitney_W_local(m, tg, lam);

                for (int m = 0; m < 3; ++m)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        C_loc[m][j] += w * (W[m].x * gradN[j].x + W[m].y * gradN[j].y);
                    }
                }
            }
        }

        for (int m = 0; m < 3; ++m)
        {
            int eid = te.e[m];
            int sgn = te.sgn[m];
            int I = ed.edge_to_dof[(size_t)eid];
            if (I < 0)
                continue;

            for (int j = 0; j < 3; ++j)
            {
                int node = t.v[j];
                int J = node_to_dof[(size_t)node];
                if (J < 0)
                    continue;

                const double cij = (double)sgn * C_loc[m][j];
                full(I, nt + J) += coeff_top_right * cij;
                full(nt + J, I) += coeff_bottom_left * cij;
            }
        }
    }
}

/******************************************************************************/
/* FUNCAO: assemble_wavenumber_system_closed_form                             */
/* DESCRICAO: Monta diretamente o sistema acoplado A x = k0^2 B x usando as   */
/* formas locais fechadas da Secao 2.2.3, isto e, Eq. (120) a (125). Esta     */
/* rotina faz a ponte entre esses blocos locais e a assembleia global da Eq.  */
/* (119), no arranjo validado do repositorio.                                 */
/* ENTRADA: out: CoupledWaveNumberSystem &; mesh: const Mesh2D &; beta:       */
/* double; eps_r_tri: const std::vector<double> &; mu_r_tri: const            */
/* std::vector<double> &.                                                     */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void assemble_wavenumber_system_closed_form(
    CoupledWaveNumberSystem &out,
    const Mesh2D &mesh,
    double beta,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri)
{
    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const Tri &t = mesh.tris[(size_t)tid];
        const TriGeomEdge tg = tri_geom_edge(mesh, t);
        const TriEdges &te = out.edge.ed.tri_edges[(size_t)tid];
        const auto blk = explicit_tri2d::tri2d_wavenumber_rearranged_closed_form_eq_119(
            tg,
            beta,
            eps_r_tri[(size_t)tid],
            mu_r_tri[(size_t)tid]);

        // O helper local ja retorna os blocos no formato do sistema global:
        //   A x = k0^2 B x
        // com o seguinte mapeamento didatico:
        //   A_tt <- Eq. (120)
        //   A_tz <- Eq. (121)
        //   A_zt <- Eq. (122)
        //   A_zz <- Eq. (123)
        //   B_tt <- Eq. (124)
        //   B_zz <- Eq. (125)

        // Eq. (120) e Eq. (124): blocos Et-Et em A e B.
        for (int m = 0; m < 3; ++m)
        {
            const int edge_m = te.e[m];
            const int sgn_m = te.sgn[m];
            const int I = out.edge.ed.edge_to_dof[(size_t)edge_m];
            if (I < 0)
                continue;

            for (int n = 0; n < 3; ++n)
            {
                const int edge_n = te.e[n];
                const int sgn_n = te.sgn[n];
                const int J = out.edge.ed.edge_to_dof[(size_t)edge_n];
                if (J < 0)
                    continue;

                const double sign = static_cast<double>(sgn_m * sgn_n);
                out.A(I, J) += sign * blk.A_tt[m][n];
                out.B(I, J) += sign * blk.B_tt[m][n];
            }
        }

        // Eq. (123) e Eq. (125): blocos Ez-Ez em A e B.
        for (int i = 0; i < 3; ++i)
        {
            const int I = out.scal.dof_map[(size_t)t.v[i]];
            if (I < 0)
                continue;

            for (int j = 0; j < 3; ++j)
            {
                const int J = out.scal.dof_map[(size_t)t.v[j]];
                if (J < 0)
                    continue;

                out.A(out.nt + I, out.nt + J) += blk.A_zz[i][j];
                out.B(out.nt + I, out.nt + J) += blk.B_zz[i][j];
            }
        }

        // Eq. (121): bloco Et-Ez.
        for (int m = 0; m < 3; ++m)
        {
            const int edge_m = te.e[m];
            const int sgn_m = te.sgn[m];
            const int I = out.edge.ed.edge_to_dof[(size_t)edge_m];
            if (I < 0)
                continue;

            for (int j = 0; j < 3; ++j)
            {
                const int J = out.scal.dof_map[(size_t)t.v[j]];
                if (J < 0)
                    continue;

                out.A(I, out.nt + J) += static_cast<double>(sgn_m) * blk.A_tz[m][j];
            }
        }

        // Eq. (122): bloco Ez-Et.
        for (int i = 0; i < 3; ++i)
        {
            const int I = out.scal.dof_map[(size_t)t.v[i]];
            if (I < 0)
                continue;

            for (int n = 0; n < 3; ++n)
            {
                const int edge_n = te.e[n];
                const int sgn_n = te.sgn[n];
                const int J = out.edge.ed.edge_to_dof[(size_t)edge_n];
                if (J < 0)
                    continue;

                out.A(out.nt + I, J) += static_cast<double>(sgn_n) * blk.A_zt[i][n];
            }
        }
    }
}

/******************************************************************************/
/* FUNCAO: assemble_beta_system_closed_form                                   */
/* DESCRICAO: Monta diretamente o sistema rearranjado P x = beta^2 Q x usando */
/* os blocos locais closed-form coerentes com a Eq. (136). Esta rotina faz a  */
/* ponte entre as Eq. (137) a (142), calculadas localmente por triangulo, e   */
/* a montagem global rearranjada do problema em beta^2.                       */
/* ENTRADA: out: CoupledBetaSystem &; mesh: const Mesh2D &; k0: double;       */
/* eps_r_tri: const std::vector<double> &; mu_r_tri: const std::vector<double>*/
/* &.                                                                         */
/* SAIDA: sem retorno explicito (void).                                       */
/******************************************************************************/
void assemble_beta_system_closed_form(
    CoupledBetaSystem &out,
    const Mesh2D &mesh,
    double k0,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri)
{
    for (int tid = 0; tid < (int)mesh.tris.size(); ++tid)
    {
        const Tri &t = mesh.tris[(size_t)tid];
        const TriGeomEdge tg = tri_geom_edge(mesh, t);
        const TriEdges &te = out.edge.ed.tri_edges[(size_t)tid];
        const auto blk = explicit_tri2d::tri2d_beta_rearranged_closed_form_eq_136(
            tg,
            k0,
            eps_r_tri[(size_t)tid],
            mu_r_tri[(size_t)tid]);

        // O helper local ja retorna os blocos na forma do sistema global:
        //   P x = beta^2 Q x
        // com o seguinte mapeamento didatico:
        //   P_tt <- Eq. (137)
        //   P_zz <- parcela k0^2 eps_r M_z extraida da Eq. (142)
        //   Q_tt <- forma rearranjada validada para a Eq. (136)
        //   Q_tz <- Eq. (138)
        //   Q_zt <- Eq. (139)
        //   Q_zz <- Eq. (140)

        // Blocos Et-Et de P e Q.
        for (int m = 0; m < 3; ++m)
        {
            const int edge_m = te.e[m];
            const int sgn_m = te.sgn[m];
            const int I = out.edge.ed.edge_to_dof[(size_t)edge_m];
            if (I < 0)
                continue;

            for (int n = 0; n < 3; ++n)
            {
                const int edge_n = te.e[n];
                const int sgn_n = te.sgn[n];
                const int J = out.edge.ed.edge_to_dof[(size_t)edge_n];
                if (J < 0)
                    continue;

                const double sign = static_cast<double>(sgn_m * sgn_n);
                out.P(I, J) += sign * blk.P_tt[m][n];
                out.Q(I, J) += sign * blk.Q_tt[m][n];
            }
        }

        // Blocos Ez-Ez de P e Q.
        for (int i = 0; i < 3; ++i)
        {
            const int I = out.scal.dof_map[(size_t)t.v[i]];
            if (I < 0)
                continue;

            for (int j = 0; j < 3; ++j)
            {
                const int J = out.scal.dof_map[(size_t)t.v[j]];
                if (J < 0)
                    continue;

                out.P(out.nt + I, out.nt + J) += blk.P_zz[i][j];
                out.Q(out.nt + I, out.nt + J) += blk.Q_zz[i][j];
            }
        }

        // Bloco Et-Ez de Q.
        for (int m = 0; m < 3; ++m)
        {
            const int edge_m = te.e[m];
            const int sgn_m = te.sgn[m];
            const int I = out.edge.ed.edge_to_dof[(size_t)edge_m];
            if (I < 0)
                continue;

            for (int j = 0; j < 3; ++j)
            {
                const int J = out.scal.dof_map[(size_t)t.v[j]];
                if (J < 0)
                    continue;

                out.Q(I, out.nt + J) += static_cast<double>(sgn_m) * blk.Q_tz[m][j];
            }
        }

        // Bloco Ez-Et de Q.
        for (int i = 0; i < 3; ++i)
        {
            const int I = out.scal.dof_map[(size_t)t.v[i]];
            if (I < 0)
                continue;

            for (int n = 0; n < 3; ++n)
            {
                const int edge_n = te.e[n];
                const int sgn_n = te.sgn[n];
                const int J = out.edge.ed.edge_to_dof[(size_t)edge_n];
                if (J < 0)
                    continue;

                out.Q(out.nt + I, J) += static_cast<double>(sgn_n) * blk.Q_zt[i][n];
            }
        }
    }
}
} // namespace

/******************************************************************************/
/* FUNCAO: build_coupled_wavenumber_system_E                                  */
/* DESCRICAO: Monta o sistema acoplado A x = k0^2 B x para k0 dado beta.      */
/* Corresponde ao problema da Secao 2.2.3. Do ponto de vista do sistema       */
/* global ja montado sobre a malha, esta rotina materializa a Eq. (119),      */
/* usando os blocos locais das Eq. (120)-(125).                               */
/* ENTRADA: mesh: const Mesh2D &; beta: double; eps_r_tri: const              */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &.              */
/* SAIDA: CoupledWaveNumberSystem.                                            */
/******************************************************************************/
CoupledWaveNumberSystem build_coupled_wavenumber_system_E(
    const Mesh2D &mesh,
    double beta,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend)
{
    if (!std::isfinite(beta))
        throw std::runtime_error("beta deve ser finito.");

    auto ctx = build_context_E(mesh, eps_r_tri, mu_r_tri, backend);

    CoupledWaveNumberSystem out;
    out.edge = std::move(ctx.edge);
    out.scal = std::move(ctx.scal);
    out.nt = ctx.nt;
    out.nz = ctx.nz;

    const int N = out.nt + out.nz;
    out.A = DenseMat(N);
    out.B = DenseMat(N);

    const double beta2 = beta * beta;

    if (backend == ElementAssemblyBackend::ClosedForm)
    {
        // Caminho explicitamente rastreavel pelas Eq. (120)-(125).
        assemble_wavenumber_system_closed_form(out, mesh, beta, eps_r_tri, mu_r_tri);
        return out;
    }

    const DenseMat &Bt = ctx.mt_muinv;

    // Sistema global final da Secao 2.2.3:
    //   [S_tt + beta^2 Mt^(1/mu)      beta^2 C ] [Et] = k0^2 [T_tt      0      ] [Et]
    //   [beta^2 C^T                   beta^2 Gz] [Ez]        [0         beta^2 Tz] [Ez]
    //
    // onde C = int(1/mu) W.grad(N) e Gz = grad-grad escalar com 1/mu.
    // Na numeracao do artigo, este e o sistema global da Eq. (119).

    // Eq. (120): Sel(tt) = (1/mu)curlcurl + beta^2 (1/mu) massa_t.
    add_block_scaled(out.A, 0, 0, out.edge.S, +1.0);
    add_block_scaled(out.A, 0, 0, Bt, +beta2);

    // Eq. (124): Tel(tt) = eps_r * massa_t.
    add_block_scaled(out.B, 0, 0, out.edge.T, +1.0);

    // Eq. (123): Sel(zz) = beta^2 (1/mu) * grad-grad.
    add_block_scaled(out.A, out.nt, out.nt, out.scal.S, +beta2);

    // Eq. (125): Tel(zz) = beta^2 * eps_r * massa_z.
    add_block_scaled(out.B, out.nt, out.nt, out.scal.T, +beta2);

    // Eq. (121)-(122): Sel(tz) e Sel(zt) com acoplamento +beta^2/mu.
    assemble_coupling_block_C(
        out.A,
        out.nt,
        +beta2,
        +beta2,
        mesh,
        out.edge.ed,
        out.scal.dof_map,
        mu_r_tri,
        backend);

    return out;
}

/******************************************************************************/
/* FUNCAO: build_coupled_beta_system_E                                        */
/* DESCRICAO: Monta o sistema acoplado P x = beta^2 Q x para beta dado k0.    */
/* Corresponde ao problema da Secao 2.2.4. Do ponto de vista do sistema       */
/* global ja montado sobre a malha, esta rotina materializa a Eq. (136),      */
/* usando os blocos locais das Eq. (137)-(142).                               */
/* ENTRADA: mesh: const Mesh2D &; k0: double; eps_r_tri: const                */
/* std::vector<double> &; mu_r_tri: const std::vector<double> &.              */
/* SAIDA: CoupledBetaSystem.                                                  */
/******************************************************************************/
CoupledBetaSystem build_coupled_beta_system_E(
    const Mesh2D &mesh,
    double k0,
    const std::vector<double> &eps_r_tri,
    const std::vector<double> &mu_r_tri,
    ElementAssemblyBackend backend)
{
    if (!std::isfinite(k0))
        throw std::runtime_error("k0 deve ser finito.");

    auto ctx = build_context_E(
        mesh,
        eps_r_tri,
        mu_r_tri,
        backend);
    

    CoupledBetaSystem out;
    out.edge = std::move(ctx.edge);
    out.scal = std::move(ctx.scal);
    out.nt = ctx.nt;
    out.nz = ctx.nz;

    const int N = out.nt + out.nz;
    out.P = DenseMat(N);
    out.Q = DenseMat(N);

    const double k02 = k0 * k0;

    if (backend == ElementAssemblyBackend::ClosedForm)
    {
        // Caminho explicitamente rastreavel pela Eq. (136) e pelos blocos
        // locais Eq. (137)-(142), no rearranjo validado do repositorio.
        assemble_beta_system_closed_form(out, mesh, k0, eps_r_tri, mu_r_tri);
        return out;
    }

    // Conjunto de equacoes (126)-(127), rearranjado como:
    //   P x = beta^2 Q x
    //
    // com
    //   P_tt = (1/mu)curlcurl - k0^2 eps*Mt
    //   P_zz = k0^2 eps*Mz
    //
    //   Q_tt = -(1/mu)Mt
    //   Q_tz = +(1/mu)C
    //   Q_zt = +(1/mu)C^T
    //   Q_zz = (1/mu)Gz
    //
    // onde Mt e a massa de aresta, Mz a massa escalar, Gz o grad-grad escalar
    // e C a matriz de acoplamento vetorial-escalar.
    //
    // Forma matricial usada:
    //   P = [ S_tt - k0^2 T_tt      0          ]
    //       [ 0                      k0^2 T_zz ]
    //   Q = [ -Mt^(1/mu)            +C         ]
    //       [ +C^T                  +Gz        ]
    //
    // Na numeracao do artigo, este e o sistema global da Eq. (136).

    // Eq. (137), rearranjada para P_tt:
    //   (1/mu)curlcurl - k0^2 eps_r Mt.
    add_block_scaled(out.P, 0, 0, out.edge.S, +1.0);
    add_block_scaled(out.P, 0, 0, out.edge.T, -k02);

    // Eq. (142), termo de massa escalar com k0^2 eps_r, levado para P_zz.
    add_block_scaled(out.P, out.nt, out.nt, out.scal.T, +k02);

    // Eq. (141), rearranjada para Q_tt com sinal negativo.
    add_block_scaled(out.Q, 0, 0, ctx.mt_muinv, -1.0);

    // Eq. (140), rearranjada para Q_zz = (1/mu) grad-grad.
    add_block_scaled(out.Q, out.nt, out.nt, out.scal.S, +1.0);

    // Eq. (138)-(139): blocos mistos Q_tz e Q_zt.
    assemble_coupling_block_C(
        out.Q,
        out.nt,
        +1.0,
        +1.0,
        mesh,
        out.edge.ed,
        out.scal.dof_map,
        mu_r_tri,
        backend);

    return out;
}
