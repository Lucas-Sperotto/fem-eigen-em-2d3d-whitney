/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/explicit/tri2d_coupled_explicit.hpp                           */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 1.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Formas fechadas locais do sistema acoplado vetorial+escalar 2D  */
/* para o problema de k0 dado beta (Secao 2.2.3).                             */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Eq. (120)-(125).    */
/*****************************************************************************/
/* Observacao: Este cabecalho reaproveita as formas fechadas escalares e de   */
/* aresta ja implementadas nas secoes anteriores.                             */
/*****************************************************************************/

#pragma once

#include "explicit/tri2d_edge_explicit.hpp"
#include "explicit/tri2d_scalar_explicit.hpp"

namespace explicit_tri2d
{

struct Tri2DCoupledWaveNumberLocalBlocks
{
    double Sel_tt[3][3] = {{0.0}};
    double Sel_tz[3][3] = {{0.0}};
    double Sel_zt[3][3] = {{0.0}};
    double Sel_zz[3][3] = {{0.0}};
    double Tel_tt[3][3] = {{0.0}};
    double Tel_zz[3][3] = {{0.0}};
};

struct Tri2DCoupledWaveNumberRearrangedLocalBlocks
{
    double A_tt[3][3] = {{0.0}};
    double A_tz[3][3] = {{0.0}};
    double A_zt[3][3] = {{0.0}};
    double A_zz[3][3] = {{0.0}};
    double B_tt[3][3] = {{0.0}};
    double B_zz[3][3] = {{0.0}};
};

struct Tri2DCoupledBetaLocalBlocks
{
    double Sel_tt[3][3] = {{0.0}};
    double Tel_tz[3][3] = {{0.0}};
    double Tel_zt[3][3] = {{0.0}};
    double Sel_zz[3][3] = {{0.0}};
    double Tel_tt[3][3] = {{0.0}};
    double Tel_zz[3][3] = {{0.0}};
};

struct Tri2DCoupledBetaRearrangedLocalBlocks
{
    double P_tt[3][3] = {{0.0}};
    double P_zz[3][3] = {{0.0}};
    double Q_tt[3][3] = {{0.0}};
    double Q_tz[3][3] = {{0.0}};
    double Q_zt[3][3] = {{0.0}};
    double Q_zz[3][3] = {{0.0}};
};

/******************************************************************************/
/* FUNCAO: tri2d_edge_scalar_coupling_entry_base                              */
/* DESCRICAO: Calcula a integral base C_mj = int_T W_m.grad(N_j) dA, sem os   */
/* fatores beta^2 e 1/mu_r. Esta expressao e a base fechada das Eq. (121) e   */
/* (122) da Secao 2.2.3 e tambem das Eq. (138) e (139) da Secao 2.2.4. Em     */
/* outras palavras: a mesma forma explicita local e reaproveitada nos dois     */
/* problemas acoplados, mudando apenas o fator global e o lado do sistema.    */
/* ENTRADA: data: const Tri2DEdgeClosedFormData &; j: int; m: int.            */
/* SAIDA: double.                                                             */
/******************************************************************************/
inline double tri2d_edge_scalar_coupling_entry_base(
    const Tri2DEdgeClosedFormData &data,
    int m,
    int j)
{
    const double scale = data.edge_len[m] / (8.0 * data.area * data.area);
    // Base algebrica comum das Eq. (121) e (122), antes dos fatores beta^2/mu_r.
    return scale *
           (data.b[j] * (data.Acoef[m] + data.Bcoef[m] * data.ytri) +
            data.c[j] * (data.Ccoef[m] + data.Dcoef[m] * data.xtri));
}

/******************************************************************************/
/* FUNCAO: tri2d_wavenumber_closed_form_eq_120_125                            */
/* DESCRICAO: Monta os seis sub-blocos locais da formulacao 2.2.3 para obter  */
/* k0 dado beta. A funcao codifica diretamente as Eq. (120) a (125) e         */
/* explicita, por reuso, os mesmos nucleos locais depois reaproveitados no    */
/* problema de 2.2.4 (Eq. (137) a (142)), especialmente os termos base de     */
/* aresta, grad-grad escalar e acoplamento misto.                             */
/* ENTRADA: tg: const TriGeomEdge &; beta: double; eps_r: double; mu_r:       */
/* double.                                                                    */
/* SAIDA: Tri2DCoupledWaveNumberLocalBlocks.                                  */
/******************************************************************************/
inline Tri2DCoupledWaveNumberLocalBlocks tri2d_wavenumber_closed_form_eq_120_125(
    const TriGeomEdge &tg,
    double beta,
    double eps_r,
    double mu_r)
{
    Tri2DCoupledWaveNumberLocalBlocks blk;

    const double beta2 = beta * beta;
    const double inv_mu = 1.0 / mu_r;
    const Tri2DEdgeClosedFormData edge = tri2d_edge_closed_form_data(tg);

    double scalar_sel[3][3];
    double scalar_tel[3][3];
    tri2d_scalar_closed_form_eq_30_33(tg.g.A, tg.g.b, tg.g.c, scalar_sel, scalar_tel);

    for (int m = 0; m < 3; ++m)
    {
        for (int n = 0; n < 3; ++n)
        {
            // Eq. (120): Sel(tt) = (1/mu) curl-curl + beta^2 (1/mu) massa_t.
            blk.Sel_tt[m][n] =
                inv_mu * tri2d_edge_curlcurl_entry_eq_66(edge, m, n) +
                beta2 * inv_mu * tri2d_edge_mass_entry_eq_67(edge, m, n);

            // Eq. (124): Tel(tt) = eps_r * massa_t.
            blk.Tel_tt[m][n] = eps_r * tri2d_edge_mass_entry_eq_67(edge, m, n);
        }
    }

    for (int m = 0; m < 3; ++m)
    {
        for (int j = 0; j < 3; ++j)
        {
            const double c_mj = tri2d_edge_scalar_coupling_entry_base(edge, m, j);

            // Eq. (121): Sel(tz).
            blk.Sel_tz[m][j] = beta2 * inv_mu * c_mj;

            // Eq. (122): Sel(zt) = Sel(tz)^T na notacao local.
            blk.Sel_zt[j][m] = beta2 * inv_mu * c_mj;
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            // Eq. (123): Sel(zz) = beta^2 (1/mu) grad-grad escalar.
            blk.Sel_zz[i][j] = beta2 * inv_mu * scalar_sel[i][j];

            // Eq. (125): Tel(zz) = beta^2 eps_r * massa escalar.
            blk.Tel_zz[i][j] = beta2 * eps_r * scalar_tel[i][j];
        }
    }

    return blk;
}

/******************************************************************************/
/* FUNCAO: tri2d_wavenumber_rearranged_closed_form_eq_119                     */
/* DESCRICAO: Reorganiza os blocos locais da Secao 2.2.3 no mesmo arranjo     */
/* usado pelo sistema global A x = k0^2 B x. Esta funcao deixa explicito o    */
/* elo didatico entre as Eq. (120) a (125) e a Eq. (119).                     */
/* ENTRADA: tg: const TriGeomEdge &; beta: double; eps_r: double; mu_r:       */
/* double.                                                                    */
/* SAIDA: Tri2DCoupledWaveNumberRearrangedLocalBlocks.                        */
/******************************************************************************/
inline Tri2DCoupledWaveNumberRearrangedLocalBlocks tri2d_wavenumber_rearranged_closed_form_eq_119(
    const TriGeomEdge &tg,
    double beta,
    double eps_r,
    double mu_r)
{
    Tri2DCoupledWaveNumberRearrangedLocalBlocks blk119;
    const auto blk = tri2d_wavenumber_closed_form_eq_120_125(tg, beta, eps_r, mu_r);

    for (int m = 0; m < 3; ++m)
    {
        for (int n = 0; n < 3; ++n)
        {
            blk119.A_tt[m][n] = blk.Sel_tt[m][n];
            blk119.B_tt[m][n] = blk.Tel_tt[m][n];
            blk119.A_tz[m][n] = blk.Sel_tz[m][n];
            blk119.A_zt[m][n] = blk.Sel_zt[m][n];
            blk119.A_zz[m][n] = blk.Sel_zz[m][n];
            blk119.B_zz[m][n] = blk.Tel_zz[m][n];
        }
    }

    return blk119;
}

/******************************************************************************/
/* FUNCAO: tri2d_beta_closed_form_eq_137_142                                  */
/* DESCRICAO: Monta os seis sub-blocos locais da formulacao 2.2.4 para obter  */
/* beta dado k0. A funcao codifica de modo explicito as Eq. (137) a (142),    */
/* usando os mesmos nucleos closed-form ja reutilizados em 2.2.3.             */
/* ENTRADA: tg: const TriGeomEdge &; k0: double; eps_r: double; mu_r: double. */
/* SAIDA: Tri2DCoupledBetaLocalBlocks.                                        */
/******************************************************************************/
inline Tri2DCoupledBetaLocalBlocks tri2d_beta_closed_form_eq_137_142(
    const TriGeomEdge &tg,
    double k0,
    double eps_r,
    double mu_r)
{
    Tri2DCoupledBetaLocalBlocks blk;

    const double k02 = k0 * k0;
    const double inv_mu = 1.0 / mu_r;
    const Tri2DEdgeClosedFormData edge = tri2d_edge_closed_form_data(tg);

    double scalar_sel[3][3];
    double scalar_tel[3][3];
    tri2d_scalar_closed_form_eq_30_33(tg.g.A, tg.g.b, tg.g.c, scalar_sel, scalar_tel);

    for (int m = 0; m < 3; ++m)
    {
        for (int n = 0; n < 3; ++n)
        {
            // Eq. (137): Sel(tt) = (1/mu) curl-curl - k0^2 eps_r massa_t.
        
            blk.Sel_tt[m][n] = inv_mu * tri2d_edge_curlcurl_entry_eq_66(edge, m, n) - k02 * eps_r * tri2d_edge_mass_entry_eq_67(edge, m, n);

            // Eq. (141): Tel(tt) = (1/mu) * massa_t.
            blk.Tel_tt[m][n] = inv_mu * tri2d_edge_mass_entry_eq_67(edge, m, n);
        }
    }

    for (int m = 0; m < 3; ++m)
    {
        for (int j = 0; j < 3; ++j)
        {
            const double c_mj = tri2d_edge_scalar_coupling_entry_base(edge, m, j);

            // Eq. (138): Tel(tz).
            blk.Tel_tz[m][j] = inv_mu * c_mj;

            // Eq. (139): Tel(zt) = Tel(tz)^T na notacao local.
            blk.Tel_zt[j][m] = inv_mu * c_mj;
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            // Eq. (140): Sel(zz) = (1/mu) grad-grad escalar.
            blk.Sel_zz[i][j] = inv_mu * scalar_sel[i][j];

            // Eq. (142): Tel(zz) = (1/mu) grad-grad - k0^2 eps_r massa_z.
            blk.Tel_zz[i][j] =
                inv_mu * scalar_sel[i][j] -
                k02 * eps_r * scalar_tel[i][j];
        }
    }

    return blk;
}

/******************************************************************************/
/* FUNCAO: tri2d_beta_rearranged_closed_form_eq_136                           */
/* DESCRICAO: Monta os sub-blocos locais no formato rearranjado usado pelo    */
/* codigo para a Secao 2.2.4, isto e, P x = beta^2 Q x. A funcao reaproveita  */
/* as Eq. (137) a (142) quando isso e algebricamente direto e preserva a      */
/* forma validada do repositorio para o termo Q_tt = -(1/mu_r) M_t.           */
/* ENTRADA: tg: const TriGeomEdge &; k0: double; eps_r: double; mu_r: double. */
/* SAIDA: Tri2DCoupledBetaRearrangedLocalBlocks.                              */
/******************************************************************************/
inline Tri2DCoupledBetaRearrangedLocalBlocks tri2d_beta_rearranged_closed_form_eq_136(
    const TriGeomEdge &tg,
    double k0,
    double eps_r,
    double mu_r)
{
    Tri2DCoupledBetaRearrangedLocalBlocks blk136;
    const auto blk = tri2d_beta_closed_form_eq_137_142(tg, k0, eps_r, mu_r);

    for (int m = 0; m < 3; ++m)
    {
        for (int n = 0; n < 3; ++n)
        {
            // Eq. (137), ja no lado esquerdo do sistema rearranjado.
            blk136.P_tt[m][n] = blk.Sel_tt[m][n];

            // Forma validada do repositorio para a Eq. (136):
            // Q_tt = -(1/mu_r) * M_t
            blk136.Q_tt[m][n] = -blk.Tel_tt[m][n];
        }
    }

    for (int m = 0; m < 3; ++m)
    {
        for (int j = 0; j < 3; ++j)
        {
            blk136.Q_tz[m][j] = blk.Tel_tz[m][j];
            blk136.Q_zt[j][m] = blk.Tel_zt[j][m];
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            // No rearranjo validado:
            //   Q_zz <- Sel(zz) = (1/mu_r) grad-grad
            blk136.Q_zz[i][j] = blk.Sel_zz[i][j];

            // Isola apenas a parcela k0^2 eps_r M_z para P_zz:
            //   Sel(zz) - Tel(zz) = k0^2 eps_r M_z
            blk136.P_zz[i][j] = blk.Sel_zz[i][j] - blk.Tel_zz[i][j];
        }
    }

    return blk136;
}

} // namespace explicit_tri2d
