/*****************************************************************************/
/* PROJETO: TP3485 FEM Eigen EM                                               */
/*****************************************************************************/
/* Arquivo: src/helm10/main_field_reconstruction_example.cpp                 */
/* Autor: Prof. Lucas Kriesel Sperotto                                        */
/* E-mail: speroto@unemat.br                                                  */
/* Versao: 2.0 | Ano: 2026                                                    */
/*****************************************************************************/
/* Descricao: Exemplo minimo e didatico do uso da reconstrucao de campos      */
/* transversais do HELM10, sem depender do solve modal completo.              */
/*****************************************************************************/
/* Artigo base: NASA 19950011772. Referencias principais: Eq. (35)-(41).     */
/*****************************************************************************/

#include "helm10/field_reconstruction.hpp"

#include <iomanip>
#include <iostream>
#include <vector>

/******************************************************************************/
/* FUNCAO: main                                                               */
/* DESCRICAO: Demonstra, em uma malha minima, como reconstruir gradientes e   */
/* campos transversais TE/TM a partir de um potencial longitudinal nodal.     */
/* ENTRADA: argc: int; argv: char **.                                         */
/* SAIDA: int.                                                                */
/******************************************************************************/
int main()
{
    using helm10::field_reconstruction::LongitudinalScalarKind;
    using helm10::field_reconstruction::make_homogeneous_medium_above_kc;
    using helm10::field_reconstruction::reconstruct_transverse_fields;

    // Malha minima com dois triangulos. A ideia aqui e puramente didatica:
    // queremos um dominio simples no qual o potencial escalar varie ao longo
    // do plano, para que dpsi/dx e dpsi/dy nao sejam triviais.
    Mesh2D mesh;
    mesh.nodes = {
        {0.0, 0.0, true},
        {1.0, 0.0, true},
        {1.0, 1.0, true},
        {0.0, 1.0, true},
    };
    mesh.tris = {
        {{0, 1, 2}},
        {{0, 2, 3}},
    };

    // Potencial longitudinal nodal escolhido apenas para ilustrar a passagem
    // matematica das Eq. (35)-(41).
    const std::vector<double> psi = {0.0, 1.0, 1.5, 0.5};

    // Politica padrao do HELM10: se a frequencia nao for informada, escolhemos
    // um valor automatico acima do maior kc que sera exportado, para manter os
    // modos propagantes na interpretacao fisica.
    const auto medium = make_homogeneous_medium_above_kc(3.0, 1.0, 1.0);

    // Exemplo TE: a variavel escalar e associada a Hz.
    const auto te = reconstruct_transverse_fields(
        mesh,
        psi,
        LongitudinalScalarKind::TE_Hz,
        3.0,
        medium,
        false);

    // Exemplo TM: a variavel escalar e associada a Ez. Mesmo com frequencia
    // disponivel, o padrao atual do projeto salva Ex e Ey sem multiplicar por
    // Ztm. Beta e Ztm ficam como informacao fisica auxiliar.
    const auto tm = reconstruct_transverse_fields(
        mesh,
        psi,
        LongitudinalScalarKind::TM_Ez,
        3.0,
        medium,
        false);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "[field reconstruction example]\n";
    std::cout << "freq_hz=" << medium.frequency_hz
              << " k=" << medium.k << "\n";
    std::cout << "TE status: " << te.status_label << " | " << te.status_message << "\n";
    std::cout << "TM status: " << tm.status_label << " | " << tm.status_message << "\n";
    std::cout << "TM beta(info)=" << tm.beta
              << " Ztm(info)=" << tm.ztm << "\n";

    std::cout << "\nnode  x  y   psi   dpsi_dx   dpsi_dy   Ex_TE   Ey_TE   Ex_TM_semZ   Ey_TM_semZ   Ex_TM_comZ   Ey_TM_comZ\n";
    for (size_t node_id = 0; node_id < mesh.nodes.size(); ++node_id)
    {
        const auto &node = mesh.nodes[node_id];
        std::cout << node_id << "  "
                  << node.x << "  "
                  << node.y << "  "
                  << psi[node_id] << "  "
                  << te.dpsi_dx[node_id] << "  "
                  << te.dpsi_dy[node_id] << "  "
                  << te.ex[node_id] << "  "
                  << te.ey[node_id] << "  "
                  << tm.ex[node_id] << "  "
                  << tm.ey[node_id] << "  "
                  << tm.ex_with_ztm[node_id] << "  "
                  << tm.ey_with_ztm[node_id] << "\n";
    }

    return 0;
}
