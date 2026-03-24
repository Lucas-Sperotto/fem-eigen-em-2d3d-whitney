# Handoff para outro repositorio: guias dieletricos de canal retangulares

Este arquivo resume o que deve ser enviado para outro chat/repositorio para
implementar, em uma base nova, uma extensao das Secoes 2.2.3 e 2.2.4 do paper
NASA TP-3485 para guias dieletricos de canal retangulares com tres indices
`n1`, `n2` e `n3`.

## 1) O que enviar para o outro chat

Envie:

1. O artigo novo com:
   - geometria exata do guia de canal,
   - definicao de `n1`, `n2`, `n3`,
   - equacoes fracas,
   - condicoes de contorno,
   - tabelas/figuras de validacao.

2. Estes arquivos deste repositorio, como referencia de montagem acoplada:
   - `src/helmvec2/README.md`
   - `src/helmvec3/README.md`
   - `src/helmvec2/helmvec2_coupled_system.hpp`
   - `src/helmvec2/helmvec2_coupled_system.cpp`
   - `src/helmvec2/helmvec23_shared.hpp`
   - `src/helmvec2/main_helmvec2_rect.cpp`
   - `src/helmvec3/main_helmvec3_rect.cpp`

3. Estes arquivos de infraestrutura FEM 2D que a montagem acoplada usa:
   - `src/edge/edge_assembly.hpp`
   - `src/edge/edge_assembly.cpp`
   - `src/edge/edge_basis.hpp`
   - `src/edge/edge_basis.cpp`
   - `src/edge/edge_dofs.hpp`
   - `src/edge/edge_dofs.cpp`
   - `src/core/helm10_scalar_system.hpp`
   - `src/core/helm10_scalar_system.cpp`
   - `src/core/fem_scalar.hpp`
   - `src/core/fem_scalar_helm10.cpp`
   - `src/core/mesh2d.hpp`
   - `src/core/mesh2d_rect.hpp`
   - `src/core/mesh2d_rect.cpp`
   - `src/core/dense.hpp`
   - `src/core/block_ops.hpp`
   - `src/core/lapack_eig.hpp`

4. Opcional, mas util:
   - `README.md`
   - capturas das tabelas do artigo novo,
   - qualquer resultado de referencia que voce ja tenha.

## 2) O que o outro chat precisa entender antes de codificar

O modelo novo precisa confirmar, a partir do artigo novo:

1. Se o problema continua sendo 2D vetorial+escalar no plano transversal
   com incognitas `Et/Ez`.

2. Se `eps_r(x,y) = n(x,y)^2` e `mu_r = 1` em todo o dominio.

3. Se o guia e fechado por PEC em um dominio truncado, ou se e um guia
   dieletrico aberto. Isto e critico:
   - no repositorio atual, as Secoes 2.2.3 e 2.2.4 assumem contorno PEC;
   - em guia dieletrico aberto, o tratamento de contorno pode mudar
     completamente (caixa truncada, contorno absorvente, PML, etc.).

4. Se os tres indices `n1`, `n2`, `n3` representam:
   - nucleo / substrato / cobertura,
   - tres regioes laterais,
   - ou outra decomposicao geometrica.

5. Qual e o observavel do artigo novo:
   - `k0` para `beta` dado,
   - `beta` para `k0` dado,
   - `n_eff = beta / k0`,
   - curvas de dispersao,
   - ou todos.

6. Quais casos base devem ser reproduzidos e em que formato:
   - tabelas,
   - curvas,
   - perfis de campo,
   - erros percentuais.

## 3) Prompt pronto para o outro chat

Use o texto abaixo, junto com o artigo novo e os arquivos listados acima.

```text
Quero implementar neste repositorio uma extensao das Secoes 2.2.3 e 2.2.4 do
meu repositorio anterior para guias dieletricos de canal retangulares com tres
indices `n1`, `n2` e `n3`.

Contexto:
- Estou anexando o artigo base novo com geometria, equacoes e casos de teste.
- Tambem estou anexando arquivos de referencia do repositorio anterior que
  implementam os problemas acoplados `k0(beta dado)` e `beta(k0 dado)` para
  guia retangular parcialmente preenchido em PEC.
- Quero reutilizar ao maximo a arquitetura existente, mas adaptando a
  formulacao e as condicoes de contorno ao artigo novo.

Arquivos de referencia anexados:
- `src/helmvec2/README.md`
- `src/helmvec3/README.md`
- `src/helmvec2/helmvec2_coupled_system.hpp`
- `src/helmvec2/helmvec2_coupled_system.cpp`
- `src/helmvec2/helmvec23_shared.hpp`
- `src/helmvec2/main_helmvec2_rect.cpp`
- `src/helmvec3/main_helmvec3_rect.cpp`
- `src/edge/edge_assembly.hpp`
- `src/edge/edge_assembly.cpp`
- `src/edge/edge_basis.hpp`
- `src/edge/edge_basis.cpp`
- `src/edge/edge_dofs.hpp`
- `src/edge/edge_dofs.cpp`
- `src/core/helm10_scalar_system.hpp`
- `src/core/helm10_scalar_system.cpp`
- `src/core/fem_scalar.hpp`
- `src/core/fem_scalar_helm10.cpp`
- `src/core/mesh2d.hpp`
- `src/core/mesh2d_rect.hpp`
- `src/core/mesh2d_rect.cpp`
- `src/core/dense.hpp`
- `src/core/block_ops.hpp`
- `src/core/lapack_eig.hpp`

Quero que voce trabalhe assim:

1. Leia primeiro o artigo novo e os arquivos anexados.
2. Resuma, antes de editar codigo, estas decisoes:
   - geometria transversal exata;
   - significado de `n1`, `n2`, `n3`;
   - relacao entre indice e permissividade relativa;
   - incognitas do problema;
   - forma fraca;
   - tipo de problema de autovalor;
   - condicoes de contorno;
   - casos de validacao.
3. Identifique claramente o que pode ser reaproveitado do repositorio anterior
   e o que precisa ser refeito.
4. Implemente primeiro a infraestrutura geometrica/material:
   - perfil `eps_r(x,y)` por regioes do canal;
   - utilitarios para `n1`, `n2`, `n3`;
   - geracao de casos de teste.
5. Depois implemente:
   - o problema equivalente a 2.2.3 (`k0` dado `beta`);
   - o problema equivalente a 2.2.4 (`beta` dado `k0`);
   - pos-processamento fisico (`beta`, `k0`, `n_eff`, filtragem de raizes).
6. Gere saidas reprodutiveis:
   - tabela CSV,
   - resumo em Markdown,
   - se o artigo pedir, curvas/figuras.
7. Comente tudo em portugues e mantenha o codigo didatico.
8. Reutilize a estrutura das funcoes do repositorio anterior sempre que fizer
   sentido, mas sem copiar cegamente as BCs de PEC se o artigo novo nao usar
   PEC.

Requisitos tecnicos:
- assumir `mu_r = 1` apenas se o artigo permitir;
- usar `eps_r = n^2` apenas se o artigo confirmar isso;
- manter separacao entre montagem edge, montagem escalar e sistema acoplado;
- deixar claro no codigo o equivalente matematico dos blocos `A/B` e `P/Q`;
- validar contra as tabelas/curvas do artigo novo com erro relativo percentual.

Entrega esperada:
- codigo completo no novo repositorio;
- comentarios em portugues;
- README do modulo novo;
- comandos de execucao;
- comparacao numerica com o artigo.

Antes de comecar a editar, quero um resumo tecnico do artigo novo e um plano de
implementacao em etapas.
```

## 4) Instrucoes adicionais para voce enviar junto

Se quiser mandar um complemento curto no chat novo, use isto:

```text
Importante: nao assuma que o problema novo usa o mesmo contorno PEC do
repositorio anterior. O objetivo e reaproveitar a arquitetura de montagem
acoplada de 2.2.3/2.2.4, nao forcar a mesma fisica sem verificar o artigo.
```

## 5) Estrategia recomendada de portabilidade

No repositorio novo, a ordem recomendada e:

1. copiar apenas a infraestrutura minima de algebra e malha;
2. portar `edge_assembly` e `helm10_scalar_system`;
3. portar `helmvec2_coupled_system`;
4. substituir os perfis `eps_step_x/eps_step_y` por um perfil generico de
   canal com `n1`, `n2`, `n3`;
5. criar drivers novos especificos do artigo;
6. adicionar validacao e exportacao.

## 6) Observacao importante

Se o artigo novo tratar guia dieletrico aberto, o outro chat deve ser
pressionado a explicar claramente:

- como o dominio aberto sera truncado;
- como a radiacao sera tratada;
- por que a BC escolhida e consistente com o artigo.

Sem isso, existe alto risco de "reusar" 2.2.3/2.2.4 de forma estruturalmente
parecida, mas fisicamente errada.
