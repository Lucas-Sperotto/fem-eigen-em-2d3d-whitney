# 1. Introdução

Uma boa introdução não serve apenas para "abrir" um artigo. Ela deve ensinar ao leitor qual problema merece ser levado a sério, por que certos caminhos falharam no passado e qual será a trilha de ideias que tornará o estudo fecundo. Esta introdução cumpre exatamente esse papel.

## Texto-base traduzido

O método dos elementos finitos (FEM) tem sido amplamente utilizado como ferramenta de análise e projeto em muitas disciplinas de engenharia, como estruturas e dinâmica computacional. Embora o FEM tenha sido aplicado ao eletromagnetismo, permaneceu principalmente confinado a máquinas elétricas e ímãs (ref. 1). Nos últimos 20 anos, tem havido grande interesse na aplicação deste método a componentes de micro-ondas, como guias de onda e antenas. Porém, por muitos anos, seu uso foi restrito devido às chamadas soluções espúrias associadas aos elementos de aresta (ref. 2).

Muito recentemente, foram empregados os chamados “elementos de aresta” com sucesso para superar essas soluções espúrias. Nos últimos anos, o uso desses elementos no FEM reavivou o interesse em aplicar o FEM a problemas de engenharia de micro-ondas (ref. 3). Isso, em combinação com os avanços em hardware e software computacional, ajudou a tornar o FEM uma ferramenta atrativa para eletromagnetismo. Além disso, existe uma variedade de ferramentas comerciais de modelagem geométrica para modelar com precisão qualquer geometria tridimensional e gerar a malha necessária com elementos do tipo triângulos e tetraedros (refs. 4 e 5).

Neste artigo, as ferramentas de FEM para análise de problemas de autovalores em eletromagnetismo foram descritas. Este artigo é dividido em duas partes: a seção 2 trata dos problemas bidimensionais; a seção 3, dos problemas tridimensionais. Ao longo deste artigo, elementos triangulares são utilizados para modelar problemas bidimensionais e elementos tetraédricos são utilizados para modelar problemas tridimensionais.

Na seção 2.1, uma formulação FEM escalar é utilizada para guias de onda bidimensionais de forma arbitrária. Elementos triangulares com funções de base nodais são utilizados para formular as matrizes do FEM. Os autovalores para diferentes tipos de guias de onda são obtidos e os gráficos de intensidade de campo são apresentados para vários modos de guia de onda.

Na seção 2.2, um FEM vetorial é introduzido com elementos de aresta bidimensionais para analisar guias de onda não homogêneos. Para maior clareza na formulação, a seção 2.2 é dividida em quatro subseções. A seção 2.2.1 fornece a solução de guias de onda homogêneos com campos vetoriais transversais de duas componentes. A seção 2.2.2 apresenta o cálculo de autovalores para guias de onda não homogêneos utilizando campos vetoriais de três componentes. Combinações de funções de base de aresta e nodais foram utilizadas para campos vetoriais e componentes de campo longitudinal, respectivamente. As seções 2.2.3 e 2.2.4 estendem a formulação da seção 2.2.2 para determinar o número de onda ou a constante de propagação para guias de onda preenchidos de forma não homogênea, quando um deles é especificado.

Na seção 3.1, é descrita a formulação do FEM vetorial tridimensional. Funções de base de aresta para elementos tetraédricos são introduzidas para formular as matrizes de elementos finitos para cavidades tridimensionais preenchidas com material não homogêneo.

Nas seções 2.1, 2.2 e 3.1, exemplos numéricos são apresentados para demonstrar a validade da análise e dos programas computacionais desenvolvidos. Para todos os exemplos, os resultados do FEM apresentam boa precisão. Todos os exemplos numéricos foram verificados quanto à convergência numérica. Em virtude do FEM, os códigos computacionais apresentados neste artigo podem lidar com qualquer geometria de forma arbitrária preenchida com materiais não homogêneos, salvo indicação em contrário.

## O problema histórico que a introdução realmente anuncia

Se eu estivesse explicando esta abertura a um aluno muitos anos depois da publicação do artigo, eu diria o seguinte: o tema central aqui não é apenas "usar elementos finitos em eletromagnetismo". Isso já seria pouco. O verdadeiro problema é usar o FEM de modo fiel à física do campo eletromagnético, sem introduzir soluções matematicamente admissíveis, mas fisicamente falsas.

É por isso que a questão dos modos espúrios domina silenciosamente a introdução inteira. Toda a arquitetura posterior do artigo nasce dessa ferida conceitual.

## A virada intelectual do texto

A grande mudança anunciada nesta introdução é a passagem de formulações vetoriais inadequadas para formulações apoiadas em elementos de aresta. Em linguagem mais simples: os autores perceberam que não basta discretizar a equação; é preciso discretizar também o tipo correto de continuidade do campo.

Esse é um ponto que merece ser guardado com cuidado:

- em problemas eletromagnéticos, o campo vetorial não "vive" no espaço discreto de qualquer interpolação nodal;
- a continuidade fisicamente relevante para o campo elétrico tangencial é diferente da continuidade natural de um campo escalar;
- por isso, a escolha do espaço discreto é parte da física do problema, e não mero detalhe de programação.

## Como um leitor sério deve interpretar a estrutura anunciada

Esta introdução também oferece o mapa completo do artigo. Ela nos diz que o texto progride em três grandes degraus:

### 1. Primeiro degrau: o caso escalar 2D

O estudo começa no terreno mais transparente. Isso é pedagogicamente sábio. No problema escalar, a passagem da forma forte à forma fraca e depois à matriz global é mais fácil de enxergar.

### 2. Segundo degrau: o caso vetorial 2D

Aqui começa a verdadeira prova de maturidade da formulação. O artigo primeiro trata o caso vetorial homogêneo, depois o caso não homogêneo com três componentes, e por fim os problemas em que `k_0` ou `\beta` passam a ser a incógnita espectral.

### 3. Terceiro degrau: o caso 3D

O problema tridimensional não surge como ruptura, mas como extensão natural do raciocínio já consolidado em 2D. O estudante atento perceberá que o artigo não salta diretamente para 3D; ele prepara cuidadosamente o caminho.

## Interpretação didática da introdução

Há uma lição de método científico nesta página. Os autores não começam pelo problema mais impressionante; começam pelo problema mais instrutivo. Isso revela maturidade intelectual. Quem deseja deixar uma obra duradoura não tenta apenas exibir resultados: procura construir uma sequência de ideias que outras pessoas possam reaprender e reutilizar.

Neste sentido, a introdução já é quase um manifesto pedagógico. Ela diz, ainda que indiretamente:

- não comece pelo mais complicado se você ainda não entendeu o essencial;
- não confunda sofisticação algébrica com fidelidade física;
- não trate a geometria e o espaço discreto como acessórios secundários.

## Pontos de consistência conceitual

- A passagem que menciona soluções espúrias precisa ser lida com cuidado histórico: o problema clássico estava em formulações vetoriais inadequadas, em especial quando o espaço discreto não respeitava a estrutura do campo. Os elementos de aresta aparecem no texto como remédio para isso, não como causa.
- A escolha de triângulos e tetraedros não é apenas conveniência geométrica. Ela prepara o terreno para tratar domínios arbitrários sem sacrificar a lógica variacional do método.
- O artigo se organiza em torno de problemas de autovalor porque eles fornecem um campo de teste extremamente exigente para a formulação numérica.

## Conselhos de leitura para o restante do estudo

- Ao entrar na Seção 2.1, o leitor deve observar com calma como nasce a forma fraca.
- Ao entrar na Seção 2.2, deve prestar atenção especial à mudança do espaço discreto.
- Ao entrar na Seção 3.1, deve tentar reconhecer o que foi preservado da filosofia 2D e o que precisou ser ampliado.

## Fecho de leitura

Uma introdução como esta merece ser relida depois que todo o artigo tiver sido estudado. Só então se percebe que ela não era apenas uma apresentação: era uma promessa de método, uma declaração de princípios e um roteiro intelectual para todo o restante da obra.
