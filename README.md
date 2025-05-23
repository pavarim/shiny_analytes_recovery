# Interface Shiny para Simulação e Recuperação de Curvas com Splines

Este projeto implementa em Shiny e R, um simulador e modelador funcional para curvas agregadas geradas a partir de concentrações aleatórias e curvas constituintes predefinidas (nesse caso, das funções Logito e Spahet).

---
## Motivação

Baseado no modelo de dados funcionais linear (especificamente a Lei de Beer-Lambert) com expansão em bases B-splines, o objetivo é estimar curvas constituintes (e.g, curvas de absorbância ou transmitância) a partir de curvas agregadas simuladas com ruído, controlável e seguindo uma combinação linear entre as curvas componentes e um erro aleatório.
O artigo "Análise Não-Paramétrica de Dados Funcionais: Uma Aplicação à Quimiometria" de Marley Apolinario Saraiva (2009) serviu de base e inspiração para esse projeto.

---
##  Funcionalidades

- Recuperação de curvas constituintes (logito, spahet)
- Manipulação de parâmetros como grau da Spline e número de Knots em tempo real para estudar casos de ajuste (como underfitting ou overfitting)
- Controle do ruído das curvas agregadas simuladas
- Simulação de curvas agregas (controlando quantidade, tamanho da curva e ruído)
- Visualização dos erros de estimação (EQM) nas recuperações
- Interação gráfica do usuário via Shiny
- Botão para download das estimações
---
## Exemplo de estimação
Um exemplo de recuperação das curvas abaixo:
![Exemplo](exemplo.png)

## Tecnologias
- Linguagem R
- ShinyApp
- pacote fda (functional data analysis)
- pacote dplyr
---
