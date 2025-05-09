#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(fda)


# Define UI for application that draws a histogram
ui <- fluidPage(
  
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("snr",
                        "Razão Sinal-Ruído:",
                        min = 1,
                        max = 20,
                        value = 3),
            textInput("ivqtd", "Quantidade de pontos:", value = 50),
            textInput("samplen", "Número de amostras:", value = 100),
            textInput("splinedegree", "Grau da spline:", value = 3),
            textInput("knotsn", "Número de Knots:", value = 4),
            h5("Funções constituintes", align="center"),
            plotOutput("analytesplot")
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("curve"),
           plotOutput("esm"),
           downloadButton("baixarPlot", "Baixar gráfico das estimações")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
    ids = reactive({
      req(input$ivqtd)
      
      seq(0,1, length = as.numeric(input$ivqtd))
    })
    
    logit = reactive({ 1 / (1 + exp(-20*(ids()-0.5))) })
    spahet = reactive({sqrt(ids()*(1-ids()))*sin(2*pi*(1+2^(-0.6)) /(ids() + 2^(-0.6)))})
    
    output$analytesplot <- renderPlot({
        par(mfrow = c(2,2))
        plot(ids(),logit(),t="l", ylab = "Logito", xlab = "Id")
        plot(ids(), spahet(),t = "l", ylab = "Spahet", xlab = "Id")
    })
    
    concentracao = reactive({
      req(input$samplen)
      data.frame(y1 = runif(as.numeric(input$samplen)), y2 = runif(as.numeric(input$samplen))) %>% unname() %>% as.matrix()
    }) 
    
    absorbancia = reactive({
      req(input$snr)
      req(input$samplen)
      
      agreg = 0.6*logit() + 0.4*spahet()
      sd_e = sd(agreg) / as.numeric(input$snr)
      abs = matrix(0,nrow = as.numeric(input$samplen), ncol = length(ids()))
      for (i in 1:nrow(concentracao())) {
        erro_i = rnorm(length(ids()), mean = 0, sd = sd_e)
        abs[i,] = concentracao()[i,1]*logit() + concentracao()[i,2]*spahet() + erro_i
      }
      abs
    })
    
    output$curve = renderPlot({
      plot(ids(), absorbancia()[1,], type = "o", ylim = c(min(absorbancia()), max(absorbancia())),
           xlab = "Id", ylab = "Curva agregada", main = paste0("Simulação de ", input$samplen, " curvas agregadas"))
      for (i in 2:dim(absorbancia())[1]){
        lines(ids(), absorbancia()[i,], col = i, type="o")
      }
    })
    
    B = reactive({
      req(input$splinedegree)
      req(input$knotsn)
      req(input$ivqtd)
      req(input$samplen)
      
      k_n = as.numeric(input$knotsn)
      
      validate(
        need(k_n > as.numeric(input$splinedegree), "Número de knots deve ser maior que o grau do polinômio"),
        need(k_n <= as.numeric(input$ivqtd) - 1, "Número de intervalos não pode ultrapassar a quantidade de pontos"),
        need(as.numeric(input$samplen) > 2, "Tamanho amostral muito pequeno")
      )
      
      knots = seq(0,1,length.out = k_n)
      absbasis = create.bspline.basis(rangeval = c(0,1),norder = as.numeric(input$splinedegree), breaks = knots)
      return(getbasismatrix(ids(), basisobj = absbasis, nderiv = 0))
    })
    
    X = reactive({
      as.vector(t(absorbancia()))
    })
    
    D = reactive({
      req(input$samplen)
        n <- as.numeric(input$samplen)
        o <- length(ids())
        m <- dim(concentracao())[2]
        l <- ncol(B())
        
        d <- matrix(0, n * o, l * (m + 1)) 
        
        row_index <- 1
        for (i in 1:n) {
          y_row <- c(1, as.numeric(concentracao()[i,]))
          for (p in 1:o) {
            b_row <- B()[p, ]
            d_row <- c()
            for (j in 1:(m + 1)) {
              d_row <- c(d_row, b_row * y_row[j])
            }
            d[row_index, ] <- d_row
            row_index <- row_index + 1
          }
        }
        return(d)
    })
    
    Betahat = reactive({
      solve(t(D())%*%D())%*%t(D())%*%X()
      
  })
    
    intercepto = reactive({
      l <- ncol(B())
      B() %*% Betahat()[1:l]
    })
    
    logit_esm = reactive({
      l <- ncol(B())
      B() %*% Betahat()[(l+1):(2*l)]
    })
    
    observe({
      logit_esm()
    })
    
    spahet_esm = reactive({
      l <- ncol(B())
      B() %*% Betahat()[(2*l+1):(3*l)] 
    })
    
    output$esm = renderPlot({
      r2l = 1 - (sum(((logit() - logit_esm())^2)) / sum((logit()-mean(logit()))^2)) 
      r2s = 1 - (sum(((spahet() - spahet_esm())^2)) / sum((spahet()-mean(spahet()))^2)) 
      eqml = sum((logit_esm() - logit())^2)/ncol(absorbancia()) 
      eqms = sum((spahet_esm() - spahet())^2)/ncol(absorbancia()) 

      
      par(mfrow = c(1,2))
      plot(ids(), logit(), t = "l", ylab = "Logito x Logito estimada", xlab = "Id",
           main = paste0("EQM: ",format(eqml,digits=4,scientific=FALSE)))
      lines(ids(), logit_esm(), col = "red")
      plot(ids(), spahet(), t = "l", ylab = "Spahet x Spahet estimada", xlab = "Id",
           main = paste0("EQM: ", format(eqms,digits=4,scientific=FALSE)))
      lines(ids(), spahet_esm(), col = "red")
    })
    
    output$baixarPlot <- downloadHandler(
      filename = "estimativas.png",
      content = function(file) {
        eqml = sum((logit_esm() - logit())^2)/ncol(absorbancia()) 
        eqms = sum((spahet_esm() - spahet())^2)/ncol(absorbancia()) 
        png(file, width = 2000, height = 1000, res = 300)  # Alta resolução
        par(mfrow = c(1,2))
        plot(ids(), logit(), t = "l", ylab = "Logito x Logito estimada", xlab = "Id",
             main = paste0("EQM: ",format(eqml,digits=4,scientific=FALSE)))
        lines(ids(), logit_esm(), col = "red")
        plot(ids(), spahet(), t = "l", ylab = "Spahet x Spahet estimada", xlab = "Id",
             main = paste0("EQM: ", format(eqms,digits=4,scientific=FALSE)))
        lines(ids(), spahet_esm(), col = "red")
        dev.off()
      }
    )
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
