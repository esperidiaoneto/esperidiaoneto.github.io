# Cladograma simples a partir de matriz de presença (1) / ausência (0)
# Os números sobre os ramos indicam em qual ramo cada caractere surge.
# Instalar uma vez: install.packages(c("shiny", "ape", "phangorn", "readxl"))
#
# Formatos aceitos: .xlsx, .xls, .csv (vírgula ou ponto e vírgula), .tsv, .txt
# (separado por tabulação) e .nex/.nexus. Primeira coluna = nome dos táxons.
# Estados: 0 e 1; "?" e "-" são tratados como dado faltante.

library(shiny)
library(ape)
library(phangorn)
library(readxl)

COR_SINAPOMORFIA <- "#1f6fb2"   # azul: ramo interno
COR_AUTAPOMORFIA <- "#d9534f"   # vermelho: ramo terminal

# ------------------------------------------------------------------------------
# Leitura de tabelas em vários formatos
# ------------------------------------------------------------------------------
ler_tabela <- function(caminho, nome_arquivo) {
  ext <- tolower(tools::file_ext(nome_arquivo))
  
  if (ext %in% c("xlsx", "xls")) {
    df <- as.data.frame(read_excel(caminho, col_types = "text"),
                        stringsAsFactors = FALSE)
    nomes <- df[[1]]
    m <- as.matrix(df[, -1, drop = FALSE])
    
  } else if (ext %in% c("nex", "nexus")) {
    m <- do.call(rbind, read.nexus.data(caminho))
    nomes <- rownames(m)
    if (is.null(colnames(m))) colnames(m) <- paste("Caractere", seq_len(ncol(m)))
    
  } else {  # csv, tsv, txt
    linha1 <- readLines(caminho, n = 1, warn = FALSE)
    sep <- if (grepl("\t", linha1)) "\t" else if (grepl(";", linha1)) ";" else ","
    ler <- function(enc) read.csv(caminho, sep = sep, colClasses = "character",
                                  check.names = FALSE, na.strings = c("", "NA"),
                                  fileEncoding = enc, stringsAsFactors = FALSE)
    df <- tryCatch(ler("UTF-8"),
                   error = function(e) ler("latin1"),
                   warning = function(w) ler("latin1"))
    nomes <- df[[1]]
    m <- as.matrix(df[, -1, drop = FALSE])
  }
  
  nomes <- trimws(as.character(nomes))
  m[] <- trimws(m)
  
  # Mantém só as linhas de dados (descarta, p. ex., "Classificação Evolutiva")
  valida <- apply(m, 1, function(r) {
    r <- r[!is.na(r)]
    length(r) > 0 && all(r %in% c("0", "1", "?", "-"))
  })
  valida <- valida & !is.na(nomes) & nomes != ""
  m <- m[valida, , drop = FALSE]
  nomes <- nomes[valida]
  
  m[is.na(m)] <- "?"
  m[m == "-"] <- "?"
  
  sugerido <- grepl("outgroup|grupo externo|externo", nomes, ignore.case = TRUE)
  rownames(m) <- make.unique(trimws(sub("\\s*\\(.*\\)", "", nomes)))  # nomes curtos
  
  list(matriz = m, outgroup = rownames(m)[sugerido])
}

# ------------------------------------------------------------------------------
# Desenho do cladograma com os caracteres marcados nos ramos
# ------------------------------------------------------------------------------
desenhar <- function(arv, mk, mostrar) {
  marcar <- mostrar && length(mk$arestas) > 0
  cor_marca <- ifelse(mk$terminal, COR_AUTAPOMORFIA, COR_SINAPOMORFIA)
  
  cores_ramo <- rep("gray25", nrow(arv$edge))
  if (marcar) cores_ramo[mk$arestas] <- cor_marca
  
  op <- par(mar = c(4, 1, 1, 1))
  on.exit(par(op))
  
  plot(arv, type = "cladogram", use.edge.length = FALSE, cex = 1.3,
       edge.width = 3, edge.color = cores_ramo, label.offset = 0.05)
  
  if (marcar) {
    edgelabels(text = paste0(" ", mk$rotulos, " "), edge = mk$arestas,
               frame = "rect", bg = cor_marca, col = "white", font = 2, cex = 1)
    legend("bottom", inset = c(0, -0.1), xpd = TRUE, horiz = TRUE, bty = "n",
           fill = c(COR_SINAPOMORFIA, COR_AUTAPOMORFIA),
           legend = c("Sinapomorfia (ramo interno)", "Autapomorfia (ramo terminal)"))
  }
}

# ------------------------------------------------------------------------------
# Interface
# ------------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Cladograma a partir de matriz de caracteres"),
  sidebarLayout(
    sidebarPanel(
      fileInput("arquivo", "Matriz (.xlsx, .xls, .csv, .tsv, .txt, .nex)",
                accept = c(".xlsx", ".xls", ".csv", ".tsv", ".txt", ".nex", ".nexus"),
                buttonLabel = "Escolher...", placeholder = "Nenhum arquivo"),
      uiOutput("grupo_externo"),
      checkboxInput("mostrar", "Mostrar caracteres nos ramos", TRUE),
      downloadButton("baixar", "Baixar figura (PNG)"),
      h4("Caracteres"),
      verbatimTextOutput("legenda")
    ),
    mainPanel(
      plotOutput("cladograma", height = "550px"),
      textOutput("aviso")
    )
  )
)

# ------------------------------------------------------------------------------
# Servidor
# ------------------------------------------------------------------------------
server <- function(input, output) {
  
  lido <- reactive({
    req(input$arquivo)
    r <- tryCatch(ler_tabela(input$arquivo$datapath, input$arquivo$name),
                  error = function(e) NULL)
    validate(need(!is.null(r) && nrow(r$matriz) >= 4,
                  "Não consegui ler a matriz (são necessários ao menos 4 táxons com dados 0/1)."))
    r
  })
  
  matriz <- reactive(lido()$matriz)
  
  output$grupo_externo <- renderUI({
    taxa <- rownames(matriz())
    padrao <- if (length(lido()$outgroup) > 0) lido()$outgroup[1] else taxa[1]
    selectInput("og", "Grupo externo", choices = taxa, selected = padrao)
  })
  
  output$legenda <- renderText({
    paste(seq_len(ncol(matriz())), colnames(matriz()), sep = " = ", collapse = "\n")
  })
  
  # Árvore de máxima parcimônia (consenso estrito), enraizada no grupo externo
  arvore <- reactive({
    req(input$og)
    dados <- phyDat(matriz(), type = "USER", levels = c("0", "1"), ambiguity = "?")
    arvores <- bab(dados, trace = 0)
    arv <- if (inherits(arvores, "phylo")) arvores else consensus(arvores, p = 1)
    arv <- root(arv, outgroup = input$og, resolve.root = TRUE)
    ladderize(arv)
  })
  
  # Para cada caractere, procura o ramo que leva aos táxons com "1"
  # (assume 0 = estado ancestral e 1 = estado derivado)
  marcas <- reactive({
    arv <- arvore()
    m <- matriz()
    pp <- prop.part(arv)
    rotulos <- attr(pp, "labels")
    clados <- c(as.list(arv$tip.label), lapply(pp, function(i) rotulos[i]))
    tamanho <- lengths(clados)
    
    por_ramo <- list()
    sem_ramo <- c()
    for (j in seq_len(ncol(m))) {
      com_1 <- rownames(m)[which(m[, j] == "1")]
      com_0 <- rownames(m)[which(m[, j] == "0")]
      aresta <- NA
      if (length(com_1) > 0) {
        no <- which(sapply(clados, function(x) all(com_1 %in% x) && !any(com_0 %in% x)))
        no <- no[order(tamanho[no])]
        if (length(no) > 0) aresta <- match(no[1], arv$edge[, 2])
      }
      if (is.na(aresta)) {
        sem_ramo <- c(sem_ramo, j)
      } else {
        chave <- as.character(aresta)
        por_ramo[[chave]] <- c(por_ramo[[chave]], j)
      }
    }
    
    arestas <- as.integer(names(por_ramo))
    list(arestas = arestas,
         rotulos = unname(unlist(lapply(por_ramo, paste, collapse = ","))),
         terminal = arv$edge[arestas, 2] <= Ntip(arv),
         sem_ramo = sem_ramo)
  })
  
  output$cladograma <- renderPlot({
    desenhar(arvore(), marcas(), input$mostrar)
  })
  
  output$aviso <- renderText({
    sem <- marcas()$sem_ramo
    if (length(sem) > 0) {
      paste("Caracteres que não definem um único ramo nesta árvore:",
            paste(sem, collapse = ", "))
    } else ""
  })
  
  output$baixar <- downloadHandler(
    filename = function() "cladograma.png",
    content = function(file) {
      png(file, width = 2400, height = 1800, res = 300)
      on.exit(dev.off())
      desenhar(arvore(), marcas(), input$mostrar)
    }
  )
}

shinyApp(ui, server)