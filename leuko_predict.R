#############################
##B_ALL risk model predictor# 
##Marina Reixachs 2021#######
#############################

library(shiny)
library(shinythemes)
library(biomaRt)
library(DT)
library(stringr)
library(reshape2)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(randomForest)
library(survival)
library(survminer)
library(tidyverse)

# deploy with bioconductor packages
#library(BiocManager)
#options(repos = BiocManager::repositories())

#option to increase input table size
options(shiny.maxRequestSize=300*1024^2)


#FUNCTIONS ----
#plotly vertical line
vline <- function(x = 0, color = "red") {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color, width = 2)
  )
}

hline <- function(y = 0, color = "grey") {
  list(
    type = "line", 
    x0 = 0, 
    x1 = 1, 
    xref = "paper",
    y0 = y, 
    y1 = y, 
    line = list(color = color,width=1, dash='dash')
  )
}


# Load data and models -----

#Load model data
load("rf_model_target_v2.RData")

modelGenes = unique(rownames(leuko_rf_target$importance))


#Load TARGET data


#Gene gene name to ENSG conversion
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")

ens2gene <- getBM(attributes=c('external_gene_name','ensembl_gene_id', "chromosome_name"),
      mart = mart, filters = 'external_gene_name', values = modelGenes )
ens2gene = ens2gene[ens2gene$chromosome_name %in% c(1:21, "X", "Y"),]

#Generate palette for gene expression
n <- nrow(ens2gene)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' & !grepl("Pastel", rownames(brewer.pal.info)),]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = unique(col_vector[!grepl("FFFF", col_vector)])


#TARGET
target.lcpm = read.table("TARGET_lcpm_signature.txt", header = T, sep = "\t")
target.info = read.table("TARGET_clinical_info.txt", header = T, sep = "\t")


# Define UI  ---------

ui <- fluidPage(theme = shinytheme("sandstone"),
    navbarPage("B-ALL risk model",
    tabPanel("Main",
    
    ## Sidebar ----
    sidebarLayout(
        sidebarPanel(
          h3("Expression file (logCPMs):"),
          fileInput("upload","Choose a file to upload",buttonLabel = "Upload...",multiple = FALSE),
          tags$hr(),
          
          ### Input: Checkbox if file has header ----
          checkboxInput("header", "Header", TRUE),
          checkboxInput("asrownames", "Genes as rownames", FALSE),
        
          ### Input: Select separator ----
          radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = "\t"),
      
          ### Input: threshold slider ----
          h3("Score threshold:"),
          sliderInput("threshold", "",
                      min = 0, max = 1,
                      value = 0.7, step = 0.05),
          
         ),

        mainPanel(
          tabsetPanel(
          ### Input info ---- 
            tabPanel( "Info",
                    h2("Input options"),
                    h3("Data pre-processing"),
                    p(style = 'text-align: justify;',"Expected data is a file of normalised gene logCPMs with the gene names as rownames or as the first column of the dataframe."),
                    p(style = 'text-align: justify;',"Gene names and ENSEMBL gene ids are supported."),
                    p(style = 'text-align: justify;',"The model supports GRCh38 annotation only."),
                    h3("Data format"),
                    p(style = 'text-align: justify;',"Please indicate the following in the upload panel:"),
                    tags$ul(
                      tags$li(strong("Header:"), "wether the data contains a header. Headers with the sample identifiers are recommended."),
                      tags$li(strong("Genes as rownames:"), "if genes are the rownames of the provided dataframe (gene column does not have a column name) please tick this box. Otherwise, if the column containing gene indentifiers has a column name (i.e. Genes, GeneID, etc.) please untick this box."),
                      tags$li(strong("Separator:"), "indicate the separator in your input file.")
                       ),
                    h3("Missing genes"),
                    p(style = 'text-align: justify;',"Genes that cannot be mapped will be assigned a missing value. This can affect the performance of the model. We recomend checking gene identifiers."),
                    h3("Score threshold"),
                    p(style = 'text-align: justify;',"Threshold to be used to classify between high and low risk. By default 0.7 is selected (for further details see publication). Patients above threshold should be classified as high-risk and patients below the threshold as low-risk."),
                    hr(),
                    h2("Outputs"),
                    h3("Gene expression"),
                    p(style = 'text-align: justify;',"Gene expression tab contains logCPM gene expression in the provided data. In the visualization a dashed line is added for the mean logCPM expression in the TARGET cohort, where the model is trained."),
                    h3("Model scores"),
                    p(style = 'text-align: justify;',"Scores according to our risk model are displayed for the provided data. Below the Kaplan-Meyer and score distribution for the TARGET cohort are displayed according to the selected threshold."),
                    p(style = 'text-align: justify;',"The model is described in the publication ",  a(href = 'https://stackoverflow.com/', '[REF]', .noWS = "outside"), '.', .noWS = c("after-begin", "before-end")),
                    h3("Explore TARGET"),
                    p(style = 'text-align: justify;', "Visualization of score values and expression of genes in the model against available clinical variables for the TARGET cohort."),
                    p(style = 'text-align: justify;',"TARGET samples are available from  the TARGET data portal at National Cancer institute (NIH): ",  a(href = 'https://ocg.cancer.gov/programs/target/data-matrix', 'https://ocg.cancer.gov/programs/target/data-matrix', .noWS = "outside"), '.', .noWS = c("after-begin", "before-end")),
                    
            ),
            
          ### Output: Expression of signature genes ---- 
            tabPanel("Gene expression",  
                     br(), textOutput("myFile"), 
                     br(), textOutput("genesFound"), 
                     hr(), DT::dataTableOutput("genes"), 
                     br(), plotlyOutput("expression.plot"),br()),
            
          ### Output: Risk scores ----
            tabPanel("Model scores",
                     br(), h3("Scores in the provided cohort:"),
                     br(), textOutput("myFile2"),
                     br(), DT::dataTableOutput("scores.data"),
                     br(), plotlyOutput("scores.plot"),
                     hr(),
                     br(), h3("Kaplan-Meyer of TARGET cohort with the given score threshold:"),
                     br(), plotlyOutput("KM.target"), 
                     br(), plotlyOutput("scores.target")),
          
          ### Output: Plots with TARGET data ----- 
          tabPanel("Explore TARGET",
                   br(), selectInput("clinVar", "Clinical Variable", choices = ""),
                   br(), selectInput("value", "Values", choices = ""),
                   br(), plotlyOutput("explorePlots"))
          
        )
      )
    ) #sidebar closing
    
    ),#MAIN navbar closing
    
    tabPanel("Downloads",
    h2("TARGET lCPM"), 
    p(style = 'text-align: justify;',"Expression data in the TARGET cohort for the genes in the signature. Expression is provided as the log Counts Per Milion (logCPM)."),
    a(href="TARGET_lcpm_signature.txt", "TARGET_lcpm_signature.txt", download=NA, target="_blank"),
    h2("TARGET clinical information"), 
    p(style = 'text-align: justify;',"Clinical information displayed for the TARGET cohort."),
    a(href="TARGET_clinical_info.txt", "TARGET_clinical_info.txt", download=NA, target="_blank")
    ),
  )
)

# Define server  ---------

server <- function(input, output, session) {
  
  ### Read input  ---------
    counts = reactive({
      req(input$upload)
      inFile = input$upload
      if (is.null(inFile)){
        d = ""
      } else {
        d = read.table(inFile$datapath,header = input$header, sep = input$sep, stringsAsFactors = F)
        if (input$asrownames == F) {
          rownames(d) = d[,1]
          d = d[,-1]
        }
      }
      d
  })
    
    ### File not provided ---------
    output$myFile <- renderText({
      # Test if file is selected
      if (!is.null(input$upload$datapath)) {
        # Extract file name (additionally remove file extension using sub)
        return(NULL)
      } else {
        return("No input provided.")
      }
    })
    
    ### File not provided 2 ---------
    output$myFile2 <- renderText({
      # Test if file is selected
      if (!is.null(input$upload$datapath)) {
        # Extract file name (additionally remove file extension using sub)
        return(NULL)
      } else {
        return("No input provided.")
      }
    })
    
    
    
    ### Match annotation and select genes  ---------
    expression = reactive({
      countsTable = counts()
      if (grepl("ENSG",rownames(countsTable)[1])){
        rownames(countsTable) = str_split_fixed(rownames(countsTable),fixed("."), 2)[,1]
        countsTable = countsTable[rownames(countsTable) %in% ens2gene$ensembl_gene_id,]
        rownames(countsTable) = ens2gene[match(rownames(countsTable), ens2gene$ensembl_gene_id),]$external_gene_name
        a = countsTable
        #annot = paste(a, "out of", length(modelGenes), "signature genes found in expression data")
      } else {
        a = countsTable[rownames(countsTable) %in% modelGenes,]
        #annot = paste(a, "out of", length(modelGenes), "signature genes found in expression data")
      }
      notfound = setdiff(ens2gene$external_gene_name, rownames(a))
      if ( length(notfound) > 0 ) {
        toadd = setNames(data.frame(matrix(ncol = length(colnames(a)), nrow = length(notfound))), colnames(a))
        rownames(toadd) = notfound
        a = rbind(a, toadd)
      }
      a = a[order(rownames(a)),]
      a
   })
    
    ### Output genes found ---------
    output$genesFound = renderText({ 
      expr = expression()
      if (nrow(expr) > 0) {
        validGenes = sum(1*apply(expr, 1, function(x) !any(is.na(x))))
        c = paste(validGenes, "/", nrow(ens2gene), "genes found in expression data")
      } else {
        c = "No genes found, please check your data format and your gene identifiers"
      }
      c
    })
    

    
    ### Table gene expression  ---------
    output$genes = DT::renderDataTable({
      exprtable = expression()
      DT::datatable(exprtable, extensions = 'Buttons', filter = 'top', options = list(dom = 'Blfrtip', scrollX = TRUE, scrollY= "400px", scrollCollapse = T,  buttons = c('copy', 'csv', 'excel'), paging = F))
    })
    
    
    ### Plot gene expression  ---------
    output$expression.plot = renderPlotly({
      plotdata = expression()
      plotdata = melt(t(plotdata))
      colnames(plotdata) = c("Sample", "Gene", "Expression")
      plotdata.target = target.lcpm
      plotdata.target = rowMeans(target.lcpm)
      plotdata.target = melt(t(plotdata.target))
      colnames(plotdata.target) = c("Sample", "Gene", "Expression")
      
      plotdata$Sample = as.character(plotdata$Sample)
      plotdata.target$Sample = as.character(plotdata.target$Sample)
      plotdata = rbind(plotdata, plotdata.target)
      plotdata = plotdata[order(plotdata$Gene),]
      p <- plot_ly(plotdata, width = 11 * 96, height = 8 * 96) %>%
        add_boxplot(y = ~Expression, color = ~Gene, x= ~Gene, customdata= ~Sample, alpha = ~0.2, type = "box",colors = col_vector,
          transforms = list(
            list(
            type = 'filter',
            target = 'customdata',
            operation = '!=',
            value = '1'
          ))) %>%
       add_markers(y = ~Expression, x = ~Gene,  customdata= ~Sample, name = 'TARGET mean', line=list(color='grey', width=1, dash='dash', alpha = 0.5),
           transforms = list(
            list(
              type = 'filter',
              target = 'customdata',
              operation = '=',
              value = '1'
             )))    %>% 
        config(toImageButtonOptions = list(format = "svg"))
      p
    })
    
    ### Apply predictor ---------
    predictModel = reactive ({ 
      exprPredict = expression()
      exprPredict[is.na(exprPredict)] = 0
      predict.model = predict(leuko_rf_target, t(exprPredict), type="prob") 
      predict.model
      })
    
    ### Table of scores  ---------
    output$scores.data = DT::renderDataTable({
      scoresTable = as.data.frame(predictModel())
      colnames(scoresTable)[2] = "Score"
      scoresTable$Sample = rownames(scoresTable)
      DT::datatable(scoresTable[,c("Sample","Score")],rownames = F, extensions = 'Buttons', filter = 'top', options = list(dom = 'Blfrtip', scrollX = TRUE, scrollY= "400px", scrollCollapse = T,  buttons = c('copy', 'csv', 'excel'), paging = F))
    })
    
    ### Scores histogram wtih threshold line ---------
    output$scores.plot = renderPlotly({
      plotdata = as.data.frame(predictModel())
      colnames(plotdata)[2] = "Score"
      
      p <- plot_ly(plotdata, width = 6 * 96, height = 3 * 96) %>% 
        add_histogram(x = ~Score, name = "samples", marker = list(color = "lightgray",line = list(color = "darkgray", width = 2)))  %>%
        layout(shapes = list(vline(input$threshold))) %>% 
        config(toImageButtonOptions = list(format = "svg"))
      p 
    })
    
    ### KM ---------
    output$KM.target = renderPlotly({
        target.model = predict(leuko_rf_target, t(target.lcpm), type="prob") 
        modelTable =  merge(target.info, target.model, by.x = "Sample", by.y = "row.names")
        modelTable$class = "Low risk"
        modelTable[modelTable$Relapse >= input$threshold,]$class = "High risk"
        modelTable$First.Event.Logical = modelTable$First.Event == "Relapse"
        fit.trial = survfit(Surv(as.numeric(Event.Free.Survival.Time.in.Days),First.Event.Logical) ~ class,data=modelTable)
        p = ggsurvplot(fit.trial,data= modelTable, risk.table = TRUE, pval=TRUE, palette = c("#C93312","#899DA4"))
        ggplotly(p[[1]], width = 8 * 96)  %>% 
          config(toImageButtonOptions = list(format = "svg"))
        
    })
    
    
    exploreTable = reactive({
      lcpm = t(target.lcpm)
      target.model = predict(leuko_rf_target, t(target.lcpm), type="prob") 
      colnames(target.model)[2] = "Score"
      infotable =  merge(target.info, target.model, by.x = "Sample", by.y = "row.names")
      infotable = merge(infotable, lcpm, by.x = "Sample", by.y = "row.names")
      infotable
    })
    
    
    output$scores.target = renderPlotly({
      plotdata = as.data.frame(exploreTable())
      p <- plot_ly(plotdata, width = 6 * 96, height = 3 * 96) %>% 
        add_histogram(x = ~Score, name = "samples", marker = list(color = "lightgray",line = list(color = "darkgray", width = 2)))  %>%
        layout(shapes = list(vline(input$threshold))) %>% 
        config(toImageButtonOptions = list(format = "svg"))
      p
    })
    
    ### Update clinVar input -----
    updateSelectInput(session, "clinVar", 
                      choices = colnames(target.info %>% select(-Event.Free.Survival.Time.in.Days, -Age, -Sample, -Patient_ID)))
    
    ### Update clinVar input -----
    updateSelectInput(session, "value", 
                      choices = c("Score", modelGenes))
    
    output$explorePlots = renderPlotly({
      plotdata = as.data.frame(exploreTable())
      
      if (input$value == "Score"){
        p <- plot_ly(plotdata, width = 11 * 96, height = 8 * 96) %>%
          add_boxplot(y = ~get(input$value), color = ~get(input$clinVar), x= ~get(input$clinVar), alpha = ~0.2, type = "box",colors = col_vector) %>%
          layout(shapes = list(hline(input$threshold))) %>%
          layout(xaxis = list(title = input$clinVar, font = list(size = 14)), 
                 yaxis = list(title = input$value, font = list(size = 14))) %>% 
          config(toImageButtonOptions = list(format = "svg"))
        
        
      } else {
      p <- plot_ly(plotdata, width = 11 * 96, height = 8 * 96) %>%
        add_boxplot(y = ~get(input$value), color = ~get(input$clinVar), x= ~get(input$clinVar), alpha = ~0.2, type = "box",colors = col_vector) %>%
        layout(xaxis = list(title = input$clinVar, font = list(size = 14)), 
               yaxis = list(title = input$value, font = list(size = 14))) %>% 
        config(toImageButtonOptions = list(format = "svg"))
      }
      p
    })

}



# Run the application 
shinyApp(ui = ui, server = server)
