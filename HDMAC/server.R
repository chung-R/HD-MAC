library(glmnet)
library(data.table)
library(shinyWidgets)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(selectiveInference)
# data.table: https://www.ptt.cc/bbs/R_Language/M.1437467101.A.E6D.html

source("c_index.R", encoding = "BIG5")
source("lr_gene_penalty_cv_function_version2.R")
source("cox_gene_penalty_cv_function_version2.R")
source("toy_example.R")
source("which_col.R")
source("borc.R")

# Define server logic to read selected file ----
options(shiny.maxRequestSize=100*1024^2) 

server <- function(input, output, session) {
  #### 1. Import Data ####
  #### 1-a. Load Data                                  ####
  #### Load data : User Own Data
  inData1 <- reactive({
    inFile <- input$file1
    if(is.null(inFile))return(NULL)
    data.frame(fread(inFile$datapath, 
                     header = input$header, 
                     stringsAsFactors = input$stringsAsFactors, 
                     sep = ","))
    })
  
  inData <- reactive({
    switch(input$select,
           "toy1" = Blad_SUB_8clinical_4937mut_8024mRNA,
           "toy2" = Ova_TCGA_OS_clinical_muta_cleaned_313_13_670,
           "uploaded" = inData1())
  })
  
  #### 1-b. row.col.sentence                           ####
  output$row.col.sentence = renderText({
    if(is.null(inData())){
      return(NULL)
    }else{
      paste(c("There are ", nrow(inData()), "patients and ", ncol(inData()), "variables.")) }
  })
  
  #### 1-c. DataTable                                  ####
  # datatable tutorial : https://shiny.rstudio.com/gallery/datatables-options.html
  output$contents <- renderDataTable(
    expr    = inData()[,c(input$ncol1:input$ncol2)],
    options = list(lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')), pageLength = 10)
  )
  
  #### 2. Binary Response (Logistic regression) ####
  #### 2-a. Model                                      ####
  output$m1 <- renderUI({
    withMathJax('$$P(Y_{i}=1|X_{i})=\\frac{e^{\\beta^{T}X_{i}}}{1+e^{\\beta^{T}X_{i}}} \\Leftrightarrow 
                   log\\left (\\frac{P(Y_{i}=1|X_{i})}{P(Y_{i}=0|X_{i})}\\right ) = \\beta^{T}X_{i}$$')
  })
  
  #### 2-b. Response                                   ####
  output$response_selector <- renderUI({
    if(is.null(inData())){
      return(NULL)
    }else if(input$select == "toy1"){
      pickerInput(
        inputId  = "selected_response", "Choose the response variable",
        choices  = colnames(inData())[1:15],
        selected = colnames(inData())[3],
        multiple = F, 
        options  = list('actions-box'=TRUE, title="Click here", 'deselect-all-text'="None", 'select-all-text'="All"),
        width = '70%'
      )
    }else{
      pickerInput(
        inputId  = "selected_response", "Choose the response variable",
        choices  = colnames(inData())[1:15],
        selected = NULL,
        multiple = F, 
        options  = list('actions-box'=TRUE, title="Click here", 'deselect-all-text'="None", 'select-all-text'="All"),
        width = '70%'
      )
    }
  })
  
  # update the slider called "cov_slider" based on the maxdistance
  observe({
    updateSliderInput(session, "cov_slider", max=length(inData())) 
  })
  
  output$cov_slider_100 <- renderUI({
    if(is.null(inData())){
      return(NULL)
    }else{
      pickerInput(
        inputId  = "cov_slider_100_sel", "",
        choices  = paste((((input$cov_slider)%/%100)*100+1):(((input$cov_slider)%/%100+1)*100),colnames(inData())[(((input$cov_slider)%/%100)*100+1):(((input$cov_slider)%/%100+1)*100)]),
        selected = paste(input$cov_slider,colnames(inData())[input$cov_slider]), #NULL,
        multiple = F, 
        options  = list('actions-box'=TRUE, title="Click here", 'deselect-all-text'="None", 'select-all-text'="All"),
        width = '80%'
      )
    }
  })
  
  #### 2-c. Clinical covaristes: X conti.              ####
  output$predictors_continuous_selector <- renderUI({
    if(length(input$selected_response) == 0){
      return(NULL)
    }else{
      pickerInput(
        inputId  = "selected_continuous_predictors", "Define continuous clinical covariates (optional)", 
        choices  = (colnames(inData())[!colnames(inData()) %in% input$selected_response])[1:15],
        selected = NULL,
        multiple = T, 
        options  = list('actions-box' = TRUE, title = "Click here", 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width    = '75%'
      )
    }
  })
  
  output$predictors_factor_selector <- renderUI({
    if(length(input$selected_response) == 0){
      return(NULL)
    }else{
      pickerInput(
        inputId  = "selected_factor_predictors", "Define categorical clinical covariates (optional)", 
        choices  = (colnames(inData())[!colnames(inData()) %in% c(input$selected_response,input$selected_continuous_predictors)])[1:15],
        selected = NULL,
        multiple = T, 
        options  = list('actions-box' = TRUE, title = "Click here", 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width    = '75%'
      )
    }
  })
  
  #### 2-d. Select clinical covaristes for mixed model ####
  output$covaristes_model <- renderUI({
    if(length(input$selected_response) == 0){
      return(NULL)
    }
    else{
      pickerInput(
        inputId  = "selected_covariates_model", "Choose clinical covariates to fit model (optional):", 
        choices  = (colnames(inData())[colnames(inData()) %in% c(input$selected_factor_predictors,input$selected_continuous_predictors)]),
        selected = NULL,
        multiple = T, 
        options  = list('actions-box' = TRUE, title = "Click here", 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width    = '75%'
      )
    }
  })
  
  #### 2-e. Define response colunm                     ####
  # reaction of 2-b. Response
  response.col <- reactive({
    inResp <- input$selected_response
    if (is.null(inResp))
      return(NULL)
    inResp
  })
 
  #### 2-f. Define X contiuous columns                 ####
  predictors.col <- reactive({
    inPred <- which( colnames(inData()) %in%  input$selected_continuous_predictors )
    if (is.null(inPred)) return(NULL)
    inPred
  })
  
  # Define X categorical columns
  predictors.cat.col <- reactive({
    inPred <- which( colnames(inData()) %in%  input$selected_factor_predictors )
    if (is.null(inPred)) return(NULL)
    inPred
  })
  
  # Define X categorical columns
  covariates.lr.model <- reactive({
    inPred <- which( colnames(inData()) %in% input$selected_covariates_model )
    if (is.null(inPred)) return(NULL)
    inPred
  })
  
  lr.covariates <- reactive({
    
    inPred <- which( colnames(inData()) %in%  input$selected_covariates_model )
    
    if (is.null(inPred)) return(NULL)
    
    inPred
    
  })
  
  #### 2-g. Run the result                             ####
  output$lr.gene_button <- renderUI({
    if(is.null(inData()))
      return(NULL)
    # icon tutorial: https://fontawesome.com/icons?d=gallery&q=play
    actionButton(inputId = 'lr.gene', 
                 label = "Run!",
                 icon("play"),
                 style = 'color: #fff; background-color: #337ab7; border-color: #2e6da4padding:4px; font-size:250%')
  })
  
  #### 2-h. Errpr messages                             #### 
  `%then%` <- shiny:::`%OR%`
  
  lr.gene.cv <- eventReactive(input$lr.gene, {
    if(input$fdr){
      validate(
        need(input$selected_response != "", 
             'Error:\n(Choose the response variable)\nYou must choose a response variable.\n') %then%
        need(borc(inData()[,which( colnames(inData()) %in%  response.col() )]) == 0, 
             'Error:\n(Choose the response variable)\nResponse variable must be binary (0 or 1).\n'),
        
        need(input$col_keyin != "", 
             'Error:\n(Columns of genetic covariates)\n“Columns of genetic variables” cannot be empty\n') %then%
        need(which_col(input$col_keyin)!= "", 
             "Error:\n(Columns of genetic covariates)\nPlease enter correct columns of genetic variables\n"),
        
        need(input$fdr_threshold != "", 
             'Error:\n(p-value threshold)\nPlease enter a number greater than 0 and less than 1 in “p-value threshold” (0.05 is recommended). If “use FDR for screening” is unchecked, then FDR screening will not be performed.\n') %then%
        need((input$fdr_threshold > 0)&(input$fdr_threshold < 1), 
             'Error:\n(p-value threshold)\nPlease enter a number greater than 0 and less than 1 in “p-value threshold” (0.05 is recommended). If “use FDR for screening” is unchecked, then FDR screening will not be performed.\n'),
        
        need(input$num.cv.fold != "", 
             'Error:\n(Number of CV folds)\nPlease enter a positive integer in “Number of CV folds” (5 is recommended). If 1 is entered, then cross validation method for estimating the prediction power will not be performed.') %then%
        need(is.integer(input$num.cv.fold)&input$num.cv.fold > 0 , 
             'Error:\n(Number of CV folds)\nPlease enter a positive integer in “Number of CV folds” (5 is recommended). If 1 is entered, then cross validation method for estimating the prediction power will not be performed.')
      )
    }
    else{
      validate(
        need(input$selected_response != "", 
             'Error:\n(Choose the response variable)\nYou must choose a response variable.\n') %then%
        need(borc(inData()[,which( colnames(inData()) %in%  response.col() )]) == 0, 
             'Error:\n(Choose the response variable)\nResponse variable must be binary (0 or 1).\n'),
        
        need(input$col_keyin != "", 
             'Error:\n(Columns of genetic covariates)\n“Columns of genetic variables” cannot be empty\n') %then%
        need(which_col(input$col_keyin)!= "", 
             "Error:\n(Columns of genetic covariates)\nPlease enter correct columns of genetic variables\n"),
             
        need(input$num.cv.fold != "", 
             'Error:\n(Number of CV folds)\nPlease enter a positive integer in “Number of CV folds” (5 is recommended). If 1 is entered, then cross validation method for estimating the prediction power will not be performed.') %then%
        need(is.integer(input$num.cv.fold)&input$num.cv.fold > 0 ,
             'Error:\n(Number of CV folds)\nPlease enter a positive integer in “Number of CV folds” (5 is recommended). If 1 is entered, then cross validation method for estimating the prediction power will not be performed.')
      )
    }
    
    if(is.na(input$col_keyin)){
      gene.gene.col = NULL
    }
    else{
      if(is.numeric(which_col(input$col_keyin))){
        gene.gene.col = which_col(input$col_keyin)
      } 
      else{
        gene.gene.col = NULL
        for (i in 1:length(which_col(input$col_keyin))) {
          gene.gene.col = c(gene.gene.col, which(names(inData())==which_col(input$col_keyin)[i]))
        }
        sort(gene.gene.col)
      }
    }
    lr.list = lr.gene.penalty.cv.result(user.data = inData(), 
                                        response.col  = response.col(), 
                                        gene.type     = list(gene.var = colnames(inData())[gene.gene.col]),
                                        padj.method   = if(input$fdr){"fdr"}else{"none"}, 
                                        fdr.threshold = input$fdr_threshold,
                                        penalty.alpha = input$penalty_select , 
                                        nfold         = input$num.cv.fold, 
                                        force.add.cov = if(length(lr.covariates()) == 0){NULL}else{
                                          lr.covariates()
                                        })
                                        #force.add.cov = covariates.lr.model())
    if (is.null(inData())) return(NULL)
    lr.list
  })
  
  #### 2-i. Sen, Spe, Accuracy, AUC, Gene List         ####
  output$lr.gene.cv.result.sen      = renderPrint({ lr.gene.cv()$sen })
  output$lr.gene.cv.result.spe      = renderPrint({ lr.gene.cv()$spe })
  output$lr.gene.cv.result.acy      = renderPrint({ lr.gene.cv()$accuracy })
  output$lr.gene.cv.result.auc      = renderPrint({ lr.gene.cv()$mean.auc })
  output$lr.gene.cv.result.lamd     = renderPrint({ lr.gene.cv()$ld })
  output$lr.gene.cv.result.coefandp = renderPrint({ lr.gene.cv()$coef.and.p })
  
  #### 3. Survival Response (Cox PH regression)
  #### 3-a. Model                                      ####
  output$m2 <- renderUI({
    withMathJax('$$h(t)=h_{0}(t)\\times e^{\\beta^{T}X_{i}}$$')
  })

  #### 3-b. cox_event: Status covariate                ####
  output$cox_event <- renderUI({
    if(is.null(inData())){
      return(NULL)
    }else if (input$select == "toy2") {
      pickerInput(
        inputId = "cox_selected_event", "Choose the event variable",
        choices = colnames(inData())[1:15],
        selected = colnames(inData())[10],
        multiple = F, options = list('actions-box' = TRUE, title = "Click here"
                                     , 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width = '70%'
      )
    }else{
      pickerInput(
        inputId = "cox_selected_event", "Choose the event variable",
        choices = colnames(inData())[1:15],
        selected = NULL,
        multiple = F, options = list('actions-box' = TRUE, title = "Click here"
                                     , 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width = '70%'
      )
    }
  })
  
  #### 3-c. cox_time: Time covariate                   ####
  output$cox_time <- renderUI({
    if(is.null(inData())){
      return(NULL)
    }else if (input$select == "toy2") {
      pickerInput(
        inputId = "cox_selected_time", "Choose the time variable",
        choices = colnames(inData())[1:15],
        selected = colnames(inData())[9],
        multiple = F, options = list('actions-box' = TRUE, title = "Click here"
                                     , 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width = '70%'
      )
    }else{
      pickerInput(
        inputId = "cox_selected_time", "Choose the time variable",
        choices = colnames(inData())[1:15],
        selected = NULL,
        multiple = F, options = list('actions-box' = TRUE, title = "Click here"
                                     , 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width = '70%'
      )
    }
  })
  
  #### 3-d. Clinical covaristes: X conti.              ####
  output$cox_predictors_selector <- renderUI({
    if(length(input$cox_selected_event) == 0|length(input$cox_selected_time) == 0){
      return(NULL)
    }else{
      pickerInput(
        inputId = "cox_selected_predictors", "Define continuous clinical covariates (optional)", 
        choices = (colnames(inData())[!colnames(inData()) %in% c(input$cox_selected_event, input$cox_selected_time) ])[1:10],
        selected = NULL,
        multiple = T, options = list('actions-box' = TRUE, title = "Click here"
                                     , 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width = '75%'
      )
    }
  })
  
  #### 3-e. Clinical covaristes: X categorical         ####
  output$cox_predictors_factor_selector <- renderUI({
    if(length(input$cox_selected_event) == 0|length(input$cox_selected_time) == 0){
      return(NULL)
    }else{
      pickerInput(
        inputId = "cox_selected_factor_predictors", "Define categorical clinical covariates (optional)", 
        choices = (colnames(inData())[!colnames(inData()) %in% c(input$cox_selected_event,input$cox_selected_time,input$cox_selected_predictors)])[1:10],
        selected = NULL,
        multiple = T, options = list('actions-box' = TRUE, title = "Click here"
                                     , 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width = '75%'
      )
    }
  })
  
  #### 3-f. Select clinical covaristes for mixed model ####
  output$cox_covaristes_model <- renderUI({
    if(length(input$cox_selected_event) == 0|length(input$cox_selected_time) == 0){
      return(NULL)
    }else{
      pickerInput(
        inputId = "cox_selected_covariates_model", "Choose the covariates to fit in model (optional):", 
        choices = (colnames(inData())[colnames(inData()) %in% c(input$cox_selected_factor_predictors,input$cox_selected_predictors)]),
        selected = NULL,
        multiple = T, options = list('actions-box' = TRUE, title = "Click here"
                                     , 'deselect-all-text' = "None", 'select-all-text' = "All"),
        width = '75%'
      )
    }
  })
  
  #### 3-g. Define response(time,event) col            ####
  cox.time.col <- reactive({
    inResp <- which( colnames(inData()) %in% input$cox_selected_time  )
    if (is.null(inResp))return(NULL)
    inResp
  })
  cox.event.col <- reactive({
    inResp <- which( colnames(inData()) %in% input$cox_selected_event  )
    if (is.null(inResp))return(NULL)
    inResp
  })
  
  #### 3-h. Define X conti col                         ####
  cox.predictors.col <- reactive({
    
    inPred <- which( colnames(inData()) %in%  input$cox_selected_predictors )
    
    if (is.null(inPred)) return(NULL)
    
    inPred
    
  })
  cox.predictors.cat.col <- reactive({
    
    inPred <- which( colnames(inData()) %in%  input$cox_selected_factor_predictors )
    
    if (is.null(inPred)) return(NULL)
    
    inPred
    
  })
  
  # update the slider called "cov_slider2" based on the maxdistance
  observe({
    updateSliderInput(session, "cov_slider2", max=length(inData())) 
  })
  output$cov_slider_100_2 <- renderUI({
    if(is.null(inData())){
      return(NULL)
    }else{
      pickerInput(
        inputId  = "cov_slider_100_sel_2", "",
        choices  = paste((((input$cov_slider2)%/%100)*100+1):(((input$cov_slider2)%/%100+1)*100),colnames(inData())[(((input$cov_slider2)%/%100)*100+1):(((input$cov_slider2)%/%100+1)*100)]),
        selected = paste(input$cov_slider2,colnames(inData())[input$cov_slider2]), #NULL,
        multiple = F, 
        options  = list('actions-box'=TRUE, title="Click here", 'deselect-all-text'="None", 'select-all-text'="All"),
        width = '80%'
      )
    }
  })
  
  #### 3-i. cox variable in mixed model                ####
  cox.covariates <- reactive({
    
    inPred <- which( colnames(inData()) %in%  input$cox_selected_covariates_model )
    
    if (is.null(inPred)) return(NULL)
    
    inPred
    
  })
  
  #### 3-j. Print columns number                       ####
  # Mutation col and Expression col
  output$gene.col <- renderText({ paste(c("Genetic covariates columns from",input$gene.col,"to",input$gene.col.end))  })
  
  #### 3-k. Run the result                             ####
  output$cox.gene_button <- renderUI({
    if(is.null(inData()))
      return(NULL)
    actionButton('cox.gene', 
                 "Run!",
                 icon("play"),
                 style = 'color: #fff; background-color: #337ab7; border-color: #2e6da4padding:4px; font-size:250%')
    
  })
  
  cox.gene.cv <- eventReactive(input$cox.gene, {
    
    if(input$fdr_cox){
    validate(
      need(input$cox_selected_time != "", 
           "Error:\n(Choose the time variable)\nYou must choose a time variable.\n") %then%
      need(borc(inData()[,which( colnames(inData()) %in%  input$cox_selected_time )]) == 1, #cox.time.col()
           'Error:\n(Choose the time variable)\nTime variable must be a continuous variable.\n'),
      
      need(input$cox_selected_event != "", 
           "Error:\n(Choose the event variable)\nYou must choose a event variable\n") %then%
      need(borc(inData()[,which( colnames(inData()) %in%  input$cox_selected_event )]) == 0, #cox.event.col()
           'Error:\n(Choose the event variable)\nEvent variable must be binary (0 or 1; 0 is censored).\n'),
      
      need(input$col_keyin2 != "", 
           'Error:\n(Columns of genetic covariates)\n“Columns of genetic variables” cannot be empty\n') %then%
      need(which_col(input$col_keyin2)!= "", 
           "Error:\n(Columns of genetic covariates)\nPlease enter correct columns of genetic variables\n"),
  
      need(input$fdr_threshold_cox != "", 
           'Error:\n(p-value threshold)\nPlease enter a number greater than 0 and less than 1 in “p-value threshold” (0.05 is recommended). If “use FDR for screening” is unchecked, then FDR screening will not be performed.\n') %then%
      need((input$fdr_threshold_cox > 0)&(input$fdr_threshold_cox < 1), 
           'Error:\n(p-value threshold)\nPlease enter a number greater than 0 and less than 1 in “p-value threshold” (0.05 is recommended). If “use FDR for screening” is unchecked, then FDR screening will not be performed.\n'),

      need(input$num.cv.fold_cox != "", 
           'Error:\n(Number of CV folds)\nPlease enter a positive integer in “Number of CV folds” (5 is recommended). If 1 is entered, then cross validation method for estimating the prediction power will not be performed.') %then%
      need(is.integer(input$num.cv.fold_cox)&input$num.cv.fold_cox > 0 , 
           'Error:\n(Number of CV folds)\nPlease enter a positive integer in “Number of CV folds” (5 is recommended). If 1 is entered, then cross validation method for estimating the prediction power will not be performed.')
    )}
    else{
      validate(
        need(input$cox_selected_time != "", 
             "Error:\n(Choose the time variable)\nYou must choose a time variable.\n") %then%
        need(borc(inData()[,which( colnames(inData()) %in%  input$cox_selected_time )]) == 1, #cox.time.col()
             'Error:\n(Choose the time variable)\nTime variable must be a continuous variable.\n'),
        
        need(input$cox_selected_event != "", 
             "Error:\n(Choose the event variable)\nYou must choose a event variable\n") %then%
        need(borc(inData()[,which( colnames(inData()) %in%  input$cox_selected_event )]) == 0, #cox.event.col()
             'Error:\n(Choose the event variable)\nEvent variable must be binary (0 or 1; 0 is censored).\n'),
        
        need(input$col_keyin2 != "", 
             'Error:\n(Columns of genetic covariates)\n“Columns of genetic variables” cannot be empty\n') %then%
        need(which_col(input$col_keyin2)!= "", 
             "Error:\n(Columns of genetic covariates)\nPlease enter correct columns of genetic variables\n"),
        
        need(input$num.cv.fold_cox != "", 
             'Error:\n(Number of CV folds)\nPlease enter a positive integer in “Number of CV folds” (5 is recommended). If 1 is entered, then cross validation method for estimating the prediction power will not be performed.') %then%
        need(is.integer(input$num.cv.fold_cox)&input$num.cv.fold_cox > 0 , 
             'Error:\n(Number of CV folds)\nPlease enter a positive integer in “Number of CV folds” (5 is recommended). If 1 is entered, then cross validation method for estimating the prediction power will not be performed.')
      )
    }
    
    if(is.na(input$col_keyin2)){
      gene.gene.col.cox = NULL
    }
    else{
      if(is.numeric(which_col(input$col_keyin2))){
        gene.gene.col.cox = which_col(input$col_keyin2)
      } 
      else{
        gene.gene.col.cox = NULL
        for (i in 1:length(which_col(input$col_keyin2))) {
          gene.gene.col.cox = c(gene.gene.col.cox, which(names(inData())==which_col(input$col_keyin2)[i]))
        }
        sort(gene.gene.col.cox)
      }
    }
    
    cox.list = cox.gene.penalty.cv.result(user.data = inData(),
                                          time.col  = cox.time.col(), 
                                          event.col = cox.event.col(), 
                                          gene.type =  list(gene.var = colnames(inData())[gene.gene.col.cox]),
                                          padj.method   = if(input$fdr_cox){"fdr"}else{"none"}, 
                                          fdr.threshold = input$fdr_threshold_cox,
                                          penalty.alpha = input$penalty_select_cox , 
                                          nfold = input$num.cv.fold_cox, 
                                          force.add.cov = if(length(cox.covariates()) == 0){NULL}else{
                                            cox.covariates()
                                          }
    )
    if (is.null(inData())) return(NULL)
    cox.list
  })
  
  #### 3-l. c-index, Gene List                         ####
  output$cox.gene.cv.result.cindex = renderPrint({
    cox.gene.cv()$c_index
  })
  output$cox.gene.cv.result.genelist = renderPrint({
    cox.gene.cv()$selected.gene
  })
  
  output$cox.gene.cv.result.coefandp = renderPrint({ cox.gene.cv()$coef.and.p})
  
}