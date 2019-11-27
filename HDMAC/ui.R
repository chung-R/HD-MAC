library(shiny)
library(glmnet)
library(data.table)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(shinyBS)
library(selectiveInference)
# high-dimensional analysis of molecular alterations in cancer

shinyUI(
  navbarPage("HD-MAC",
             #### table1: Upload Data         ####
             tabPanel("Upload data",
                      
                      # Sidebar layout with input and output definitions 
                      sidebarLayout(
                        # Sidebar panel for inputs
                        sidebarPanel(
                          
                        #### Input: Select toyexample or upload your own data ####
                        selectInput(inputId = "select",
                                    label = "Select toy example or upload your own data",
                                    choices = c(Toy_example_logistic = "toy1", Toy_example_survival = "toy2", Upload_data = "uploaded"),
                                    selected = "toy1"),
                          
                        #### Panel of toy example : ####
                        conditionalPanel(condition = "input.select == 'toy1'",
                                           
                                         # show colunms
                                         numericInput("ncol1", "Show columns of data from ", 1 , min = 1, max = 27380),
                                         numericInput("ncol2", "to", 10, min = 1, max = 27380)
                                         ),
                        conditionalPanel(condition = "input.select == 'toy2'",
                                         
                                         # show colunms
                                         numericInput("ncol1", "Show columns of data from ", 1 , min = 1, max = 27380),
                                         numericInput("ncol2", "to", 10, min = 1, max = 27380)
                        ),  
                        #### Panel of uploaded data : ####
                        conditionalPanel(condition = "input.select == 'uploaded'",
                                         # Input: Select a file
                                         fileInput("file1", "Upload your own csv file",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv")),
                                         
                                         # annotation:
                                         helpText("1. Data size limit is 50 MB. Please upload .csv file."),
                                         helpText("2. Your data should be consisted of response variable, genetic covariates (necessary) and clinical ones (optional)."),
                                         helpText("3. The clinical covariates (including response) should precede the genetic."),
                                         tags$hr(),
                                         
                                         # Input: Checkbox if file has header
                                         checkboxInput("header", label = "Header", value = TRUE),
                                         checkboxInput("stringsAsFactors", label = "stringsAsFactors", value = TRUE),  
                                         
                                         # Horizontal line
                                         tags$hr(),
                                           
                                         # show colunms
                                         numericInput("ncol1", "Show columns of data from ", 1 , min = 1, max = 27380),
                                         numericInput("ncol2", "to", 10, min = 1, max = 27380)
                                         ),
                        helpText('1. If you select the toy_example_logistic, you can run the data by clicking the "Binary Reponse: Logistic Regression" tab and hitting "Run" on the page. If you select the toy_example_survival, you can run the data by clicking the "Survival Response: CoxPH Regression" tab and hitting "Run" on the page. Otherwise, choose "Upload_data" to upload your own data.'),
                        helpText("2. For those who want to use TCGA dataset, you can download it through Broad's FireBrowse portal website (http://firebrowse.org/), use R package TCGAlinks (a Bioconductor package) or use our R code (filename: Data Download From TCGA.R ; link: https://github.com/chung-R/HD-MAC/blob/master/Data%20Download%20From%20TCGA.R).")
                        ),
                        
                        #### Main panel for displaying outputs ####
                        mainPanel(
                          # server: 67 rows
                          textOutput("row.col.sentence"),
                          
                          # Output: Data file
                          dataTableOutput('contents')
                          # tableOutput("contents")
                        )
                      )
             ), # tabPanel: "Upload Data" -- End
             
             
             #### table2: Logistic Regression ####
             tabPanel("Binary Response : Logistic Regression Model",
                      tags$head(
                        tags$style(HTML(".shiny-output-error-validation {color: red;}"))
                      ),
                      #### model ####
                      h3("Logistic Regression Model:"),
                      uiOutput('m1'),
                      hr(),
                      
                      #### Input: Select response ####
                      h3("Response"),
                      uiOutput("response_selector"),
                      hr(),
                      
                      #### Input: Covariates(genetic and clinical) ####
                      h3("Covariates"),
                      fluidRow(
                        column(4,
                               # genetic covariate column
                               textInput('col_keyin', 'Columns of genetic covariates', NULL, value = "4948-12971", width = '90%'), #placeholder = "100-120,150-155,197-200"),
                               sliderInput("cov_slider", " ", min=1, max=8000, value=4948, width = '550px'),
                               helpText("In case you need to know which column corresponds to what feature you can use this slider. 
                                                         Drag the slider to any column and the adjacent columns with column numbers and variable names will show on the right.")),
                        column(2,
                               helpText("e.g.1:100-120"),
                               helpText("e.g.2:100-120,150,197-200"),
                               hr(),
                               uiOutput("cov_slider_100")
                               #uiOutput("tt") #### 1122
                               ),
                        column(3,
                               uiOutput("predictors_continuous_selector"),
                               uiOutput("predictors_factor_selector"),
                               uiOutput("covaristes_model"),
                               helpText("PS:Binary covariates should be represented by 0 and 1."))
                      ),
                      hr(),
                      
                      #### Input: Regression penalty ####
                      h3("Choose regression penalty for the model"),
                      selectInput('penalty_select', 'Regression penalty', c("Ridge" = 0,"Lasso" = 1,"Adaptive Lasso" = 2), selected = 2),
                      hr(),
                      
                      #### Input: Screening ####
                      h3("Screening (optional)"),
                      checkboxInput("fdr", "Use FDR for screening", T),
                      numericInput('fdr_threshold', 'p-value threshold', 0.05, min = 0, max = 1, step = 0.01),
                      hr(),
                     
                      #### Input: number of cross-validation folds ####
                      h3("Cross-Validation for prediction power (optional)"),
                      numericInput('num.cv.fold', 'Number of CV folds', 1, min = 1, max = 10, step = 1),
                      hr(),
                      
                      #### Run Button: Method ####
                      fluidRow(
                        column(4, uiOutput('lr.gene_button'))
                      ),
                      hr(),
                      
                      #### Output: Final Result ####
                      h3("Final Result"),
                      hr(),
                      
                      tabsetPanel(
                        tabPanel("Gene List, estimated coefficients and p-values", 
                                 h4("Gene list, estimated coefficients and p-values"),
                                 verbatimTextOutput("lr.gene.cv.result.coefandp")
                                 ),
                        
                        tabPanel("Prediction Power", 
                                 h4("Sensitivity"),
                                 # helpText("Sensitivity (also called the true positive rate or probability of detection) measures the proportion of actual positives that are correctly identified."),
                                 verbatimTextOutput("lr.gene.cv.result.sen"),
                                 
                                 h4("Specificity"),
                                 # helpText("Specificity (also called the true negative rate) measures the proportion of actual negatives that are correctly identified."),
                                 verbatimTextOutput("lr.gene.cv.result.spe"),
                                 
                                 h4("Accuracy"),
                                 # helpText("Overall accuracy."),
                                 verbatimTextOutput("lr.gene.cv.result.acy"),
                                 
                                 h4("AUC (%)"),
                                 # helpText("Area Under the ROC Curve. Evaluate the performance of regression model."),
                                 verbatimTextOutput("lr.gene.cv.result.auc")
                                 )
                      )

             ),
             
             
             #### table3: coxph regression    ####
             tabPanel("Survival Response : CoxPH Model",
                      tags$head(
                        tags$style(HTML(".shiny-output-error-validation {color: red;}"))
                      ),
                      h3("Cox PH Model:"),
                      uiOutput('m2'),
                      hr(),
                      
                      # Input: Select number of rows to display
                      h3("Response"),
                      uiOutput("cox_time"),
                      uiOutput("cox_event"),
                      hr(),
                      
                      # Input: Covariates(genetic and clinical)
                      h3("Covariates"),
                      fluidRow(
                        column(4,
                               # genetic covariate column
                               textInput('col_keyin2', 'Columns of genetic covariates', NULL, value = "14-683", width = '90%'), #placeholder = "100-120,150-155,197-200"),
                               sliderInput("cov_slider2", "", min=1, max=8000, value=14, width = '550px'),
                               helpText("In case you need to know which column corresponds to what feature you can use this slider. 
                                                         Drag the slider to any column and the adjacent columns with column numbers and variable names will show on the right.")),
                        column(2,
                               helpText("e.g.1:100-120"),
                               helpText("e.g.2:100-120,150,197-200"),
                               hr(),
                               uiOutput("cov_slider_100_2")),
                        column(3,
                               uiOutput("cox_predictors_selector"),
                               uiOutput("cox_predictors_factor_selector"),
                               uiOutput("cox_covaristes_model"),
                               helpText("PS:Binary covariates should be represented by 0 and 1."))
                      ),
                      hr(),
                      
                      # Input: Regression penalty
                      h3("Choose regression penalty for the model"),
                      selectInput('penalty_select_cox', 
                                  'Regression penalty', 
                                  c("Ridge" = 0,"Lasso" = 1,"Adaptive Lasso" = 2),
                                  selected = 1),
                      hr(),
                      
                      # Input: Screening
                      h3("Screening (optional)"),
                      checkboxInput("fdr_cox", 
                                    "Use FDR for screening", 
                                    FALSE),
                      numericInput('fdr_threshold_cox', 
                                   'p-value threshold', 0.05, min = 0, max = 1,
                                   step = 0.01),
                      hr(),
                      
                      # Input: number of cross-validation folds
                      h3("Cross-Validation for prediction power (optional)"),
                      numericInput('num.cv.fold_cox', 
                                   'Number of cv folds', 1, min = 1, max = 10,
                                   step = 1),
                      hr(),
                      
                      # Run Button: Method
                      fluidRow(
                        column(4, uiOutput('cox.gene_button')),
                        column(4, textOutput("mut.col_cox")),
                        column(4,  textOutput("exp.col_cox"))),
                      hr(),
                      
                      # Output: Final Result
                      h3("Final Result"),
                      hr(),
                      
                      tabsetPanel(
                        tabPanel("Gene List, estimated coefficients and p-values", 
                                 # h4("Gene List"),
                                 # verbatimTextOutput("cox.gene.cv.result.genelist"),
                                 
                                 h4("Estimated coefficients and p-values"),
                                 verbatimTextOutput("cox.gene.cv.result.coefandp")
                        ),
                        
                        tabPanel("Prediction Power", 
                                 h4("C-index"),
                                 helpText("Concordance index. C-index it measures how well the model discriminates between different responses, i.e., is your predicted response low for low observed responses and high for high observed responses."),
                                 verbatimTextOutput("cox.gene.cv.result.cindex")
                        )
                      )
                        
             )             
             
  )
)
