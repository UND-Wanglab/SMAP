tabItem(tabName = "input",
        tabsetPanel(tabPanel(p(id = "element","Variant Peptide Table"),
                             fluidRow(br(),
                                      box(title = "Upload Variant Peptide Table",
                                          radioButtons(inputId = "data", label = "Data Type", 
                                                       choices = c("User Data" = "user",
                                                                   "Example Data" = "example"),
                                                       selected = "example"),
                                          radioButtons(inputId = "sep", label = "Separator",
                                                       choices = c(Comma = ",",
                                                                   Semicolon = ";",
                                                                   Tab = "\t"),
                                                       selected = "\t"),
                                          radioButtons(inputId = "normalized", label = "Is Input Data Log2 Normalized?",
                                                       choices = c(Yes = "yes",
                                                                   No = "no"),
                                                       selected = "no"),
                                          fileInput(inputId = "inputFile", label = "Choose input file", multiple = F,accept = c("text/csv",
                                                                                                                                "text/comma-separated-values,text/plain",
                                                                                                                                ".csv")),
                                          
                                          HTML(paste(p(strong("Input file MUST follow this column specification:")),
                                                     "Peptide", HTML('&emsp;'), "Gene/Protein", HTML('&emsp;'), "PSM", HTML('&emsp;'), "SNP",
                                                     HTML('&emsp;'), "Sample 1 ... Sample N")), solidHeader = T, status = "primary")),
                             
                             fluidRow(
                               box(title = "Raw Data Preview (If Applicable):", DT::dataTableOutput(outputId = "rawhead"), width = 12, solidHeader = T, status = "primary")
                             ),
                             
                             fluidRow(
                               box(title = "Log2 Normalized Data Preview:", DT::dataTableOutput(outputId = "normhead"), width = 12, solidHeader = T, status = "primary")
                             )),
                    
                    tabPanel(p(id = "element","Genotype Data"),
                             fluidRow(br(),
                                      box(title = "Upload Genotype Data",
                                          radioButtons(inputId = "gdata", label = "Data Type",
                                                       choices = c("User Data" = "user",
                                                                   "Example Data" = "example"),
                                                       selected = "example"),
                                          fileInput(inputId = "genotypeFile", label = "Choose Genotype File", accept = ".vcf", multiple = F),
                                          HTML(paste(p(strong("Input file MUST follow VCF file format")),
                                                     p(strong("Missing genotypes MUST be coded as \"./.\"")))), hr(),
                                          HTML(paste(p(strong("**ID Column MUST match SNP column in variant peptide table**")))), solidHeader = T, status = "primary")),
                             
                             fluidRow(
                               box(title = "Genotype Data Preview:", DT::dataTableOutput(outputId = "vcfhead"), width = 12, solidHeader = T, status = "primary")
                             )
                    )))