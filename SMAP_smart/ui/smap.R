tabItem(tabName = "smap",
        fluidRow(
          box(title = "Run SMAP", hr(),
              "Before hitting run: Ensure all filters are set to the desired values in previous tab", br(), br(),
              actionButton(inputId = "runSMAP", label = "Run SMAP"), solidHeader = T, status = "primary", height = 200),
          box(title = "Download Results",
              textInput("result.filename", "Enter name for download file:", value = "SMAP_output"),
              downloadButton(outputId = "result.download", label = "Download"), solidHeader = T, status = "primary", height = 200)
        ),
        fluidRow(
          box(title = "SMAP Output", DT::dataTableOutput(outputId = "result.table") %>% withSpinner(type = 7), width = 6, solidHeader = T, status = "primary"),
          box(title = "Cscore and Delta Cscore Distribution", br(),br(),
              plotOutput(outputId = "cscore") %>% withSpinner(type = 7), width = 6, height = 544,solidHeader = T, status = "primary")
        )
)