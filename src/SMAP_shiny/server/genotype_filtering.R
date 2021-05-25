
## Make missing genotype distribution
output$na.hist <- renderPlot({
  df <- vcf.missing()
  plot <- hist(df$n, breaks = input$na.bins, main = "", xlab = "Number of Missing Genotypes")
})




## Write table for filtered vcf file
vcf.filt <- reactive({
  df <- vcf.df()
  # if (input$filt.overlap == "filter") {
  #   df <- df[which(df$ID %in% filtered()$SNP),]
  # } else {df <- df}
  if (input$filtna == "filter") {
    missing <- vcf.missing()[which(vcf.missing()$n <= input$na.allowed),"ID"]
    df <- df[which(df$ID %in% missing),]
  } else {df <- df}
  return(df)
})


# Create allele-frequency dataframe from filtered vcf
af <- reactive({
  df <- vcf.missing()[which(vcf.missing()$ID %in% vcf.filt()$ID),]
  return(df)
})


## Create table for head of filtered data
output$vcf.filt.table <- DT::renderDataTable({
  df <- vcf.filt()
  return(head(df, n = 100))
}, options = list(scrollX = TRUE))


## Filter summary values for vcf
output$filt.snps <- renderText({paste0("Total SNPs After Filtering: ", nrow(vcf.filt()))})
output$filt.missing <- renderText({paste0("SNPs with At Least One Missing Genotype After Filtering: ", 
                                          nrow(vcf.missing()[which(vcf.missing()$n != 0 & vcf.missing()[,"ID"] %in% vcf.filt()$ID),]))})


## Downloading filtered genotype data
output$vcf.download <- downloadHandler(
  filename = function() {
    if (input$vcf.sep == "\t") {
      paste0(input$vcf.filename, ".txt", sep = "")
    } else if (input$vcf.sep == "csv") {
      paste0(input$vcf.filename, ".csv", sep = "")
    }
  },
  content = function(file) {
    df <- vcf.filt()
    write.table(df, file, row.names = F, quote = F, sep = input$vcf.sep)
  }
)

## Allele frequency distribution plots
output$ref.hist <- renderPlot({
  req(af())
  df <- af()
  plot <- hist(df$ref.freq, breaks = input$ref.bins, main = "", xlab = "Reference Allele Frequency")
  return(plot)
})
output$alt.hist <- renderPlot({
  req(af())
  df <- af()
  plot <- hist(df$alt.freq, breaks = input$alt.bins, main = "", xlab = "Alternate Allele Frequency")
  return(plot)
})

