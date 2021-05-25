
## Set default filters
default.filters <- reactive({
  df <- norm.values()
  f <- c(0,0,0)
  names(f) <- c("min","max","sn")
  mean.min <- mean(df[which(df$min != 0),"min"])
  sd.min <- sd(df[which(df$min != 0),"min"])
  f["min"] <- as.numeric(mean.min + (2*sd.min))
  mean.max <- mean(df[,"max"])
  sd.max <- sd(df[,"max"])
  f["max"] <- as.numeric(mean.max - (2*sd.max))
  f["sn"] <- as.numeric(3)
  return(f)
})



## Default minimal intensity filter
output$minrecfilter <- renderText({
  paste0("Default value (if left blank above) is ", "<b>", round(default.filters()["min"], 2),"</b>"," (Mean + 2SD)")
})

## Default maximal intensity filter
output$maxrecfilter <- renderText({
  paste0("Default value (if left blank above) is ", "<b>", round(default.filters()["max"], 2), "</b>"," (Mean - 2SD)")
})

## Default S/N ratio filter
output$snrecfilter <- renderText({
  paste0("Default value (if left blank above) is ", "<b>", round(default.filters()["sn"], 2), "</b>")
})



## Create filtered data
filtered <- reactive({
  fdf <- norm.values()
  
  # Minimal filter (maximum bound)
  if (input$minfilter == "filter") {
    if (!isTruthy(input$minfiltervalue)) {
      fdf <- fdf[which(fdf$min < default.filters()["min"]),]
    } else {
      fdf <- fdf[which(fdf$min < input$minfiltervalue),]
    }
  } else {fdf <- fdf}
  
  # Maximal filter (minimum bound)
  if (input$maxfilter == "filter") {
    if (!isTruthy(input$maxfiltervalue)) {
      fdf <- fdf[which(fdf$max > default.filters()["max"]),]
    } else {
      fdf <- fdf[which(fdf$max > input$maxfiltervalue),]
    }
  } else {fdf <- fdf}
  
  # S/N ratio filter
  if (input$snfilter == "filter") {
    if (!isTruthy(input$snfiltervalue)) {
      fdf <- fdf[which(fdf$sn > default.filters()["sn"]),]
    } else {
      fdf <- fdf[which(fdf$sn > input$snfiltervalue),]
    }
  }  else {fdf <- fdf}
  
  return(fdf[,c(1:4,n())])
  
})


## Write final filter values and summary values
output$minimal.filter <- renderText({
  if (input$minfilter == "filter") {
    if (!isTruthy(input$minfiltervalue)) {
      paste0("Minimal Filter Value: ", round(default.filters()["min"], digits = 2))
    } else {
      paste0("Minimal Filter Value: ", round(input$minfiltervalue, digits = 2))
    }
  } else {paste0("Minimal Filter Value: NO FILTER")}
})
output$maximal.filter <- renderText({
  if (input$maxfilter == "filter") {
    if (!isTruthy(input$maxfiltervalue)) {
      paste0("Maximal Filter Value: ", round(default.filters()["max"], digits = 2))
    } else {
      paste0("Maximal Filter Value: ", round(input$maxfiltervalue, digits = 2))
    }
  } else {paste0("Maximal Filter Value: NO FILTER")}
})
output$sn.filter <- renderText({
  if (input$snfilter == "filter") {
    if (!isTruthy(input$snfiltervalue)) {
      paste0("S/N Ratio Filter Value: ", round(default.filters()["sn"], digits = 2))
    } else {
      paste0("S/N Ratio Filter Value: ", round(input$snfiltervalue, digits = 2))
    }
  } else {paste0("S/N Ratio Filter Value: NO FILTER")}
})
output$filt.num <- renderText({paste0("Total Rows After Filtering: ", nrow(filtered()))})
output$filt.unique <- renderText({paste0("Unique Peptides After Filtering: ", length(unique(filtered()[,1])))})


## Create table for head of filtered data
output$filter.table <- DT::renderDataTable({
  df <- filtered()
  df[,5:ncol(df)] <- round(df[,5:ncol(df)], digits = 2)
  return(head(df, n = 100))
}, options = list(scrollX = TRUE))


## Downloading filtered peptide data
output$peptide.download <- downloadHandler(
  filename = function() {
    if (input$pep.sep == "\t") {
      paste0(input$pep.filename,".txt", sep = "")
    } else if (input$pep.sep == "csv") {
      paste0(input$pep.filename,".csv",sep = "")
    }
  },
  content = function(file) {
    df <- filtered()
    write.table(df, file, row.names = F, quote = F, sep = input$pep.sep)
  }
)
