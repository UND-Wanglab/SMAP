
## RUN SMAP


# fname <- eventReactive(input$runSMAP, {
#   return(paste0("tmp/output.",gsub(":",".",format(Sys.time(), format="%H:%M:%S"))))
# })

running <- reactiveVal("")



result <- eventReactive(input$runSMAP, {
  ## I will validate inputs first
  # Then read out inputs
  req(filtered())
  req(vcf.filt())
  
  running("yes")
  
  write.table(filtered(),"peptide_input.txt", row.names = F, quote = F, sep = "\t")
  
  v <- vcf.filt()
  colnames(v)[1] <- "#CHROM"
  write.table(v,"vcf_input.vcf", col.names = T, row.names = F, quote = F, sep = "\t")
  

  cmd <- paste0("perl new_SMAP.pl -fc 0 -nl 10000000000 --log2 1 -vf peptide_input.txt -g vcf_input.vcf -o temp")
  df <- data.frame(system(cmd, intern = T, wait = T)) %>%
    separate(colnames(.)[1], c("Sample ID","Inferred ID","Cscore","Delta Cscore"), sep = "\t")
  df <- df[-1,]
  rownames(df) <- NULL
  df$Cscore <- as.numeric(df$Cscore)
  df$`Delta Cscore` <- as.numeric(df$`Delta Cscore`)
  
  unlink(c("peptide_input.txt","vcf_input.vcf","temp", "sample_specific_genotype.snp",
           "sample_specific_genotype.vcf", "SMAP.log", "inferred_genotype.txt", "Score.txt"))
  
  running("no")
  
  return(df)
  
})



## Print results output
output$result.table <- DT::renderDataTable({
    req(result())
    df <- result()
    df$Cscore <- round(df$Cscore, digits = 2)
    df$`Delta Cscore` <- round(df$`Delta Cscore`, digits = 2)
    return(df)
    
}, options = list(scrollX = TRUE))




# ## Create Cscore plot
output$cscore <- renderPlot({
  req(result())
  df <- result()
  
  df <- df %>%
    mutate(Match = case_when(`Sample ID` == `Inferred ID` ~ "Yes",
                             `Sample ID` != `Inferred ID` ~ "No"))
  
  df$Match <- factor(df$Match, levels = c("Yes","No"))
  
  values = c("springgreen4","firebrick1")
  names(values) <- c("Yes","No")
 
  plot <- ggplot(df, aes(`Cscore`,`Delta Cscore`, label = `Inferred ID`, color = Match))+
    geom_point()+
    geom_text_repel(aes(color = Match), show.legend = F)+
    theme_bw()+
    labs(x = "Cscore", y = "Delta Cscore")+
    scale_color_manual(name = "Labels: Inferred ID\n\nMatches Sample ID:", values = values)+
    theme(axis.title = element_text(size = 15), 
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 12))
  
  return(plot)
})



## Option for downloading results
output$result.download <- downloadHandler(
  filename = function() {
    paste0(input$result.filename, ".csv", sep = "")
  },
  content = function(file) {
    df <- result()
    write.csv(df, file, row.names = F, quote = F)
  }
)
