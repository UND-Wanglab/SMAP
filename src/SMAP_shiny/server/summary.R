
## Table for raw variant peptide table
output$rawhead <- DT::renderDataTable({
  if (input$normalized == "no") {
    df <- peptide()
    df[,5:ncol(df)] <- round(df[,5:ncol(df)], digits = 0)
    return(head(df, n = 100))
  }
}, options = list(scrollX = TRUE))


## Table for normalized variant peptide table
output$normhead <- DT::renderDataTable({
  df <- normdf()
  df[,5:ncol(df)] <- round(df[,5:ncol(df)], digits = 2)
  return(head(df, n = 100))
}, options = list(scrollX = TRUE))


## Table for genotype input
output$vcfhead <- DT::renderDataTable({
  vcf <- vcf.df()
  return(head(vcf, n = 100))
}, options = list(scrollX = TRUE))





## Summary values for variant peptides
output$rows <- renderText(paste0("Total Rows: ",nrow(normdf())))
output$pep_num <- renderText(paste0("Number of Unique Peptides: ",length(unique(normdf()[,1]))))
output$gene_num <- renderText(paste0("Number of Unique Genes/Proteins: ",length(unique(normdf()[,2]))))
output$psm_num <- renderText(paste0("Number of peptide spectrum matches (PSMs): ",length(unique(normdf()[,3]))))
output$snp_num <- renderText(paste0("Number of SNPs: ",length(unique(normdf()[,4]))))

## Summary values for genotype data
output$geno.snps <- renderText(paste0("Total number of SNPs present in genotype input: ", nrow(vcf.df())))
output$geno.pepsnps <- renderText(paste0("Number of SNPs in both (unfiltered) variant peptide and genotype input: ", length(which(vcf.df()$ID %in% peptide()$SNP))))
output$geno.missing <- renderText(paste0("SNPs with at least one missing genotype: ", nrow(vcf.missing()[which(vcf.missing()$n != 0),])))
output$geno.samples <- renderText(paste0("Total number of samples present in genotype input: ", ncol(vcf.df()[,10:ncol(vcf.df())])))
output$geno.match <- renderText(paste0("Number of samples with an ID present in peptide input: ", length(which(colnames(vcf.df()) %in% colnames(peptide())))))





## Intensity plot for all peptides
output$intensity <- renderPlot({
  df <- norm.values()
  z <- as.matrix((df[,n()] - df[,"mean"]) / df[,"sd"])
  plot <- hist(z, breaks = input$ibins, main = "", xlab = "Z-score Transformed Intensity Values")
  return(plot)
})

## Minimal intensity plot
output$minimal <- renderPlot({
  df <- norm.values()
  plot <- hist(df[which(df$min != 0),"min"], breaks = input$minbins, main = "", xlab = "Log2 Normalized Intensity Values")
  return(plot)
})

## Maximal intensity plot
output$maximal <- renderPlot({
  df <- norm.values()
  plot <- hist(df$max, breaks = input$maxbins, main = "", xlab = "Log2 Normalized Intensity Values")
  return(plot)
})

## PCA plot
output$pca <- renderPlot({
  df <- norm.values()
  pca <- data.frame(prcomp(t(df[,n()]))$x)
  pca$sample <- rownames(pca)
  ggplot(pca, aes(PC1,PC2, label = sample))+
    geom_point()+
    geom_text_repel(show.legend = F)+
    theme_bw()+
    theme(legend.position = "none")
})


## S/N ratio plot
output$sn <- renderPlot({
  df <- norm.values()
  plot <- hist(df$sn, breaks = input$sbins, main = "", xlab = "Log2(Maximal/Minimal) Ratio")
  return(plot)
})



## Minor allele frequency distribution plot
output$maf.hist <- renderPlot({
  df <- vcf.missing()
  plot <- hist(df$maf, breaks = input$maf.bins, main = "", xlab = "Minor Allele Frequency")
})
