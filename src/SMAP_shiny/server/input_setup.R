
## Get variant peptide table read in to be used
peptide <- reactive({
  if (input$data == "example") {
    tbl <- read.table(file = "data/variant_peptide_table.txt", header = TRUE, sep = "\t")
    return(tbl)
  } else if (input$data == "user") {
    validate(
      need(input$inputFile != "", "Please select a file")
    )
    tbl <- read.table(file = input$inputFile$datapath, header = TRUE, sep = input$sep)
    return(tbl)
  }
})




## Normalize variant peptide data with log2 transformation
normdf <- reactive({
  if (input$normalized == "no") {
    normtbl <- peptide()
    normtbl[,5:ncol(normtbl)] <- log2(normtbl[,5:ncol(normtbl)] + 1)
    return(normtbl)
  } else if (input$normalized == "yes") {
    normtbl <- peptide() 
    return(normtbl)
  }
})

## Create value for column numbers of samples
n <- reactive({
  col <- ncol(normdf())
  samples <- c(5:col)
  return(samples)
})


## Create an updated variant peptide dataframe to use with summary values
norm.values <- reactive({
  # Mean and sd
  df <- normdf()
  df[,"mean"] <- apply(df[,n()], 1, mean)
  df[,"sd"] <- apply(df[,n()], 1, sd)
  # Minimum value
  for (i in 1:nrow(df)) {
    df[i,"min"] <- min(as.numeric(df[i,n()]))
  }
  # Maximum value
  for (i in 1:nrow(df)) {
    df[i,"max"] <- max(as.numeric(df[i,n()]))
  }
  # S/N ratio
  for (i in 1:nrow(df)) {
    if (df[i,"min"] == 0) {
      df[i,"sn"] <- 10
    } else {
      df[i,"sn"] <- df[i,"max"] - df[i,"min"]
    }
    if (df[i,"sn"] > 10) {
      df[i,"sn"] <- 10
    }
  }
  return(df)
})


## Get genotype data read in to be used
vcf.df <- reactive({
  if (input$gdata == "example") {
    genotype <- read.delim("data/genotype_table.vcf", skip = 1)
    colnames(genotype)[1] <- "CHROM"
    return(genotype)
  } else if (input$gdata == "user") {
    validate(
      need(input$genotypeFile != "", "Please select a file")
    )    
    meta <- readLines(input$genotypeFile$datapath, n = 1000)
    skip <- length(which(substr(meta, 1, 2) == "##"))
    
    genotype <- read.delim(input$genotypeFile$datapath, skip = skip, header = T)
    
    colnames(genotype)[1] <- "CHROM"
    return(genotype)
  }
})



## Create dataframe of number of missing genotypes 
vcf.missing <- reactive({
  req(vcf.df())
  df <- vcf.df()
  df$ref.hom <- rowSums(df == "1/1")
  df$het <- rowSums(df == "0/1")
  df$alt.hom <- rowSums(df == "0/0")
  df$tot.samples <- rowSums(df == "1/1" | df == "0/1" | df == "0/0"| df == "./.")
  df$tot.alleles <- df$tot.samples*2
  
  df$ref.freq <- ((df$ref.hom * 2) + df$het) / df$tot.alleles
  df$alt.freq <- ((df$alt.hom * 2) + df$het) / df$tot.alleles
  df$n <- rowSums(df == "./.")
  
  for (i in 1:nrow(df)) {
    if (df[i,"ref.freq"] < 0.5) {
      df[i,"maf"] <- df[i,"ref.freq"]
    } else {
      df[i,"maf"] <- df[i,"alt.freq"]
    }
  }
  
  return(df)
})






