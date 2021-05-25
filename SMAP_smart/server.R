

install.packages("???shinycssloaders")


options(shiny.maxRequestSize = 100*1024^2)


server <- function(input, output, session) {

  source("server/input_setup.R", local = T)
  
  source("server/summary.R", local = T)
  
  source("server/peptide_filtering.R", local = T)
  
  source("server/genotype_filtering.R", local = T)
  
  source("server/run_smap.R", local = T)
 
  
}
