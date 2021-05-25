suppressPackageStartupMessages(require(shinydashboard))
suppressPackageStartupMessages(require(shiny))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(DT))
suppressPackageStartupMessages(require(vcfR))
suppressPackageStartupMessages(require(shinycssloaders))
suppressPackageStartupMessages(require(tidyr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggrepel))




ui <- dashboardPage(
  dashboardHeader(title = "SMAP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem(text = "Introduction", tabName = "intro", icon = icon("info-circle")),
      menuItem(text = "Input Data", tabName = "input", icon = icon("table")),
      menuItem(text = "Data Summary & Filtering", tabName = "summary", icon = icon("filter")),
      menuItem(text = "Run SMAP", tabName = "smap", icon = icon("laptop-code"))
    )
  ),
  dashboardBody(
    tags$head(
    tags$style("#element {
      font-size: 17px;
    }")),
    tabItems(
      source("ui/intro.R", local = T)$value,
      source("ui/input.R", local = T)$value,
      source("ui/summary.R", local = T)$value,
      source("ui/smap.R", local = T)$value)

    )
  
)


  
