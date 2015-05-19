#ui.R

shinyUI(fluidPage(
  titlePanel("Biomass yields from native grasslands"),
  sidebarLayout(
    sidebarPanel(
      selectInput("variable", "Treatment:",
                  c("Switchgrass" = "1",
                    "Big Bluestem" = "2",
                    "Indiangrass" = "3",
                    "Canada Wild Rye" = "4",
                    "4 species grass mix" = "5",
                    "4 species legume mix" = "6",
                    "4 species forb mix" = "7",
                    "8 species grass-legume mix" = "8",
                    "8 species grass-forb mix" = "9",
                    "8 species forb-legume mix" = "10",
                    "12 species grass-forb-legume mix" = "11",
                    "24 species high-diversity mix" = "12"))
    ),
      mainPanel(
        plotOutput("yldPlot")
    )
  )
))