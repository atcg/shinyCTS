setwd("~/Box Sync/UCLA/Research/Papers/CTS_HybEx/pondsThroughTime/")



##### LIBRARIES #####
library(tigris)
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(dplyr)
library(tidyr)
library(rgbif)
library(viridis)
library(gridExtra)
library(rasterVis)
library(rnaturalearth)

library(shiny)


# Load in the data
load("countiesWithMandersHexcolors.RData")
load("manderCountiesUTMDissolved.RData")
load("calCountiesUTM.RData")

# Get the polygon coordinate and name data, which will be used upon brushing/selecting hexagons
polygonXcoords <- unlist(lapply(countiesWithMandersHexcolors@polygons, FUN=function(x){x@labpt}))[c(TRUE, FALSE)]
polygonYcoords <- unlist(lapply(countiesWithMandersHexcolors@polygons, FUN=function(x){x@labpt}))[c(FALSE, TRUE)]
polygonNames <- unlist(lapply(countiesWithMandersHexcolors@polygons, FUN=function(x){x@ID}))

polyz <- data.frame(polygonNames, polygonXcoords, polygonYcoords, stringsAsFactors=FALSE)

ui <- fluidPage(
  fluidRow(
    column(width=8,
           plotOutput("map", height=900,click="map_click", brush=brushOpts(id="map_brush"))
    )
  ),
  fluidRow(
    column(width = 6,
           h4("Points near click"),
           verbatimTextOutput("click_info")
    ),
    column(width = 6,
           h4("Brushed points"),
           verbatimTextOutput("brush_info")
   )
  )
)





server <- function(input, output) {
  output$map <- renderPlot({
    plot(manderCountiesUTMDissolved, axes=T)
    plot(countiesWithMandersHexcolors, col=rgb(countiesWithMandersHexcolors$col, maxColorValue = 255), add=T, lwd=0.001)
    plot(calCountiesUTM, add=T)
    points(polyz$polygonXcoords, polyz$polygonYcoords, pch=20)
    })
  
  output$click_info <- renderPrint({
    # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
    # were a base graphics plot, we'd need those.
    nearPoints(polyz, input$map_click, xvar="polygonXcoords", yvar="polygonYcoords", addDist = TRUE)
  })
  
  output$brush_info <- renderPrint({
    brushedPoints(polyz, input$map_brush, xvar="polygonXcoords", yvar="polygonYcoords")
  })
}

shinyApp(ui = ui, server=server)
