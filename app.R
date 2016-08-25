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
percentageBTS <- countiesWithMandersHexcolors@data$PercBTS

polyz <- data.frame(polygonNames, polygonXcoords, polygonYcoords, percentageBTS, stringsAsFactors=FALSE)

ui <- fluidPage(
  fluidRow(
    column(width=6,
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
  
  # Store the map single click information
  click_saved <- reactiveValues(singleclick = NULL)
  observeEvent(eventExpr = input$map_click, handlerExpr = { click_saved$singleclick <- input$map_click })
  selected_line <-  reactive({
    nearPoints(polyz, click_saved$singleclick, ## changed from "input$plot_click" to saved click.
               xvar="polygonXcoords", yvar="polygonYcoords",
               maxpoints = 1,
               addDist = TRUE)
  })
  
  output$click_info <- renderPrint({
    res <- selected_line()
    ## datatable(res)
    res
  })
  
  
  output$brush_info <- renderPrint({
    brushedPoints(polyz, input$map_brush, xvar="polygonXcoords", yvar="polygonYcoords")
  })
  #    if (nearPoints(polyz, input$map_click, xvar="polygonXcoords", yvar="polygonYcoords")$polygonNames > 0) {
  #      selectedPoints <- nearPoints(polyz, input$map_click, xvar="polygonXcoords", yvar="polygonYcoords")$polygonNames
  #    }
    output$map <- renderPlot({
#    selectedPoints <- vector()
#    if (nearPoints(polyz, input$map_click, xvar="polygonXcoords", yvar="polygonYcoords")$polygonNames > 0) {
#      selectedPoints <- nearPoints(polyz, input$map_click, xvar="polygonXcoords", yvar="polygonYcoords")$polygonNames
#    } else if (brushedPoints(polyz, input$map_brush, xvar="polygonXcoords", yvar="polygonYcoords")$polygonNames > 0) {
#      selectedPoints <- brushedPoints(polyz, input$map_brush, xvar="polygonXcoords", yvar="polygonYcoords")$polygonNames
#    }
    
    plot(manderCountiesUTMDissolved, axes=T)
    plot(countiesWithMandersHexcolors, col=rgb(countiesWithMandersHexcolors$col, maxColorValue = 255), add=T, lwd=0.001)
    #plot(countiesWithMandersHexcolors[match(countiesWithMandersHexcolors$id, IDZ)], lwd=2)
    plot(calCountiesUTM, add=T)
    
    #nearPointIDs <- selected_line$polygonNames
    brushPointIDs <- brushedPoints(polyz, input$map_brush, xvar="polygonXcoords", yvar="polygonYcoords")$polygonNames
    

    points(polyz$polygonXcoords[match(selected_line()$polygonNames,polyz$polygonNames)], polyz$polygonYcoords[match(selected_line()$polygonNames,polyz$polygonNames)], pch=20)
    points(polyz$polygonXcoords[match(brushPointIDs,polyz$polygonNames)], polyz$polygonYcoords[match(brushPointIDs,polyz$polygonNames)], pch=20)
    
    
    })
  

}

shinyApp(ui = ui, server=server)











