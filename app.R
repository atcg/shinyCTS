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
load("freqThroughTimeArray.RData")
load("genotypes.RData")
load("hexSampsPerYear.RData")

### Functions
# Modified from the R cookbook to only fill internal gaps
fillInternalNAs <- function(x) {
  ## NA's in a vector or factor are replaced with last non-NA values
  ## If firstBack is TRUE, it will fill in leading NA's with the first
  ## non-NA value. If FALSE, it will not change leading NA's.

  # If it's a factor, store the level labels and convert to integer
  lvls <- NULL
  if (is.factor(x)) {
    lvls <- levels(x)
    x    <- as.integer(x)
  }
  goodIdx <- !is.na(x)
  
  if (sum(goodIdx) != 30 & sum(goodIdx) != 0) {
    firstGood <- min(which(goodIdx==TRUE))
    lastGood <- max(which(goodIdx==TRUE))
    
    for (i in firstGood:lastGood) {
      if (is.na(x[i])==TRUE) {
        x[i] = x[i-1]
      }
    }
  
    # If it was originally a factor, convert it back
    if (!is.null(lvls)) {
      x <- factor(x, levels=seq_along(lvls), labels=lvls)
    }
  }
  x
}


# Get the polygon coordinate and name data, which will be used upon brushing/selecting hexagons
polygonXcoords <- unlist(lapply(countiesWithMandersHexcolors@polygons, FUN=function(x){x@labpt}))[c(TRUE, FALSE)]
polygonYcoords <- unlist(lapply(countiesWithMandersHexcolors@polygons, FUN=function(x){x@labpt}))[c(FALSE, TRUE)]
polygonNames <- unlist(lapply(countiesWithMandersHexcolors@polygons, FUN=function(x){x@ID}))
percentageBTS <- countiesWithMandersHexcolors@data$PercBTS

polyz <- data.frame(polygonNames, polygonXcoords, polygonYcoords, percentageBTS, stringsAsFactors=FALSE)

ui <- fluidPage(
  fluidRow(
    column(width=5,
           plotOutput("map", height=900,click="map_click", brush=brushOpts(id="map_brush"))
    ),
    column(width=1,
           # Put the checkboxes here for what kind of alleles to light up in color
           checkboxGroupInput("checkGroup", label = h3("Select allele type to color"),
                              choices = list("All superinvasive SNPs" = 1, "Superinvasive SNP E23C6"=2, "Superinvasive SNP E12C11"=3, 
                                             "Superinvasive SNP E6E11"=4, "All putative BTS based\non reference pures"=5),selected=1)
    ),
    column(width=6,
           plotOutput("alleleFreq", height=900))
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
  
  output$value <- renderPrint({input$checkGroup })
  

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
    plot(manderCountiesUTMDissolved, axes=T, main="Select hexagons with mouse.")
    plot(countiesWithMandersHexcolors, col=rgb(countiesWithMandersHexcolors$col, maxColorValue = 255), add=T, lwd=0.001)
    #plot(countiesWithMandersHexcolors[match(countiesWithMandersHexcolors$id, IDZ)], lwd=2)
    plot(calCountiesUTM, add=T)
    
    #nearPointIDs <- selected_line$polygonNames
    brushPointIDs <- brushedPoints(polyz, input$map_brush, xvar="polygonXcoords", yvar="polygonYcoords")$polygonNames
    
    points(polyz$polygonXcoords[match(selected_line()$polygonNames,polyz$polygonNames)], polyz$polygonYcoords[match(selected_line()$polygonNames,polyz$polygonNames)], pch=20)
    points(polyz$polygonXcoords[match(brushPointIDs,polyz$polygonNames)], polyz$polygonYcoords[match(brushPointIDs,polyz$polygonNames)], pch=20)
    })
    
    output$alleleFreq <- renderPlot({
      # Plot the allele freq through time here
      pointedIDs <- as.character(selected_line()$polygonNames)
      brushPointIDs <- as.character(brushedPoints(polyz, input$map_brush, xvar="polygonXcoords", yvar="polygonYcoords")$polygonNames)
      cat(pointedIDs)
      cat(brushPointIDs)
      if (length(pointedIDs) > 0) {
        plot(1,ylim=c(0,1), xlim=c(1986,2015), main="Change in allele frequency over time in selected hexagons.\nInternal unsampled years set to previously sampled year.")
        for (allele in 1:length(freqThroughTimeArray[1,,1])) {
          # Only one hexagon so we don't need to calculate the allele frequency across hexagons
          lines(1986:2015, fillInternalNAs(freqThroughTimeArray[pointedIDs,allele,]), main=allele, type="l", col="black")
          # Now handle the colored lines
        }
      } else if (length(brushPointIDs) > 0) {
        plot(1,ylim=c(0,1), xlim=c(1986,2015), main="Change in allele frequency over time in selected hexagons.\nInternal unsampled years set to previously sampled year.")
        if (length(brushPointIDs) == 1) {
          for (allele in 1:length(freqThroughTimeArray[1,,1])) {
            lapply(brushPointIDs, FUN = function(x) {lines(1986:2015, fillInternalNAs(freqThroughTimeArray[x,allele,]), type="l", col="black")})
            # Now handle the colored lines
          } else { # More than one brushedPoint hexagon
            # Calculate frequencies across all selected hexagons:
            hexSampsPerYear[hexID,year] # This is the number of samples
            
            
          }
        }
      }
    })
}

shinyApp(ui = ui, server=server)











