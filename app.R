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
load("diagnosticGenos.RData")

diagnosticAlleles <- colnames(diagnosticGenos)

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
                              choices = list("All SI SNPs" = 1, "SI-E23C6"=2, "SI-E12C11"=3, 
                                             "SI-E6E11"=4, "All diagnostic BTS alleles based on reference pures"=5),selected=1)
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
  
  output$value <- renderPrint({input$checkGroup})
  

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
    plot(manderCountiesUTMDissolved, axes=T, main="Select hexagons with mouse.", cex.main=2)
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
      #cat(pointedIDs)
      #cat(brushPointIDs)
      
      # Save which IDs we want to print in green
      green1 <- 
      green2 <- 
      green3 <-
      green4 <-
      green5 <- diagnosticAlleles # This is a character vector of diagnostic allele names
      
      
      #(input$checkGroup)
      
      
      if (length(pointedIDs) > 0) {
        plot(1,type="n", ylim=c(0,1), xlim=c(1986,2015), main="Change in allele frequency over time in selected hexagons.\nInternal unsampled years set to previously sampled year.", xlab="Year", ylab="Allele Frequency")
        lapply(names(freqThroughTimeArray[1,,1]), FUN = function(x) {lines(1986:2015, fillInternalNAs(freqThroughTimeArray[pointedIDs,x,]), type="l", col="black")})
        # Now handle the colored lines
        lapply(green5, FUN = function(x) {lines(1986:2015, fillInternalNAs(freqThroughTimeArray[pointedIDs,x,]), type="l", lwd=2, col="green")})

        
      } else if (length(brushPointIDs) > 0) {
        plot(1, type="n", ylim=c(0,1), xlim=c(1986,2015), main="Change in allele frequency over time in selected hexagons.\nInternal unsampled years set to previously sampled year.", xlab="Year", ylab="Allele Frequency")
        
        if (length(brushPointIDs) == 1) {
          lapply(names(freqThroughTimeArray[1,,1]), FUN = function(x) {lines(1986:2015, fillInternalNAs(freqThroughTimeArray[brushPointIDs,x,]), type="l", col="black")})
          # Now handle the colored lines.
          # For instance: ((freqThroughTimeArray[1,"39457_contig100976|SNX13|3_contig100976|SNX13|3\t361",1])) is one allele
          lapply(green5, FUN = function(x) {lines(1986:2015, fillInternalNAs(freqThroughTimeArray[brushPointIDs,x,]), type="l", lwd=2, col="green")})

        } else { # More than one brushedPoint hexagon
            # Calculate and plot the frequencies across all selected hexagons:
            lapply(names(freqThroughTimeArray[1,,1]), FUN = function(x) {lines(1986:2015, fillInternalNAs(colSums(freqThroughTimeArray[brushPointIDs,x,]*hexSampsPerYear[brushPointIDs,x,], na.rm=TRUE)/colSums(hexSampsPerYear[brushPointIDs,x,], na.rm=TRUE)), type="l", col="black")})
          
            lapply(green5, FUN = function(x) {lines(1986:2015, fillInternalNAs(colSums(freqThroughTimeArray[brushPointIDs,x,]*hexSampsPerYear[brushPointIDs,x,], na.rm=TRUE)/colSums(hexSampsPerYear[brushPointIDs,x,], na.rm=TRUE)), type="l", lwd=2, col="green")})
        }
    } else { # This means that no polygons have been selected yet, so we'll just make an empty plot with axis labels
        plot(1,type="n", ylim=c(0,1), xlim=c(1986,2015), main="Change in allele frequency over time in selected hexagons.\nInternal unsampled years set to previously sampled year.", xlab="Year", ylab="Allele Frequency")
      }
    })
}

shinyApp(ui = ui, server=server)








