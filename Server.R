library(rgdal)
library(viridis)

fillmap2<-function(map, figtitle, y , leg.loc="beside", y.scl=NULL,
                   main.cex=1.5,main.line=0,map.lty=1,leg.rnd=0,
                   leg.cex=1){
  
  # 0: dark 1: light light Current shading ranges from darkest to light gray white (to distinguish with lakes).
  y.uq=sort(unique(c(y,y.scl)))
  cols<-viridis(length(y.uq),direction=-1)
  shading=y
  for (i in 1:length(y)){
    shading[i]<-cols[which(y.uq==y[i])]
  }
  
  par(mar=c(0,0,2,0))
  if (leg.loc=="beside"){
    layout(matrix(1:2,ncol=2),width=c(.8,.2))
  } else 
    if (leg.loc=="below"){
      layout(matrix(1:2,nrow=2),height=c(.6,.4))
    } else (print("leg.loc options are below or beside"))
  
  plot(map,col=shading,axes=F, lty=map.lty)
  title(main=figtitle,cex.main=main.cex,line=main.line) 
  
  par(mar=c(5, 4, 4, 2) + 0.1)
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
  cols.5=cols[seq(1,length(y.uq),length.out=5)]
  lab.5=cols.5
  for (i in 1:5){lab.5[i]=y.uq[which(cols==cols.5[i])[1]]}
  lab.5=round(as.numeric(lab.5),leg.rnd)
  par(mar=c(0,0,0,0))
  if (leg.loc=="beside"){
    legend_image <- as.raster(matrix(cols, ncol=1))
    text(x=1.6, 
         y = seq(0,length(y.uq),length.out=5)/length(y.uq),
         labels = rev(lab.5), cex=leg.cex)
    rasterImage(legend_image, 0, 0, 1,1)
  } else{
    legend_image <- as.raster(matrix(cols, nrow=1))
    text(y=-0.25, 
         x = seq(0,length(y.uq),length.out=5)/(length(y.uq)*.5),
         labels = lab.5, cex=leg.cex)
    rasterImage(legend_image, 0, 0, 2,1)
  }
}

data <- read.csv('Data\\Data4Shiny.csv')
head(data)

NHtracts=readOGR("Data\\NHtracts.shp")


#simply plot the arrest data
shinyServer(function(input,output){
  output$text <- renderText({
    if (input$adj=="None"){
      "No adjustments specified. These are counts of arrests for the selected data."
    } else
      if (input$adj=="Standardized Incidence Ratio"){
        "The standardized incidence ratio (SIR) adjustment was applied here. SIR is a method of adjusting for tract population by calculating a ratio of the observed count of arrests to the expected counts of arrests. The expected count of arrests is calculated as a rate of arrests over all tracts and years times the tract population for a given year. The population used reflects the data being displayed (e.g. Black population only when 'Black Only Arrests' is selected). Values greater than 1 suggest more observed arrests than expected."
      } else
        if (input$adj=="Poisson Regression"){
          "The Poisson regression option for adjustment was applied here. In this method of adjustment, a Poisson regression model with spatio-temporal covariate adjustment was applied (see table output). The mapped values display the residual spatial variation in arrests after adjustment where higher (darker) values indicate areas of increased risk. All estimates are transformed so that they can be interpreted as a multiplicative change in the relative risk of arrests. Tract population is indirectly adjusted."
        } else 
          if (input$adj=="As a Percent of the Population"){
            "The 'As a Percent of the Population' adjustment displays the selected arrests counts divided by the appropriate population (e.g. Black population only when 'Black Only Arrests' is selected) times 100%."
          } else{#% tot arrests
            "The 'As a Percent of Total Arrests' adjustment displays the selected arrests counts divided by the total arrest counts times 100%."
          }
  })
  
  output$map <- renderPlot({
    if (input$adj=="None"){
      if (input$data=="Total Arrests"){
        MapData=data$arrests_total[seq(input$year-2009,dim(data)[1],9)]
        MapDataScl=data$arrests_total
        Caption=paste(input$year,input$data)
      } else 
        if (input$data=="White Only Arrests"){
          MapData=data$arrests_W[seq(input$year-2009,dim(data)[1],9)]
          MapDataScl=data$arrests_W
          Caption=paste(input$year,input$data)
        } else {
          MapData=data$arrests_B[seq(input$year-2009,dim(data)[1],9)]
          MapDataScl=data$arrests_B
          Caption=paste(input$year,input$data)
        }} else
          if (input$adj=="Standardized Incidence Ratio"){
            if (input$data=="Total Arrests"){
              MapData=data$SIRtotal[seq(input$year-2009,dim(data)[1],9)]
              MapDataScl=data$SIRtotal
              Caption=paste(input$year,input$data)
            } else 
              if (input$data=="White Only Arrests"){
                MapData=data$SIRwhite[seq(input$year-2009,dim(data)[1],9)]
                MapDataScl=data$arrests_W_sir
                Caption=paste(input$year,input$data)
              } else {
                MapData=data$SIRblack[seq(input$year-2009,dim(data)[1],9)]
                MapDataScl=data$arrests_B_sir
                Caption=paste(input$year,input$data)
              }        
          } else 
            if (input$adj=="Poisson Regression") {
              if (input$data=="Total Arrests"){
                MapData=exp(data$sptRE[1:44+44*(input$year-2010)])#indexing for these is sorted by year then tract
                MapDataScl=exp(data$sptRE)
                Caption=paste(input$year,input$data)
              } else 
                if (input$data=="White Only Arrests"){
                  MapData=exp(data$sptRE_W[1:44+44*(input$year-2010)])
                  MapDataScl=exp(data$sptRE_W)
                  Caption=paste(input$year,input$data)
                } else {
                  MapData=exp(data$sptRE_B[1:44+44*(input$year-2010)])
                  MapDataScl=exp(data$sptRE_B)
                  Caption=paste(input$year,input$data)
                }        
            } else 
              if (input$adj=="As a Percent of the Population"){
                if (input$data=="Total Arrests"){
                  MapData=(data$arrests_total[seq(input$year-2009,dim(data)[1],9)]+.1)/(data$ct_pop[seq(input$year-2009,dim(data)[1],9)]+.1)*100
                  MapDataScl=(data$arrests_total+.1)/(data$ct_pop+.1)*100
                  Caption=paste(input$year,input$data)
                } else 
                  if (input$data=="White Only Arrests"){
                    MapData=(data$arrests_W[seq(input$year-2009,dim(data)[1],9)]+.1)/(data$ct_white[seq(input$year-2009,dim(data)[1],9)]+.1)*100
                    MapDataScl=(data$arrests_W+.1)/(data$ct_white+.1)*100
                    Caption=paste(input$year,input$data)
                  } else {
                    MapData=(data$arrests_B[seq(input$year-2009,dim(data)[1],9)]+.1)/(data$ct_black[seq(input$year-2009,dim(data)[1],9)]+.1)*100
                    MapDataScl=(data$arrests_B+.1)/(data$ct_black+.1)*100
                    Caption=paste(input$year,input$data)
                  }      
              } else {
                if (input$data=="Total Arrests"){
                  MapData=(data$arrests_total[seq(input$year-2009,dim(data)[1],9)]+.1)/(data$arrests_total[seq(input$year-2009,dim(data)[1],9)]+.1)*100
                  MapDataScl=seq(0,100,.1)
                  Caption=paste(input$year,input$data)
                } else 
                  if (input$data=="White Only Arrests"){
                    MapData=(data$arrests_W[seq(input$year-2009,dim(data)[1],9)]+.1)/(data$arrests_total[seq(input$year-2009,dim(data)[1],9)]+.1)*100
                    MapDataScl=seq(0,100,.1)
                    Caption=paste(input$year,input$data)
                  } else {
                    MapData=(data$arrests_B[seq(input$year-2009,dim(data)[1],9)]+.1)/(data$arrests_total[seq(input$year-2009,dim(data)[1],9)]+.1)*100
                    MapDataScl=seq(0,100,.1)
                    Caption=paste(input$year,input$data)
                  }      
              }
    
    fillmap2(NHtracts,Caption,MapData,map.lty = 0,leg.loc = "beside",y.scl = MapDataScl,leg.rnd = 2)
  })
  
   output$table <- renderTable({
     if (input$adj=="Poisson Regression"){
       if (input$data=="Total Arrests"){
         mat=exp(fe[1:6,])
         colnames(mat)<-c("Mean","95% CI Lower Bound","95% CI Upper Bound")
         rownames(mat)<-c("% Black",
                          "% Living in Poverty",
                          "% Bachelors degree or more",
                          "% Male",
                          "% Secondary Homes",
                          "% Aged 18-24")
         mat
       } else
         if (input$data=="White Only Arrests"){
           mat=exp(fe[7:12,])
           colnames(mat)<-c("Mean","95% CI Lower Bound","95% CI Upper Bound")
           rownames(mat)<-c("% Black",
                            "% Living in Poverty",
                            "% Bachelors degree or more",
                            "% Male",
                            "% Secondary Homes",
                            "% Aged 18-24")
           mat
         } else {
           mat=exp(fe[13:18,])
           colnames(mat)<-c("Mean","95% CI Lower Bound","95% CI Upper Bound")
           rownames(mat)<-c("% Black",
                            "% Living in Poverty",
                            "% Bachelors degree or more",
                            "% Male",
                            "% Secondary Homes",
                            "% Aged 18-24")
           mat
         }
     }},rownames=T,colnames=T,digits=3,width="100%")
})

