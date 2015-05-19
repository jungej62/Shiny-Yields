setwd("~/Documents/Projects/PolycultureBioenergy")
source.with.encoding('~/Documents/Projects/HelperFunctions.R', encoding='UTF-8')
pdat<-read.table(file="PolyDat.txt", header=T)
pdat$Location<-as.factor(pdat$Location); pdat$fYear<-as.factor(pdat$Year);pdat$Plot<-as.factor(pdat$Plot);
pdat$Trt<-as.factor(pdat$Trt); pdat$Nfert<-as.factor(pdat$Nfert)
TrtNames<-c("1 spp: Switchgrass", "1 spp: Big Bluestem", "1 spp: Indiangrass", "1 spp: Canada Wild Rye",
            "4 spp: Grasses", "4 spp: Legumes", "4 spp: Forbs", "8 spp: Grass + Legumes",
            "8 spp: Grass + Forbs", "8 spp: Legumes + Forbs", "12 spp: Mix", "24 spp: Mix")
pdat$TotalBiomass<-pdat$TotalBiomass*2

shinyServer(function(input, output) {
  
  

  # Generate a plot of the requested variable against mpg and
  # only include outliers if requested
  output$yldPlot <- renderPlot({
    ggplot(summarySE(subset(pdat, pdat$Trt==input$variable),
                     measurevar="TotalBiomass", groupvars="Site", na.rm=T),
           aes(x=Site, y=TotalBiomass))+
      geom_bar(position=position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=TotalBiomass-se, ymax=TotalBiomass+se),
                    width=.2,position=position_dodge(.9))+
      xlab(expression("Site " ))+
      ylab(expression("Biomass Yield " ~ (Mg ~ ha^{-1})))
    
  })
})