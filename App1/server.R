library(gstat); library(sp); library(maptools); library(ggplot2); library(scales);library(ggmap)
library(GISTools);library(devtools); library(shiny); library(maps); library(reshape)
source_gist("https://gist.github.com/jungej62/5a76b72bcd0b12a7ce8a")
pdat<-read.table(file="PolyDat_cpy.txt", header=T)
pdat$Location<-as.factor(pdat$Location); pdat$fYear<-as.factor(pdat$Year);pdat$Plot<-as.factor(pdat$Plot);
pdat$Trt<-as.factor(pdat$Trt); pdat$Nfert<-as.factor(pdat$Nfert)
TrtNames<-c("1 spp: Switchgrass", "1 spp: Big Bluestem", "1 spp: Indiangrass", "1 spp: Canada Wild Rye",
           "4 spp: Grasses", "4 spp: Legumes", "4 spp: Forbs", "8 spp: Grass + Legumes",
           "8 spp: Grass + Forbs", "8 spp: Legumes + Forbs", "12 spp: Mix", "24 spp: Mix")
pdat$TotalBiomass<-pdat$TotalBiomass*2
pdat2<-pdat[,c(2:8,12)]
pdat2$Site2<-pdat2$Site
pdat2$Site2<-factor(pdat2$Site, labels=c("Becker", "Crookston", "Fargo", "Lamberton", "Mahnoman", "Roseau", "Red Lake Falls", "Saint Paul", "Waseca"))
dd1<-cast(pdat2, Lat+Lon+Site2~Trt, value="TotalBiomass", mean, na.rm=T)
dd1<-dd1[order(dd1$Site2),]


shinyServer(function(input, output) {
  selectedData<-reactive({
    data.frame(dd1[1:3], "Yld"=dd1[,input$trt])
  })
  output$yldPlot <- renderPlot({
  ggmap(get_googlemap(center=c(lon = -93.920710, lat = 46.52), zoom=6, maptype="hybrid"))+
    geom_point(data=selectedData(), aes(x=Lon, y=Lat, color=Yld), size=6)+
    geom_text(data=selectedData(), aes(x=Lon+c(0.8,-0.2,0.8,0.8,1.2,0.9,1.5,0.8,0.8),
                                          y=Lat+c(0.2, -0.2, -0.1, 0.3, -0.2, -0.2, 0.2, -0.25, -0.25),
                                          label=Site2), size=6,  fontface="bold")+
    scale_color_gradient2(low = "red", mid = "white", high = "blue", 
                          midpoint = mean(subset(pdat, Trt==input$trt)$TotalBiomass, na.rm=T), space = "rgb", na.value = "grey50", guide = "colourbar")
    
  })
})