
source("R/PrepareDataForFigures.R")

# Figure1 -----------------------------------------------------------------
plot1<-pco.plot(map,species,color = "Type", shape="Location" ,plot.title = "Species",show.legend = FALSE)
plot2<-pco.plot(map.denver,species.denver,color = "Type", shape="Location" ,plot.title = "Species, Denver",show.legend = FALSE)
plot3<-pco.plot(map.unc,species.unc,color = "Type", shape="Location" ,plot.title = "Species, UNC",show.legend = FALSE)
plot4<-pco.plot(map,pathway,color = "Type", shape="Location" ,plot.title = "Pathways",show.legend =FALSE)
plot5<-pco.plot(map.denver,pathway.denver,color = "Type", shape="Location" ,plot.title = "Pathways, Denver",show.legend = FALSE)
plot6<-pco.plot(map.unc,pathway.unc,color = "Type", shape="Location" ,plot.title = "Pathways, UNC",show.legend = FALSE)

pdf("./output/manuscript_figures/Figure1.pdf",width = 10,height = 7)
grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,nrow=2,ncol=3)
dev.off()

pco.plot(map.denver,genus.denver,color = "Type", shape="Location" ,plot.title = "Species, Denver",show.legend = FALSE)
