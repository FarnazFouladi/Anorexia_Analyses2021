#Save plots in a list
require(gridExtra)

save.plots<-function(plot.list,num.plots.perpage=9,file.name){
  
  if(!is.null(plot.list)){
    pdf(file.name,width = 11,height = 11)
    pages=floor(length(plot.list)/num.plots.perpage)
    index.last.page=pages*num.plots.perpage
    range<-seq(1,length(plot.list),9)
    lapply(range,function(x){ 
      if(x<=index.last.page){
        grid.arrange(plot.list[[x]],plot.list[[x+1]],plot.list[[x+2]],
                     plot.list[[x+3]],plot.list[[x+4]],plot.list[[x+5]],
                     plot.list[[x+6]],plot.list[[x+7]],plot.list[[x+8]],ncol=3,nrow=3)
      }else{
        
        num.remaining.plots<-length(plot.list)-index.last.page
        
        if(num.remaining.plots==1)
          grid.arrange(plot.list[[x]],ncol=3,nrow=3)
        else if (num.remaining.plots==2)
          grid.arrange(plot.list[[x]],plot.list[[x+1]],ncol=3,nrow=3)
        else if (num.remaining.plots==3)
          grid.arrange(plot.list[[x]],plot.list[[x+1]],plot.list[[x+2]],ncol=3,nrow=3)
        else if (num.remaining.plots==4)
          grid.arrange(plot.list[[x]],plot.list[[x+1]],plot.list[[x+2]],
                       plot.list[[x+3]],ncol=3,nrow=3)
        else if (num.remaining.plots==5)
          grid.arrange(plot.list[[x]],plot.list[[x+1]],plot.list[[x+2]],
                       plot.list[[x+3]],plot.list[[x+4]],ncol=3,nrow=3)
        else if (num.remaining.plots==6)
          grid.arrange(plot.list[[x]],plot.list[[x+1]],plot.list[[x+2]],
                       plot.list[[x+3]],plot.list[[x+4]],plot.list[[x+5]],ncol=3,nrow=3)
        else if (num.remaining.plots==7)
          grid.arrange(plot.list[[x]],plot.list[[x+1]],plot.list[[x+2]],
                       plot.list[[x+3]],plot.list[[x+4]],plot.list[[x+5]],
                       plot.list[[x+6]],ncol=3,nrow=3)
        else if (num.remaining.plots==8)
          grid.arrange(plot.list[[x]],plot.list[[x+1]],plot.list[[x+2]],
                       plot.list[[x+3]],plot.list[[x+4]],plot.list[[x+5]],
                       plot.list[[x+6]],plot.list[[x+7]],ncol=3,nrow=3)
        
      }
    })
    dev.off()
    
  } else {
    cat("No plot in the list!\n")
  }
  
}




