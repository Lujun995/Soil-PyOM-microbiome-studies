myfunc12<-function(data,Sample.Type,color,xlim1,xlim2,ylim1,ylim2){
  laymat=matrix(
    c(rep(c(1,1,2,2,2,2,2),5),
      rep(c(NA,NA,3,3,3,3,3),2)),
    nrow=7,byrow = TRUE)
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  data.tmp<-data[c(isBurnt,isMock),]
  Sample.Type.tmp<-Sample.Type[c(isBurnt,isMock)]
  color.tmp<-color[c(2,3)]
  
  ggplot(data.tmp,aes(MDS1,colour=Sample.Type.tmp,fill=Sample.Type.tmp))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(10) + theme(panel.grid = element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+
    scale_fill_manual(values = color.tmp)+
    scale_color_manual(values = color.tmp)+
    xlim(xlim1,xlim2)+labs(x = "NMDS1")->pic3
  
  ggplot(data, aes(MDS1, MDS2, colour = Sample.Type, fill = Sample.Type)) +
    geom_point(size = 1.5,aes(alpha=Sample.Type)) + 
    theme_bw(10) + theme(panel.grid = element_blank(),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         legend.position="none") + 
    scale_color_manual(values = color) + 
    scale_alpha_manual(values=c(0.3,1,1,0.3))+
    scale_fill_manual(values = color) + 
    xlim(xlim1,xlim2)+ylim(ylim1,ylim2)->pic2
  
  ggplot(data.tmp,aes(MDS2,colour=Sample.Type.tmp,fill=Sample.Type.tmp))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(10) + theme(panel.grid = element_blank(),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         axis.text.y=element_text(angle = 90,hjust = 0.5),
                         legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+labs(x = "NMDS2")+
    scale_fill_manual(values = color.tmp)+
    scale_color_manual(values = color.tmp)+coord_flip()+
    xlim(ylim1,ylim2)->pic1
  
  grid.arrange(pic1,pic2,pic3,layout_matrix=laymat)
}

myfunc13<-function(data,Sample.Type,color,xlim1,xlim2,ylim1,ylim2){
  laymat=matrix(
    c(rep(c(1,1,2,2,2,2,2),5),
      rep(c(NA,NA,3,3,3,3,3),2)),
    nrow=7,byrow = TRUE)
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  data.tmp<-data[c(isBurnt,isMock),]
  Sample.Type.tmp<-Sample.Type[c(isBurnt,isMock)]
  color.tmp<-color[c(2,3)]
  
  ggplot(data.tmp,aes(MDS1,colour=Sample.Type.tmp,fill=Sample.Type.tmp))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(10) + theme(panel.grid = element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+
    scale_fill_manual(values = color.tmp)+
    scale_color_manual(values = color.tmp)+
    xlim(xlim1,xlim2)+labs(x = "NMDS1")->pic3
  
  ggplot(data, aes(MDS1, MDS3, colour = Sample.Type, fill = Sample.Type)) +
    geom_point(size = 1.5,aes(alpha=Sample.Type)) + 
    theme_bw(10) + theme(panel.grid = element_blank(),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         legend.position="none") + 
    scale_color_manual(values = color) + 
    scale_alpha_manual(values=c(0.3,1,1,0.3))+
    scale_fill_manual(values = color) + 
    xlim(xlim1,xlim2)+ylim(ylim1,ylim2)->pic2
  
  ggplot(data.tmp,aes(MDS3,colour=Sample.Type.tmp,fill=Sample.Type.tmp))+
    geom_density(adjust=2,alpha=0.5)+
    theme_bw(10) + theme(panel.grid = element_blank(),
                         axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         axis.text.y=element_text(angle = 90,hjust = 0.5),
                         legend.position="none")+
    scale_y_continuous(labels=scaleFUN)+labs(x = "NMDS3")+
    scale_fill_manual(values = color.tmp)+
    scale_color_manual(values = color.tmp)+coord_flip()+
    xlim(ylim1,ylim2)->pic1
  
  grid.arrange(pic1,pic2,pic3,layout_matrix=laymat)
}