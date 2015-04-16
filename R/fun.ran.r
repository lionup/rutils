
#make footnote for ggplot2
#png("mwba_gdp.png"); print(g); makeFootnote("xxx");dev.off()
makeFootnote <- function(footnoteText=
                         format(Sys.time(), "%d %b %Y"),
                         size= 1, color= grey(.1))
{
   require(grid)
   pushViewport(viewport())
   grid.text(label= footnoteText ,
             x = unit(1,"npc") - unit(2, "mm"),
             y= unit(2, "mm"),
             just=c("right", "bottom"),
             gp=gpar(cex= size, col=color))
   popViewport()
}

#handeling missing value using last varialbe
imp <- function (a){
  missing <- is.na(a)
  anymis <- any(missing)
  imputed <- a
  while(anymis){
    imputed[tail(which(missing==TRUE),1)] <- imputed[tail(which(missing==TRUE),1)+1]
    missing <- is.na(imputed)
    anymis <- any( missing )
  }
  return (imputed)
}

