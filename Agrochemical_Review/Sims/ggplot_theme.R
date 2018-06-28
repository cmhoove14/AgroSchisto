require(extrafont)

theme_ms <- function(base_size=10, base_family="Times New Roman") {
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family)+
      theme(text=element_text(color="black"),
            plot.title = element_text(face="bold"),
            strip.background = element_blank(),
            strip.text.y = element_text(face = "bold", size = rel(1)),
            axis.title=element_text(size = rel(1)),
            axis.text=element_text(size = rel(1), color = "black"),
            legend.title=element_text(face="bold"),
            #legend.text=element_text(face="bold"),
            legend.background=element_rect(fill="transparent"),
            legend.key.size = unit(0.8, 'lines'),
            panel.border=element_rect(color="black",size=1)
    ))
}
