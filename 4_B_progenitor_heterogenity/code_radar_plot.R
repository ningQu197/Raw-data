# https://www.jianshu.com/p/ba8ac27281da
devtools::install_github("ricardo-bion/ggradar", dependencies=TRUE)
library(ggradar)

data = read.table("./GOBP_radar_plot.txt", header=T, sep="\t")
# 
colors = c("tomato", "tomato", "tomato", "tomato", "tomato")
p = ggradar(data, grid.min=0, grid.mid=0.5, grid.max=1, background.circle.colour="white")+scale_color_manual(values=colors)
ggsave("GOBP_radar_plot.txt", p, width=6, height=6)




