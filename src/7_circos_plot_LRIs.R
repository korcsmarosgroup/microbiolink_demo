# AIM:Visualizing condition specific ligand-receptor connection between cells
# Input: Table describing LRIs
# Output: PNG/PDF file

#Install circlize package if needed
#install.packages("circlize")

#call the package
library(circlize)

#!/usr/bin/env Rscript
require(circlize)

#Import input file
mat <- read.csv('microbiolink_demo/results/intercellular_interactions/Epithelial_2_Healthy_DC_Healthy_LRIs.csv')
rownames(mat) <- mat[,1]
mat <- as.matrix(mat[,-1])

png_name = 'LRI_circos.png'
#pdf_name = 'test.pdf'
#pdf('epi-DC_H.pdf', width = 7, height = 7)
png(png_name, width = 7, height = 7, units = 'in', res = 300)

chordDiagram(mat, annotationTrack = 'grid', preAllocateTracks = 1)

circos.trackPlotRegion(
  track.index = 1,
  ylim = c(0, .3),
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      mean(xlim),
      -ylim[1],
      sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, .5),
      cex = .6
    )
  },
  bg.border = NA
)
circos.clear()
dev.off()