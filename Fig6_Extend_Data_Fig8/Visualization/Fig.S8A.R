rm(list=ls())

library(ComplexUpset)
library(ggplot2)
library(xlsx)

Data_file = "Upsetplot_raw_data.xlsx"
example = read.xlsx(Data_file,sheetIndex = 1)
toplot= c('Clinical.stage.Subtype',	'Transcriptomic.data', 'Pathologic',	'Metabolomic','Radiologic')
example = example[,toplot]


upset(example, toplot,
      mode = "inclusive_intersection",
      #base_annotations=list('Size'=(intersection_size(counts=T))),
      sort_sets='descending',
      sort_intersections_by = c('cardinality', 'degree'), sort_intersections='descending',
      name = 'toplot', 
      intersections='all',
      width_ratio = 0.1 
)
