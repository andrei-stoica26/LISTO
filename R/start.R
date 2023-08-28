SeuratsList <- function(files, conditions, minCells = 10)
  lapply(1:length(files), function(i){
    message(str_c("Adding Seurat object: ", conditions[i])," from file: ", files[i], ".")
    CreateSeuratObject(counts = readRDS(files[i]), project = conditions[i], min.cells = minCells)
  })

CountCells <- function(seuratObj, type){
  if (is.null(seuratObj))WriteMessage("0", type)
  else{
    grouped_cells <- aggregate(nFeature_RNA ~ orig.ident, data = seuratObj@meta.data, length)[, c(2,1)]
    WriteMessage(paste(apply(grouped_cells, 1, paste,collapse=" "), collapse = " cells, "), type)
  }
}

WriteMessage <- function(result, type){
  message(ifelse (type == "in",
                  str_c("There are now ", result, " cells in this Seurat object."),
                  str_c(result, " cells were removed.")))
}

MergeSeurats <- function(files, conditions, ids, minCells = 10){
  seurats <- SeuratsList(files, conditions, minCells)
  seuratObj <- merge(seurats[[1]],seurats[c(2:length(files))],add.cell.ids = ids, project = "All conditions")
  message("Merged Seurat objects.")
  message(CountCells(seuratObj, "in"))
  message("Finished.")
  return (seuratObj)
}
