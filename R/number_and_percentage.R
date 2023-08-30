#' @export

ColumnToString <- function(x)
  as.character(enquo(x))[2]


ColStrToDesc <- function (col_name){
  return (switch(col_name,
                 "orig.ident" = "Condition",
                 "seurat_clusters" = "Cluster",
                 "Phase" = "Cell cycle phase",
                 "cell.type" = "Cell type",
                 "mt.band" = "Mitochondrial percentage",
                 "shannon.band" = "Shannon diversity score",
                 "CSC.band" = "Cancer cell stemness score",
                 "variable.CSC.band" = "Cancer cell stemness (variable genes) score",
                 "important.CSC.band" = "Cancer cell stemness (important markers) score",
                 "pseudotime.band" = "Pseudotime",
                 "test.band" = "Test score",
                 "activity.band" = "Activity score"

  ))
}

XYBarPlot <- function (piv, fill_col, x_col, x_lab, y_lab, legend_lab)
  ggplot(piv, aes(fill={{fill_col}}, y=n, x={{x_col}})) +
  geom_bar(stat="identity") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(fill= legend_lab) +
  guides(col = guide_legend(ncol = 2)) +
  #scale_fill_manual(values=wes_palette("Darjeeling1")) + ######allowed when the number of identities is <=5
  theme_bw() +
  theme(
    legend.position = 'right',
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x=element_text(angle=30, hjust=1)
  )


Merge2 <- function(list1, list2)
  mapply(function(x, y) list(x, y), list1, list2)

LowerCase <- function(desc){
  if (desc == "Shannon diversity score") return (desc)
  return (tolower(desc))
}


XYWrite <- function(df, xdesc, ydesc, cell.quantity){
  shared.substring <- str_c(cell.quantity, " of cells from each ", LowerCase(xdesc)," per ", LowerCase(ydesc), ".")
  WriteNumberedCSV(str_c(shared.substring, "csv"), df, AddCSVNumber(1))
  message(str_c("Wrote the ", LowerCase(shared.substring)))
}

GetLevels <- function(column){
  # Sorting at this stage applies to annotated cell types, which would otherwise get sorted and cause a discrepancy
  # between the order in the number figure and the order in the percentage figure
  if (is.null(levels(column))) return (sort(unique(column)))
  return (levels(unique(column)))
}


XYExplore <- function(seuratObj, x, y){
  # Convert Seurat columns to strings
  xcol <- ColumnToString({{x}})
  ycol <- ColumnToString({{y}})
  xdesc <- ColStrToDesc(xcol)
  ydesc <- ColStrToDesc(ycol)

  # Plot number of cells from X per Y
  piv <- dplyr::count(seuratObj@meta.data, {{x}}, {{y}})
  plot1 <- XYBarPlot (piv, fill_col = {{y}}, x_col = {{x}}, x_lab = xdesc, y_lab = "Number of cells", legend_lab = ydesc)

  # Save the correct order of levels to counter their undesired erasure by pivot_longer
  orderedXLevels <- GetLevels(piv[[xcol]])
  orderedYLevels <- GetLevels(piv[[ycol]])

  # Write number of cells from X per Y
  counts.df <- spread(piv,{{y}}, n)
  counts.df <- replace(counts.df, is.na(counts.df), 0)
  row.names(counts.df) <- counts.df[,1]
  counts.df[1] <- NULL
  counts.df <- adorn_totals(counts.df, "col", ... = 1:ncol(counts.df))
  XYWrite(counts.df, xdesc, ydesc, "Number")

  # Write percentage of cells from X per Y
  percentages.df <- counts.df/counts.df[,"Total"] * 100
  percentages.df <- rbind(percentages.df, sqrt(matrixStats::colVars(as.matrix(percentages.df), std = TRUE)))
  rownames(percentages.df)[length(rownames(percentages.df))] <- "Std"
  XYWrite(percentages.df, xdesc, ydesc, "Percentage")

  # Remove standard deviation row
  percentages.df <- percentages.df[1:nrow(percentages.df)-1,]

  # Make a column for row names. This is necessary for pivoting
  setDT(percentages.df, keep.rownames = TRUE)

  # Pivot the percentage matrix. This is necessary for plotting
  piv <- pivot_longer(dplyr::select(percentages.df, -Total), cols = -rn, names_to = ycol, values_to = "n")

  # Recover the correct order of labels (possibly lost after applying pivot_longer)
  piv$rn <- factor(piv$rn, levels = orderedXLevels)
  piv[[ycol]] <- factor(piv[[ycol]], levels = orderedYLevels)

  # Plot percentage of cells from X per Y
  plot2 <- XYBarPlot (piv, fill_col = {{y}}, x_col = rn, x_lab = xdesc, y_lab = "Percentage of cells",legend_lab = ydesc)
  return(list(plot1, plot2))
}

NumberAndPercentage <- function(seuratObj, x, y){
  forward <- XYExplore(seuratObj, {{x}}, {{y}})
  reverse <- XYExplore(seuratObj, {{y}}, {{x}})
  xstr <- LowerCase(ColStrToDesc(ColumnToString({{x}})))
  ystr <- LowerCase(ColStrToDesc(ColumnToString({{y}})))
  caption <- str_c("**A)** Number of cells from each ", ystr, " in each ", xstr, "; ",
                   " **B)** Number of cells from each ", xstr, " in each ", ystr, ";<br>",
                   " **C)** Percentage of cells from each ", ystr, " in each ", xstr, ";",
                   " **D)** Percentage of cells from each ", xstr, " in each ", ystr, "." )
  PDFSlide(str_c("Numbers and percentages of cells for ", xstr, "s and ", ystr, "s.pdf"),
           GrobPlot(Merge2(forward, reverse), 2, caption), AddFileNumber(1))
}

XYAnalysis <- function(seuratObj, ...){
  params <- unlist(lapply(as.list(match.call()), ColumnToString))
  colnames(params) <- NULL
  params <- params[3:length(params)]
  combs <- unlist(combn(params, 2, simplify = F))
  invisible(lapply(seq(1, (length(combs) - 1), 2),
                   function (x) NumberAndPercentage(seuratObj, !!as.name(combs[x]), !!as.name(combs[x+1]))))
}
