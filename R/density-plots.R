#' Grouped density plots for every (categorical x continuous) column pair
#'
#' For each combination of a categorical grouping column and a continuous
#' column, draws a \pkg{ggpubr} density plot (with mean line and rug) colored
#' by group, optionally saving one PNG per pair. Groups with fewer than two
#' observations for a given continuous column are dropped before plotting.
#'
#' @param PlotDT A data.table (or data.frame) containing the columns named
#'   in \code{ContinuousCols} and \code{CategoricalDependents}.
#' @param ContinuousCols Character vector of continuous column names (the x
#'   axis of each density plot).
#' @param CategoricalDependents Character vector of categorical column names
#'   (the grouping/fill variable of each density plot).
#' @param OutPath Directory to save PNGs into when \code{SavePlots = TRUE}.
#' @param SavePlots Logical; if \code{TRUE} (default), write one PNG per
#'   pair to \code{OutPath} and return \code{NULL}. If \code{FALSE}, return
#'   the list of plot objects instead of saving.
#' @return If \code{SavePlots = FALSE}, a named list of ggplot objects (one
#'   per categorical x continuous pair); otherwise \code{NULL} (called for
#'   its side effect of writing PNGs).
#' @examples
#' \dontrun{
#' load(system.file("extdata", "dfMigrant.rda", package = "dataMineR"))
#' DensityPlotCompile(dfMigrant, ContinuousCols = "Age",
#'                     CategoricalDependents = "Migrant",
#'                     OutPath = tempdir())
#' }
#' @export
DensityPlotCompile <- function(PlotDT, ContinuousCols, CategoricalDependents, OutPath, SavePlots = TRUE){
  CompletePlotList <- list()
  PlotDT <- copy(as.data.table(PlotDT))
  missing_cols <- setdiff(c(ContinuousCols, CategoricalDependents), names(PlotDT))
  if(length(missing_cols) > 0L){ stop("Columns not found in PlotDT: ", paste(missing_cols, collapse = ", ")) }
  if(SavePlots && !dir.exists(OutPath)){ dir.create(OutPath, recursive = TRUE, showWarnings = FALSE) }
  colnames(PlotDT) <- gsub(";", "_", colnames(PlotDT))
  if (anyDuplicated(names(PlotDT))){ stop("Column-name sanitization created duplicate names; rename the input columns first.") }
  ContinuousColsSub <- gsub(";", "_", ContinuousCols)
  CategoricalDependentsSub <- gsub(";", "_", CategoricalDependents)
  for(a in seq_along(CategoricalDependentsSub)){
    current_dep <- CategoricalDependentsSub[a]
    plot_sublist <- list()
    for(i in seq_along(ContinuousColsSub)){
      current_cont <- ContinuousColsSub[i]
      if(current_cont %in% colnames(PlotDT)){
        clean_dt <- PlotDT[!is.na(get(current_cont)) & !is.na(get(current_dep))]
        group_counts <- clean_dt[, .N, by = c(current_dep)]
        valid_groups <- group_counts[N >= 2, get(current_dep)]
        filtered_dt <- clean_dt[get(current_dep) %in% valid_groups]
        if(nrow(filtered_dt) == 0) next
        filtered_dt[[current_dep]] <- factor(filtered_dt[[current_dep]])
        num_colors <- length(unique(filtered_dt[[current_dep]]))
        custom_palette <- if(num_colors <= 2) {
          c("#00AFBB", "#E7B800")[1:num_colors]
        } else {
          hue_pal()(num_colors)
        }
        p_name <- paste(current_dep, current_cont, sep = "_")
        plot_sublist[[p_name]] <- ggdensity(filtered_dt, 
                                            x = current_cont, 
                                            add = "mean", 
                                            rug = TRUE, 
                                            color = current_dep, 
                                            fill = current_dep,
                                            palette = custom_palette )
        if(SavePlots){
          ggplot2::ggsave(file.path(OutPath, paste0(p_name, ".png")), plot = plot_sublist[[p_name]], width = 5, height = 5, dpi = 100)
        }
      } 
    }
    if(!SavePlots){ CompletePlotList <- c(CompletePlotList, plot_sublist) }
  }
  if(!SavePlots){
    return(CompletePlotList)
  }
}
