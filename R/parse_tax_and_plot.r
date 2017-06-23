# This file is part of wheelP.
# 
# wheelP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# wheelP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with wheelP.  If not, see <http://www.gnu.org/licenses/>.

#' Metacoder plot
#'
#' Requires you to load metacoder package, otherwise names not print
#'
#' @author From metacoder paper
#'
#' @importFrom magrittr "%>%"
#'
#' @export
parse_tax_and_plot <- function(file,col,output_file = output_file,
                               n.supertaxa = n.supertaxa,
                               num.changed = num.changed,
                               min_fdr = min_fdr){


  data <- metacoder::parse_taxonomy_table(#input = as.data.frame(bp_res),
    file_path =  file,
    taxon_col = c("class" = col),
    other_col_type = "obs_info",
    class_sep = "\\|")
  data$obs_data$FDR <- as.numeric(data$obs_data$FDR)
  data$taxon_funcs <- c(data$taxon_funcs,
                        change = function(x, subset = NULL) {
                          vapply(metacoder::obs(x),
                                 function(i) {
                                   i <- unlist(i)
                                   #obs_change <- data$obs_data[i, ]$log2FoldChange[data$obs_data[i, ]$padj <= min_p_value]
                                   #obs_change <- data$logFC[i, ]$logFC[data$obs_data[i, ]$padj <= min_p_value]
                                   obs_change <- x$obs_data[i, ]$logFC[x$obs_data[i, ]$FDR <= min_fdr]
                                   mean(obs_change, na.rm = TRUE)
                                 },
                                 numeric(1))
                        },
                        num_changed = function(x, subset = NULL) {
                          vapply(metacoder::obs(x),
                                 function(i) {
                                   i <- unlist(i)
                                   #sum(data$obs_data[i, ]$padj <= min_p_value, na.rm = TRUE)
                                   sum(x$obs_data[i, ]$FDR <= min_fdr, na.rm = TRUE)
                                 },
                                 numeric(1))
                        })

  # Define replacement
  to_replace <- matrix(ncol = 2, byrow = TRUE,
                       c("biological_process",""))


  # Plot
  data %>%
    metacoder::filter_taxa(num_changed > num.changed) %>%
    metacoder::filter_taxa(n_supertaxa <= n.supertaxa) %>%
    metacoder::mutate_taxa(name = gsub("_", " ", name),
                f_change = 2^abs(change) * sign(change)) %>%
    metacoder::mutate_taxa(short_name = vapply(name, FUN.VALUE = character(1), function(x) {
      wheelP:::mgsub(pattern = to_replace[, 1],
                     replacement =  to_replace[, 2], x, fixed = TRUE)})) %>%
    metacoder::heat_tree(#node_label = c("A","B"),
      node_label = ifelse(abs(change) > 0.5, short_name, NA),
      node_size = num_changed,
      # node_size_trans = "log10",
      node_size_range = c(0.008, 0.03),
      # node_label_size_trans = "log10",
      node_label_size_range = c(0.012, 0.02),
      # edge_size_trans = "log10",
      edge_size_range = c(0.008, 0.03) / 2,
      node_color = f_change,
      node_color_trans = "linear",
      node_color_range = metacoder::diverging_palette(),
      node_color_interval = c(-5, 5),
      edge_color_trans = "linear",
      edge_color_range = metacoder::diverging_palette(),
      edge_color_interval =  c(-5, 5),
      node_label_max = 80,
      node_color_axis_label = "log2(Fold change)",
      node_size_axis_label = "Number of genes",
      layout = "da", initial_layout = "re",
      output_file = output_file)

  return(data)
}

#' mgsub
#'
#' Internal
#'
#' @references from: http://stackoverflow.com/questions/15253954/replace-multiple-arguments-with-gsub
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}

