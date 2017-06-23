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

#' From metacoder
#'
#' @author From Metacoder paper
#'
#' @export
term_class <- function(x,current = x, all_paths = TRUE,
                       type = GO.db::GOCCPARENTS, verbose = TRUE,
                       valid_relationships = c("is_a")) {
  # x <- dat$go[1]
  # all_paths <- FALSE
  # type <- GO.db::GOBPPARENTS
  # verbose <- TRUE
  # valid_relationships <- c("is_a")
  # as.list(type[x[1]])[[1]]
  # as.list(type[[x[1]]])[[1]]

  # Get immediate children of current taxon
  parents = tryCatch({
    possible_parents <- AnnotationDbi::as.list(type[x[1]])[[1]]
    if (! is.null(valid_relationships)) {
      possible_parents <- possible_parents[names(possible_parents) %in% valid_relationships]
    }
    names(AnnotationDbi::Term(possible_parents))
  }, error = function(e) {
    c()
  })

  # only go down one path if desired
  if (! all_paths) {
    parents <- parents[1]
  }
  parents <- parents[parents != "all"]

  if (is.null(parents)) {
    return(c())
  } else if (length(parents) == 0) {
    cat(length(x))
    return(paste0(collapse = "|", AnnotationDbi::Term(x)))
  } else {
    next_x <- lapply(parents, function(y) c(y, x))

    # Run this function on them to get their output
    child_output <- lapply(next_x, term_class,
                           all_paths = all_paths, type = type)
    output <- unlist(child_output, recursive = FALSE)

    return(output)
  }
}
