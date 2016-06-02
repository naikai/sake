#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
create.label.color <- function(data, color)
{
     groupCodes <- as.factor(data)
          colorCodes <- color(length(levels(groupCodes)))
            color.idx <- match(groupCodes, levels(groupCodes))
              label.color <- colorCodes[color.idx]
                return(label.color)
}


