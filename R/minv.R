minv <- function(x, m) tapply(x, 1:length(x), function(y) min(m, y))
