maxv <- function(x, m) tapply(x, 1:length(x), function(y) max(m, y))
