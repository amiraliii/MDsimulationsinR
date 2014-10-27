multiColnames <- function(...,cname){
    # reads in the argument names are the objects
    objnames <- c(...)
    # if all passes, then we go on with assignment of colnames
    lapply(objnames,function(x) {obj <- get(x);colnames(obj) <- cname;assign(x=x,obj,envir = globalenv())})
    return(invisible())
    }
