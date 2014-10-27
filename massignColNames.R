multiColnames <- function(objlist,cname){
    # reads in the argument names are the objects
    lapply(objlist,function(x) {obj <- get(x);colnames(obj) <- cname;assign(x=x,obj,envir = globalenv())})
    return(invisible())
    }
