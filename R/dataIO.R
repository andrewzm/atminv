#' Check md5 sum in cache directory
#'
#' @description This function takes an md5 and checks whether it is in the default cache directory ~/cache/.
#' @param md5 the md5 checksum
#' @param path the directory where the cache is stored
#' @return True or False, depending on whether the file exists in the cache.
#' @export
#' @examples
#' require(digest)
#' md5 <- digest(2)
#' check_md5(md5,".")
check_md5 <- function(md5,path) {
    return(file.exists(file.path(path,paste0(md5,".rda"))))
}

#' Generate md5 checksum from function
#'
#' Creates an md5 checksum of a function by dumping it and then digesting the file.
#' @param fn the function to digest
#' @param path the directory where dumping will be carried out
#' @details Note that a ~/cache folder needs to be present for this to work.
#' @keywords md5, digest
#' @export
#' @examples
#' myfun1 <- function(x) x^2
#' x <- md5_fun(myfun1,".")
md5_fun <- function(fn,path) {
    path <- file.path(path,"fun_check.R")
    dump(ls(pattern = 'fn'),path)
    #save(fn,file=path)
    fun_digest <- digest::digest(file=path)
    unlink(path)
    return(fun_digest)
}

#' md5 checksum function
#'
#' Creates a wrapper for generating the md5 checksum for the calling function and its arguments. If a file with the md5 as its filename exists in the specified folder then this is loaded. Otherwise the function is evaluated and the results are stored in the specified folder.
#' @param path the path to where the files will be cached. 
#' @return An md5 wrapper for intensive computations. This has arguments
#' \itemize{
#' \item \code{fn}: the function being called
#' \item \code{...}: the arguments to be passed to the function
#' \item \code{print_output}: if T, details of md5 digest are outputted to screen.
#' }
#' @details md5 checksums on the called function and any function-arguments are generated by creating a text file in the specified folder and digesting 
#' the text file before deleting it. md5 checksums on character arguments are carried out on the file (if the file exists) or on the character as appropriate.
#' @keywords md5, digest
#' @export
#' @examples
#' myfun1 <- function(x) x^2
#' myfun2 <- function(y,f) y^3 + f(y)
#' md5_wrapper <- md5_cache(".")
#' x <- md5_wrapper(myfun2,2,myfun1)
md5_cache <- function(path) {
    stopifnot(is.character(path))
    dir.create(file.path(path), showWarnings = FALSE)
    md5 <- NA
    
    function(fn,...,print_output=F) {
        args <- list(...)
        if (length(args)==0) stop("Nothing to check md5 against")
        md5 <<- md5_fun(fn,path)
        lapply(args,function(x) {  
            if(is.character(x) & length(x) == 1) {
                if(file.exists(x)) {
                    md5_temp <- digest(file=x)
                } else {
                    md5_temp <- digest(x)
                }
            } else {
                if(class(x) == "function") {
                    md5_temp <- md5_fun(x,path)
                } else {
                    md5_temp <- digest(x)
                }
            }
            md5 <<- digest(c(md5,md5_temp))
        })
        
        if(check_md5(md5,path)) {
            if(print_output) cat("Found existing MD5 for data. Loading...",sep="\n")
            load(file=file.path(path,paste(md5,".rda",sep=""))) 
        } else {
            cat("Data not found. Evaluating expression...",sep="\n")
            flush.console()
            X <- fn(...)
            save(X,file=file.path(path,paste0(md5,".rda")))
        }
        if(print_output) cat(paste("md5 for this result is ",md5,sep=""),sep="\n")
        return(X)  
    }
}

#' Polygon-attribution
#' 
#' Takes a data frame \code{df} with fields \code{x} and \code{y} and a shape in data frame format, \code{shape_table}, and determines whether, and in which, polygon each point in \code{df} lies.
#' @param df data frame to which the polygon id is attributed. Must contain fields \code{x} and \code{y}.
#' @param shape_df the polygon-containing data frame. Must contain fields \code{x}, \code{y} and \code{id}. Multiple polygons/land-masses can be treated simultaneously through \code{id}.
#' @return a vector of length \code{nrow(df)} indicating the id in which each point in \code{df} lies.
#' @export
#' @examples
#' df <- data.frame(x=c(0.2,1),y=c(0.2,1))
#' shape_df <- data.frame(x=c(0,0,0.5,0.5),
#' y = c(0,0.5,0.5,0),
#' id = 1)
#' df$id <- attribute_polygon(df,shape_df)
attribute_polygon <- function(df,shape_df) {
    
    df$n <- 1:nrow(df)
    vals <- rep(0,nrow(df))
    
    if (length(unique(shape_df$id)) == 1) {
        df2 <- subset(df,x > min(shape_df$x) & x < max(shape_df$x) & y > min(shape_df$y) & y < max(shape_df$y)) # Find enclosing box
        myind <- df2[which(pnt.in.poly(cbind(df2$x,df2$y),shape_df[c("x","y")])$pip == 1),]$n
        vals[myind] <- T
    } else {
        for(i in unique(shape_df$id)){
            my_sub <- subset(shape_df,id==i)
            df2 <- subset(df,x > min(my_sub$x) & x < max(my_sub$x) & y > min(my_sub$y) & y < max(my_sub$y)) # Find enclosing box
            if(nrow(df2) > 0) {
                myind <- df2[which(pnt.in.poly(cbind(df2$x,df2$y),my_sub[c("x","y")])$pip == 1),]$n
                vals[myind] <- i
            }
        }
        
    }
    return(vals)
}
