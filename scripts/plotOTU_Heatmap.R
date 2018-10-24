#'
#' plot a heatmap of samples based on relative abundance of OTUs
#'
#' This function consumes an OTU table and a rank, as well as some optional parameters, 
#' and creates a heatmap showing the abundance of the OTUs at the given taxonomic rank 
#' for each sample

library("gplots")
library("RColorBrewer")


############################################################################################
.valid.data <- function(data, is.OTU=TRUE, 
                        export.sampleIDs=FALSE) {
  if ( class(data) != "list" ) {
      stop("please provide ecology data sets as list, see ?RAM.input.formatting for details")
  } 
  if ( length(data) != length(names(data)) ) {
        stop("each datasets in the list should have a given name, see ?RAM.input.formatting for details")
  }
  samples <- list()
  for ( i in 1:length(data) ) {
    label <- names(data)[i]
    elem <- data[[i]]
    if ( is.null(elem) ) { break }
    if ( is.OTU ) {
      valid.OTU(elem)
      samples[[label]] <- names(elem)[-ncol(elem)]
    } else {
      samples[[label]] <- rownames(elem)
      if ( "taxonomy" %in% colnames(elem) ) {
         stop("are you sure this is NOT an OTU table?")
      }
    }
  }
  sampleids <- samples[[1]]
  if ( length(data) != 1 ) {
    for ( i in 2:length(data) ) {
      if ( !identical(sampleids, samples[[i]]) ) {
          stop("samples in each data set don't match, check sample# and order in your data sets before performing any analysis")
      }
    }
  } 
  if (export.sampleIDs) { return(sampleids) }
}


valid.OTU <- function(otu1, otu2=NULL) {
  
  given.both <- !is.null(otu2)

  # tests for otu1
  # ensure object is a data.frame
  if (!is.data.frame(otu1)) {
    stop("the given object for otu1 is not a data frame.")
  }
  
  tax <- dim(otu1)[2]
  
  # ensure the last column contains taxonomic data
  if (names(otu1)[tax] != "taxonomy") {
    stop("no 'taxonomy' column detected in otu1. Please ensure the last column of the table is titled 'taxonomy' (capitalization matters) and try again.")
  }
  
  # check if all data is numeric other than tax columns
  if (!all(apply(otu1[ ,-tax, drop=FALSE], is.numeric, MARGIN=1))) {
    stop("OTU data other than taxonomic information is not numeric in otu1. Please fix the data and try again.")
  }
  
  # check that all counts are non-negative
  if (any(otu1[ , -tax] < 0)) {
    stop("negative counts detected in otu1; cannot continue. Please ensure all counts are positive and try again.")
  }
  
  if (any(colSums(otu1[ ,-tax, drop=FALSE]) == 0)) {
    warning("some samples in otu1 have zero counts. Some functions may behave bizarrely as a result.")
  }
  
  missing.tax <- which(otu1[ , tax, drop=FALSE] == "")
  
  if (length(missing.tax) != 0) {
    missing.tax.percent <- 100 * (length(missing.tax) / dim(otu1)[1]) 
    warning(paste("taxonomic data is missing for ", format(missing.tax.percent, digits=3),
                  "% of OTUs in otu1.", sep=""))
  }
  
  # tests for otu2 data
  if (given.both) {
    
    if (!identical(names(otu1), names(otu2))) {
      stop("the samples for otu1 and otu2 do not match. Please ensure you have matching otu1 and otu2 data.")
    }
    
    # ensure object is a data.frame
    if (!is.data.frame(otu2)) {
      stop("the given object for otu2 is not a data frame.")
    }
    
    tax <- dim(otu2)[2]
    
    # ensure the last column contains taxonomic data
    if (names(otu2)[tax] != "taxonomy") {
      stop("no 'taxonomy' column detected in otu2. Please ensure the last column of the table is titled 'taxonomy' (capitalization matters) and try again.")
    }
    
    # check if all data is numeric other than tax columns
    if (!all(apply(otu2[ ,-tax, drop=FALSE], is.numeric, MARGIN=1))) {
      stop("OTU data other than taxonomic information is not numeric in otu2. Please fix the data and try again.")
    }
    
    # check that all counts are non-negative
    if (any(otu2[ , -tax, drop=FALSE] < 0)) {
      stop("negative counts detected in otu2; cannot continue. Please ensure all counts are positive and try again.")
    }
    
    if (any(colSums(otu2[ ,-tax, drop=FALSE]) == 0)) {
      warning("some samples in otu2 have zero counts. Some functions may behave bizarrely as a result.")
    }
    
    missing.tax <- which(otu2[ , tax, drop=FALSE] == "")
    
    if (length(missing.tax) != 0) {
      missing.tax.percent <- 100 * (length(missing.tax) / dim(otu2)[1]) 
      warning(paste("taxonomic data is missing for ", format(missing.tax.percent, digits=3),
                    "% of OTUs in otu2.", sep=""))
    }
  }
}

.valid.factors <- function(meta, factors, min.factors=1, max.factors=Inf) {
  
  if (max.factors < min.factors || max.factors < 1 || min.factors < 0) {
    stop("something has gone wrong internally; please contact the package maintainer.")
  }

  if (!is.data.frame(meta)) {
    stop("'meta' must be a data frame.")
  }
  
  if (!is.character(factors)) {
    stop("'factors' must be a character vector; see ?RAM.factors for help.")
  }
  
  factors.len <- length(factors)
  
  if (factors.len < min.factors) {
    stop(paste("only", factors.len, "factor(s) were given where", min.factors, "were needed."))
  }
  
  if (length(names(factors)) != factors.len) {
    warning("not all factors are named; converting the names of ALL factors to column names; see ?RAM.factors for help.")
    names(factors) <- factors
  }
  
  if (factors.len > max.factors) {
    warning(paste(factors.len, "factors were given, this function supports up to",
                  max.factors, "factor(s); ignoring the others."))
    factors <- factors[1:max.factors]
    names(factors) <- names(factors)[1:max.factors]
  }
  
  if (!all(factors %in% names(meta))) {
    stop("not all factors were found in 'meta', please check your spelling and try again.")
  }
  
  output <- meta[ , factors, drop=FALSE]
  names(output) <- names(factors)
  
  output
}

.valid.plot.settings <- function(file=NULL, ext=NULL) {
  # ignore NULL values (the default when unused)
  if (is.null(file) && is.null(ext)) {
    return()
  }
  
  # if only one of {ext, file} is specified, raise an error
  if (xor(is.null(ext), is.null(file))) {
    stop("please specify either both 'file' and 'ext', or neither.")
  }
  
  # validate the extension
  .valid.ext(ext)
  
  # validate file path; prepend the device with period (for valid file name)
  file <- .ensure.filepath(file, ext)
  
  file
}

.valid.ext <- function(ext) {
  valid.exts <- c("pdf", "png", "tiff", "svg", "bmp", "jpeg")
  
  if (!ext %in% valid.exts) {
    stop("invalid extension given. See ?RAM.plotting for more details.")
  }
}

.valid.rank <- function(rank) {
  if (!is.character(rank) || !identical(length(rank), 1L)) {
    stop("rank must be a character vector of length 1; see ?RAM.rank.formatting for help.")
  }
  
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax.classes.short <- c("k", "p", "c", "o", "f", "g", "s")
  tax.classes.pattern <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  tax.classes.all <- c(tax.classes, tax.classes.short, tax.classes.pattern)
  
  if (!(rank %in% tax.classes.all)) {
    stop("invalid format for rank. See ?RAM.rank.formatting for help.")
  }
  
  invisible()
}

.valid.meta <- function(otu1, otu2=NULL, meta) {
  valid.OTU(otu1, otu2)
  
  # break when meta == NULL, since meta is an optional argument in some 
  # functions
  if (is.null(meta)) { return(invisible(TRUE)) }
  
  if (!is.data.frame(meta)) {
    stop("'meta' must be a data frame.")
  }
  
  otu1.samples <- names(otu1)[-dim(otu1)[2]]
  meta.samples <- rownames(meta)
  
  if (!identical(otu1.samples, meta.samples)) {
    in.otu1 <- setdiff(otu1.samples, meta.samples)
    in.meta <- setdiff(meta.samples, otu1.samples)
    
    error.msg <- paste("Sample names for otu1 and meta differ. They may be out of order.\nIn otu1, not in meta:\n", 
                       paste(in.otu1, collapse=" "), "\nIn meta, not in otu1:\n", 
                       paste(in.meta, collapse=" "))
    
    stop(error.msg)
  }
  
  if (!is.null(otu2)) {
    otu2.samples <- names(otu2)[-dim(otu2)[2]]
    
    if (!identical(otu2.samples, meta.samples)) {
      in.otu2 <- setdiff(otu2.samples, meta.samples)
      in.meta <- setdiff(meta.samples, otu2.samples)
      
      
      error.msg <- paste("Sample names for otu2 and meta differ. They may be out of order.\nIn otu2, not in meta:\n", 
                         paste(in.otu2, collapse=" "), "\nIn meta, not in otu2:\n", 
                         paste(in.meta, collapse=" "))
      
      stop(error.msg)
    }
  }
  
  invisible(TRUE)
}

.valid.factor <- function(meta, meta.factor){
  if ( !any(meta.factor %in% names(meta)) ) {
     stop("none of the variables are in the metadata")
  } else {
    vec <- vector()
    for ( i in meta.factor) {
      if ( !i %in% names(meta) ) {
        vec <- c(vec, i)
      }
    }
  }
  vec <- unique(vec)
  if ( length(vec) !=0  ) {
    warning(paste("the following variables are not in the metadata: ", paste(vec, collapse=", "), sep=""))
    factors <- meta.factor[-which(meta.factor %in% vec)]
  } else {
    factors <- meta.factor
  }
  return(factors)
}     

# open the appropriate device, given an extension
.get.dev <- function(file, ext, height=8, width=16, dpi=1000) {
  .valid.ext(ext)
  file <- .ensure.filepath(file, ext)
  
  args <- list(file=file, height=height, width=width, res=dpi, bg="white",
               units="in")
  
  if (ext == "tiff") {
    args <- c(args, compression="lzw")
  }
  
  if (ext == "pdf" || ext == "svg") {
    args$res <- NULL
    args$units <- NULL
  }
  
  # this looks strange, but since the functions to open the devices are just called
  # pdf, jpeg, ... themselves, we simply parse our ext value to get the appropriate 
  # function.
  expr <- parse(text=ext)
  
  # evaluate the expression (which gives a function); call with other arguments
  do.call(eval(expr), args)
}

.get.rank.ind <- function(rank) {
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax.classes.plural <- c("kingdoms", "phyla", "classes", "orders", "families", "genera", "species")
  tax.classes.short <- c("k", "p", "c", "o", "f", "g", "s")
  tax.classes.pattern <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  tax.classes.all <- c(tax.classes, tax.classes.plural, tax.classes.short, tax.classes.pattern)
  
  # we do not validate the rank here, as some methods call this 
  # in a for loop when rank=NULL
  
  # convert to upper case as ignore.case in grepl
  # the length call below should return 7 (barring massive taxonomical discoveries...)

   val <- unique(which(toupper(rank)==toupper(tax.classes.all)) %% length(tax.classes))
  # since we took the index mod 7, all species values were given index 0
  # however, we want them to have index 7:
   if (val == 0) {
      val <- 7;
    }
 
  
  return(val)
}


##################################################################################

tax.split <- function(otu1, otu2=NULL, rank=NULL) {
  valid.OTU(otu1, otu2)
  single.otu <- is.null(otu2)
  single.rank <- !is.null(rank)
  tax.classes <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # if given otu and otu2, call tax.split for both
  if (!single.otu) {
    
    output <- list()
    
    output$otu1 <- tax.split(otu1=otu1, rank=rank)
    output$otu2 <- tax.split(otu1=otu2, rank=rank)
    
    return(output)
  }
  
  if (single.rank) {
    # get the index for rank (also does data validation for rank)
    .valid.rank(rank)
    tax.ind <- .get.rank.ind(rank)
  } 
  
  # split OTU table taxonomy information into individual columns
  suppressWarnings(otu1.split <- col.splitup(otu1, col="taxonomy", sep="; ", 
                   max=length(tax.classes), names=tax.classes, drop=TRUE))
  
  # columns from 1 to max contain the numeric data, the other taxonomic
  max <- dim(otu1)[2] - 1
  
  # we need seven taxonomy columns; if one is missing (because nothing classified
  # at that level), add empty column
  
  # while we have less columns than we should...
  #while (dim(otu1.split)[2] < max + length(tax.classes)) {
    # ...add a column filled with empty strings
  #  otu1.split <- cbind(otu1.split, rep("", times=dim(otu1.split)[1]))
  #}
  
  # strip the 'formatting' bits of taxonomic info
  otu1.split[ ,-(1:max)] <- gsub("k__|p__|c__|o__|f__|g__|s__|;", "", 
        as.matrix(otu1.split[ ,-(1:max)]))
  
  if (single.rank) {
    # add the single taxonomic column to the numeric data, return that data frame
    return(otu1.split[ , names(otu1.split) %in% c(names(otu1)[1:max], tax.classes[tax.ind])])
    
  } else {
    # set up list for output
    otu1.taxa <- vector("list", length(tax.classes))
    names(otu1.taxa) <- tax.classes 
    
    for (i in 1:length(tax.classes)) {
      # creating a list of data frames is surprisingly difficult in R;
      # if you do not use the list(df1, df2, ...) syntax, you get a list composed
      # of the first data frame in its entirety, followed be the individual columns
      # of the other data frames. Instead we wrap each data frame in a list itself 
      # before we add it, then we can access them with [[]]
      otu1.taxa[i] <- list(otu1.split[ , names(otu1.split) %in% c(names(otu1)[1:max],
                   tax.classes[i])])
    }
    
    return(otu1.taxa)
  }
}

tax.abund <- function(otu1, otu2=NULL, rank=NULL,            
                      drop.unclassified=FALSE, 
                      top=NULL, count=TRUE, mode="number") {
  
  single.otu <- is.null(otu2)
  valid.OTU(otu1, otu2)
  single.rank <- !is.null(rank)
  # data validation for top is done later in the function (when the dimensions of 
  # the taxonomy matrix are known)
  filter <- !is.null(top)
  
  if (!is.logical(drop.unclassified) || length(drop.unclassified) != 1L) {
    stop("argument 'drop.unclassified' must be a logical vector of length one.")
  }
  
  if (!is.logical(count) || length(count) != 1L) {
    stop("argument 'count' must be a logical vector of length one.")
  }
  
  if (!is.character(mode) || !any(mode %in% c("number", "percent"))) {
    stop("argument 'mode' must be one of 'number' or 'percent'.")
  }
  
  # if given otu and otu2, call tax.abund for both
  if (!single.otu) {
    
    drop <- drop.unclassified
    output <- list()
    
    output$otu1 <- tax.abund(otu1, rank=rank, drop.unclassified=drop, top=top,
         count=count, mode=mode)
    output$otu2 <- tax.abund(otu2, rank=rank, drop.unclassified=drop, top=top,
         count=count, mode=mode)
    
    return(output)
  }
  
  # get the OTU table in proper format
  if (single.rank) {
    .valid.rank(rank)
    tax.list <- list(tax.split(otu1, rank=rank))
  } else {
    tax.list <- tax.split(otu1)
  }
  
  for (i in seq(along=tax.list)) { 
    
    # update taxonomy label to "taxonomy"
    names(tax.list[[i]])[dim(tax.list[[i]])[2]] <- "taxonomy" 
    # aggregate the otu table by taxonomy rank names 
    tax.list[[i]] = stats::aggregate(tax.list[[i]][ , -dim(tax.list[[i]])[2]],  
          by = list(tax.list[[i]]$taxonomy), FUN = .sumcol) 
    # change row names to header provided by aggregate
    rownames(tax.list[[i]]) <- tax.list[[i]][ , 1] 
    # remove first column (that information is now in the row names)
    tax.list[[i]] <- tax.list[[i]][ , -1]
    # transpose table (the as.data.frame generates the V1 heading) 
    tax.list[[i]] <- as.data.frame(t(tax.list[[i]])) 
    
    # if count is false, return relative abundance
    if (!count) {
      tax.list[[i]] <- vegan::decostand(tax.list[[i]], method="total")
    }
    
    # order the table by column sums
    tax.list[[i]] <- tax.list[[i]][ , order(colSums(tax.list[[i]]), 
            decreasing = TRUE), drop=FALSE] 
    
    # remove all zero entries
    tax.list[[i]] <- tax.list[[i]][ , colSums(tax.list[[i]]) > 0, drop=FALSE]
    
    # rename the "V1" column
    if (!single.rank) {
      names(tax.list[[i]])[names(tax.list[[i]]) == "V1"] <- 
    paste("unclassified_below_", .get.rank(i - 1), sep="")
    } else {
      # if we are only processing one element, we cannot use the i index 
      # to get the correct rank
      names(tax.list[[i]])[names(tax.list[[i]]) == "V1"] <- 
    paste("unclassified_below_", .get.rank(.get.rank.ind(rank) - 1), sep="")
    }
    
    # remove unclassified columns
    if (drop.unclassified) {
      # this selects all columns NOT containing in the blacklist
      remove.pat <- gsub(.get.rank.pat(.get.rank(i)), "", paste0(.blacklist(.get.rank(i)), "|no_taxonomy"))
      
      tax.list[[i]] <- tax.list[[i]][ , !grepl(remove.pat, names(tax.list[[i]]),
               ignore.case=TRUE), 
             drop=FALSE]
    }
    
    # keep only the 'top' most abundant groups, where top is user-given 
    if (filter) {
      tax.list[[i]] <- .select.top(tax.list[[i]], top, count, mode)
    }
  }

  tax.out <- list()
  index <- 1
  for ( i in 1:length(tax.list) ) {
    tax <- tax.list[[i]]
    if ( is.null(tax) ) { break }
    lab <- names(tax.list)[i]
    if (is.null(lab)) {
      lab1 <-index
    } else {
       lab1 <- lab
    }
    names(tax) <- gsub(" ", "_", names(tax))
    tax.out[[lab1]] <- tax
    index <- index + 1
  }
  
  if (single.rank) {
    return(tax.out[[1]])
    
  } else {
    return(tax.out)
  }
}

col.splitup <- function(df, col="", sep="", max=NULL, names=NULL, drop=TRUE) { 
    # validate all inputs
    if ( sep == "" ) {
        stop(paste("\n", "    Error: separator ?  check ?col.splitup  ", "\n", sep=""))
    }
    if ( col == "" )  {
       stop(paste("\n", "    Error: column to split? check ?col.splitup  ", "\n", sep=""))
    }
    if ( !(col %in% names(df)) ) {
        stop(paste("\n", "    Error: column to be split is not in the dataframe", "\n", sep=""))
    }  
    
    # col position in df
    if ( col %in% names(df) ) {
        num.col <- which( names(df) %in% col )
    }
    
    # split the column to list;
    list <- strsplit(df[[col]], sep, fixed=FALSE); 
    vec <- vector();

    # determine max number of split columns to be remained in output
    for (i in 1:length(list) ) {
      vec<-c( vec, length(list[[i]]) );
    }

    def.max <- max(c(max, length(names)))
    maximum <- c( max, max(vec), length(names) )
    max <- max(maximum)
    
    if ( max(vec) > def.max ) {
#        warning(paste("\n", "    ", col, " can be split to: ", max(vec), " coloumns; ", 
#                 "\n", "    required to be split to ", def.max, " columns; ", "\n", 
#                  "    Will KEEP all ", max(vec), " and IGNORE user defined ",  def.max, 
#                 " columns", sep=""))
     } else if ( max(vec) < def.max )  {
            #warning(paste("\n", "    ", col, " can be split to: ", max(vec), 
#                " columns; ", "\n", "    required to be split to: ", max, 
 #               " columns; ", "\n", "    column names provided: ", length(names), "\n", 
 #               "    Will fill empty strings in the additional ", def.max-max(vec), 
 #              " column(s)", sep=""))
      } else {
        #warning(paste("\n", "    ", col, " can be split to: ", max(vec), " coloumns; ", 
#                 "\n", "    required to be split to ", def.max, " columns; ", "\n", 
#                  "    ", col, " will be split to: ", max(vec), " columns", sep=""))
     }
   

    for ( i in 1:length(list) ) {
      # since max is equal to or larger than length(list[[i]]) as defined by section above
      if ( length(list[[i]]) < max ) {
        x=rep( "", max-length(list[[i]]) );
        list[[i]] <- c(list[[i]], x)
      } else {
        list[[i]] <- list[[i]]
      }
    } 
    
    # rbind to form a data frame of split columns    
    new<-as.data.frame( do.call("rbind", list) )
    
    # names of new columns
    if ( is.null(names) ) {
       new.name <- colnames(new)
    } else {
      if ( length(names) == max ) {
        new.name <- names
      } else if ( length(names) < max ) {
      #warning(paste("\n", "    ", col, " being split to: ", max, 
      #          " columns;", "\n", "    column names provided: ", length(names), "; ", "\n",
      #          "    will only change the first ", max, " of split columns", sep=""))
        new.name <- c(names, colnames(new)[(length(names)+1):max])
      # } 
      # since max is the maximum of length(names), pre-defined max and max(vec), so 
      # following condition is not possible
      #  else if ( length(names) > max ) {
      #     warning(paste("\n", "    remained number of split columns from ", col, " is: ", max, 
      #          ";", "\n", " number of names provided: ", length(names), "; ", "\n",
      #         "     will ignore the last ", length(names) -max, " of names", sep=""))
      #  new.name <- names[1:max]
      } else {
        new.name <- colnames(new)
      }
    }
    
    colnames(new) <- new.name    
    # for data.table dt[, setdiff(colnames(dt),col), with=FALSE] 
    if ( ! drop ) {
        #warning(paste("    Keep ", col, " column in output!", sep=""))
        df.new <- cbind(df, new)
    } else {
        #warning(paste("    Drop ", col, " column in output!", sep=""))
        df.new <- cbind(df[, setdiff(colnames(df),col)], new)
    }
    return(df.new)
}



plot.OTU.heatmap <- function(data, is.OTU=TRUE, meta=NULL, rank="g", row.factor=NULL,
                                 top=NULL, count=FALSE, drop.unclassified=FALSE,
                                 dendro="none", file=NULL, ext=NULL, 
                                 width=9, height=8, leg.x=-0.08, leg.y=0) {

  save <- !is.null(file)
  if ( is.OTU && !is.null(rank) ) {
    valid.OTU(data)
    .valid.meta(otu1=data, meta=meta)
    .valid.rank(rank)
    data.tax <- tax.abund(data, rank=rank, drop.unclassified=drop.unclassified,
                         count=count, top=top)  
  } else if ( !is.OTU || is.null(rank) ) {
    rank <- NULL    
    if ( count ) {
      stand.method=NULL
    } else {
      stand.method="total"
    }    
    data.tax <- data.revamp(data=list(df=data), is.OTU=is.OTU, 
                          stand.method=stand.method, 
                          ranks=rank, top=top)[[1]]
  }
  
  data.tax <- data.tax[match(rownames(meta), rownames(data.tax)),]
  if ( !identical(rownames(meta), rownames(data.tax)) ) {
    stop("samples are not the same in metadata and community data")
  }
 
  # we create a NULL col.factor for compatability with the internal
  # .valid.factor and .factor.settings functions
  # Note: a column factor is not valid here, as the factors are for samples
  # which are on the rows (the aggregated taxonomy is on the columns)
  col.factor <- NULL
  given.rfactor <- !is.null(row.factor)
  given.cfactor <- !is.null(col.factor)
  
  if (is.null(row.factor)) {
    rfactor <- NULL
  }
  
  if (given.rfactor) {
    rfactor <- .valid.factors(meta, row.factor, max.factors = 1)
    
    # if we have a row factor, but the length doesn't match the number of samples,
    # stop
    if (nrow(rfactor) != nrow(data.tax)) {
      stop(paste("The given data has", nrow(data.tax) - 1, 
                 "samples, but the row.factor has", nrow(rfactor), 
                 "points. These two values must match to continue."))
    }
  }
  
  if (!is.null(meta) && !all(rownames(data.tax) %in% rownames(meta))) {
    stop("the rownames for 'data' and 'meta' do not match.")
  }
  
  # general approach for this function: we build up a list of the arguments
  # we wish to pass to heatmap.2 and legend based on arguments supplied by 
  # the user and some default settings. Once we have created a list with all
  # the arguments, we pass that to do.call (with heatmap.2 or legend) to 
  # create the plot.
  
  if (save) {
    # does data validation for file
    .get.dev(file, ext, height=height, width=width)
  }
  
  # generate a palette with a large gradient
  cols <- RColorBrewer::brewer.pal(3, "YlOrRd") 
  gradient <- colorRampPalette(cols)(100)
  
  longest.row.label <- max(sapply(rownames(data.tax), FUN=nchar))
  longest.col.label <- max(sapply(colnames(data.tax), FUN=nchar))
  longest.factor.label <- max(sapply(levels(factor(meta[[row.factor]])), FUN=nchar))
  
  hmap.args <- list(x=as.matrix(data.tax), dendrogram=dendro, trace = "none",
                    key = TRUE, keysize = 1, density.info = c("none"), 
                    col=gradient, xpd=NA, 
                    # we set the margins based on the length of the labels
                    ### (this probably needs more tweaking!)
                    margins=c(0.6 * longest.col.label, 0.7 * longest.row.label ))
  
  hmap.args$dendrogram <- dendro
  # order based on metadata, if appropriate
  if (given.rfactor) {
    # ensure that the order of the samples matches the order of the metadata
    # exactly, after grouping all identical samples together
    data.tax <- data.tax[order(rfactor[ ,1]), ]
    hmap.args$x <- as.matrix(data.tax)
    
    if (dendro == "row" || dendro == "both") {
      # warning("if metadata is provided, clustering will not occur and no dendrogram will be shown for the rows.")
      leg.args <- list(x="right", inset=c(leg.x, leg.y), xpd=NA, cex=0.7)
      if ( dendro == "row" ) {
        hmap.args <- c(hmap.args, Rowv=TRUE, Colv=FALSE)
      } else {
        hmap.args <- c(hmap.args, Rowv=TRUE, Colv=TRUE)
      }
      #hmap.args$dendrogram <- "none"
      #leg.args <-list(x=locator(), xpd=NA)
      hmap.args$margins <- c(0.6 * longest.col.label, 0.7 * longest.row.label +longest.factor.label)
    } else {
      hmap.args <- c(hmap.args, Colv=TRUE, Rowv=FALSE)
      #leg.args <-list(x=locator(), xpd=NA)
      leg.args <- list(x="left", inset=c(leg.x, leg.y), xpd=NA, cex=0.7)
      hmap.args$margins <- c(0.6 * longest.col.label, 0.7 * longest.row.label )
    }
  }
  
  #leg.args <-list(x=locator(), xpd=NA)
  # leg.args <- list(x="right", inset=c(leg.x, leg.y), xpd=NA, cex=0.7)
  
  args <- .factor.settings(rfactor, NULL, hmap.args, leg.args)
  hmap.args <- args[[1]]
  leg.args <- args[[2]]
  
  do.call(heatmap.2, hmap.args)
  
  if (given.rfactor) {
    do.call(legend, leg.args)
  }
  
  if (save) {
    dev.off()
  }
  
  invisible()
}

Args <-commandArgs(T)
otu_data <- read.table(Args[1],header=T)
colnames(otu_data)[length(otu_data[1,])] <- "taxonomy"
rownames(otu_data) <- otu_data[,1]
otu_data <- otu_data[,-1]
plot.OTU.heatmap(otu_data,top=10, drop.unclassified=TRUE, dendro="both", file="OTU_heatmap",ext="png", width=4800,height=4800)
