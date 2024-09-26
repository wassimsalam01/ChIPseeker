merge_two_si = function(x1, x2){
  if (length(unique(gsub("^[0-9]+","",c(x1, x2)))) == 1){
    return(paste0(gsub("[^0-9]*$","",x1), "-", x2))
  } else {
    return(paste0(x1, "-", x2))
  }
}

generate_break_lbs = function(breaks) {
  lbs <- c()
  
  # break labels
  break_labels = scales::label_number(scale_cut = scales::cut_si(unit = "b"))(breaks)
  break_labels = gsub(" b$"," bp", break_labels)

  # category labels
  for (i in 2:length(breaks)) {
    if (i == length(breaks)) {
      lbs <- c(lbs, paste0(">", break_labels[i-1]))
    } else {
      lbs <- c(lbs, merge_two_si(break_labels[i-1], break_labels[i]))
    }
  }
  
  return(lbs)
}

generate_colors <- function(palette, n) {
  if (length(palette) == 1){
    brewer_cols = RColorBrewer::brewer.pal(
      name = palette, 
      n = min(RColorBrewer::brewer.pal.info[palette, "maxcolors"], n)
    )
  } else if (length(palette) > 1){
    brewer_cols = palette
  } else {
    brewer_cols = c("#9ecae1", "#3182bd", "#C7A76C", "#86B875", "#39BEB1", "#CD99D8")
  }
  
  color_func <- grDevices::colorRampPalette(rev(brewer_cols))
  cols <- color_func(n)
  
  return(cols)
}

##' plot feature distribution based on the distances to the TSS
##'
##'
##' @title plotDistToTSS.data.frame
##' @param peakDist peak annotation
##' @param distanceColumn column name of the distance from peak to nearest gene
##' @param distanceBreaks default is 'c(0, 1000, 3000, 5000, 10000, 100000)'
##' @param palette palette name for coloring different distances. Run `RColorBrewer::display.brewer.all()` to see all applicable values.
##' @param xlab x label
##' @param ylab y lable
##' @param title figure title
##' @param categoryColumn category column, default is ".id"
##' @return bar plot that summarize distance from peak to
##' TSS of the nearest gene.
##' @importFrom magrittr %<>%
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_fill_brewer
##' @importFrom ggplot2 scale_fill_hue
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 geom_text
##' @examples
##' \dontrun{
##' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##' peakfile <- system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
##' peakAnno <- annotatePeak(peakfile, TxDb=txdb)
##' plotDistToTSS(peakAnno)
##' }
##' @seealso \code{\link{annotatePeak}}
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
plotDistToTSS.data.frame <- function(peakDist,
                                     distanceColumn="distanceToTSS",
                                     distanceBreaks=c(0, 1000, 3000, 5000, 10000, 100000),
                                     palette = NULL,
                                     xlab="",
                                     ylab="Binding sites (%) (5'->3')",
                                     title="Distribution of transcription factor-binding loci relative to TSS",
                                     categoryColumn = ".id") {

    distanceBreaks = sort(distanceBreaks)
    hasZero = sum(distanceBreaks == 0)
    if (!hasZero) distanceBreaks = c(0, distanceBreaks)
    hasInf = sum(is.infinite(distanceBreaks))
    if (!hasInf) distanceBreaks = c(distanceBreaks, Inf)
    lbs = generate_break_lbs(distanceBreaks)
    peakDist$Feature = cut(abs(peakDist[[distanceColumn]]), 
                           breaks = distanceBreaks,
                           labels = lbs,
                           include.lowest = TRUE)

    ## sign containing -1 and 1 for upstream and downstream
    peakDist$sign <- sign(peakDist[,distanceColumn])

    ## count frequencies
    if (categoryColumn == 1) {
      peakDist = peakDist |> 
        summarise(freq = length(Feature), .by = c("Feature", "sign")) |> 
        mutate(freq = freq/sum(freq) * 100)
    } else {
      peakDist = peakDist |> 
        summarise(freq = length(Feature), .by = c(categoryColumn,"Feature","sign")) |> 
        mutate(freq = freq/sum(freq) * 100, .by = categoryColumn)
    }

    if (any(peakDist$sign == 0)) {
        zeroDist <- peakDist[peakDist$sign == 0,]
        zeroDist$freq <- zeroDist$freq/2
        zeroDist$sign <- -1
        peakDist[peakDist$sign == 0,] <- zeroDist
        zeroDist$sign <- 1
        peakDist <- rbind(peakDist, zeroDist)
    }

    if (categoryColumn == 1) {
        peakDist %<>% group_by(Feature, sign) %>%
            summarise(freq = sum(freq))

        totalFreq <- peakDist %>% group_by(sign) %>%
            summarise(total = sum(freq))
    } else {
        peakDist %<>% group_by(.id, Feature, sign) %>%
            summarise(freq = sum(freq))
        totalFreq <- peakDist %>% group_by(.id, sign) %>%
            summarise(total = sum(freq))
    }


    ## preparing ylim and y tick labels
    ds = max(totalFreq$total[totalFreq$sign == 1])
    dslim = ceiling(ds/10) * 10
    us = max(totalFreq$total[totalFreq$sign == -1])
    uslim = ceiling(us/10) * 10
    ybreaks <- seq(-uslim, dslim, by=10)
    ylbs <- abs(ybreaks)
    ylbs[ylbs == 0] <- "TSS"

    peakDist$Feature <- factor(peakDist$Feature, levels=rev(levels(peakDist$Feature)))
    if (categoryColumn == 1) {
        p <- ggplot(peakDist, aes(x=1, fill=Feature))
    } else {
        p <- ggplot(peakDist, aes_string(x=categoryColumn, fill="Feature"))
    }

    p <- p + geom_bar(data=subset(peakDist, sign==1), aes(y=freq), stat="identity") +
        geom_bar(data=subset(peakDist, sign==-1), aes(y=-freq), stat="identity")

    p <- p + geom_hline(yintercept = 0, colour = "black") +
        coord_flip() + theme_bw() +
            scale_y_continuous(breaks=ybreaks,labels=ylbs)

    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

    if (categoryColumn == 1) {
        p <- p + scale_x_continuous(breaks=NULL)
    }

    cols <- generate_colors(palette_name = palette, n = length(lbs))
    p <- p + scale_fill_manual(values=rev(cols), guide=guide_legend(reverse=TRUE))

    return(p)
}
