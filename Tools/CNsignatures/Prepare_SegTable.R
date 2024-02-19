my.path <-  paste(getwd())
filePaths <- list.files(my.path, "\\.txt$", full.names = TRUE)
segment_tables <- lapply(filePaths, read.delim)

for (i in seq_along(segment_tables))
{ cn <- data.frame(segment_tables[[i]])
  segTable <- c()
  for (c in unique(cn$chromosome))
  { snfilt <- cn[cn$chromosome==c,]
    sn.rle <- rle(snfilt[,"segVal"])
    starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
    ends <- cumsum(sn.rle$lengths)
    lapply(1:length(sn.rle$lengths), function(s) {
      from <- snfilt$start[starts[s]]
      to <- snfilt$end[ends[s]]
      segValue <- sn.rle$values[s]
      c(snfilt$chromosome[starts[s]], from, to, segValue)
    }) -> segtmp
    segTableRaw <- data.frame(matrix(unlist(segtmp), ncol = 4, byrow = T), stringsAsFactors = F)
    colnames(segTableRaw) <- c("chromosome", "start", "end", "segVal")
    segTable <- rbind(segTable, segTableRaw)
  }  
  if (nrow(cn)>0)
    segTable[,1] <- paste0("chr", segTable[,1])
    segment_tables[[i]] <- segTable
}


#Final version
##collapse equal value neighboring segments
segment_tables <-list()
for(i in SampleSheet$Seq.ID)
{
  if(file.exists(paste0(i, "_segments.tsv")))
  {
    cn <- read.table(paste0(i, "_segments.tsv"),header=T,
                   sep="\t",stringsAsFactors = F)
                   
  segTable <- c()
  for (c in unique(cn$chromosome))
  { snfilt <- cn[cn$chromosome==c,]
    sn.rle <- rle(snfilt[,"segVal"])
    starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
    ends <- cumsum(sn.rle$lengths)
    lapply(1:length(sn.rle$lengths), function(s) {
      from <- snfilt$start[starts[s]]
      to <- snfilt$end[ends[s]]
      segValue <- sn.rle$values[s]
      c(snfilt$chromosome[starts[s]], from, to, segValue)
    }) -> segtmp
    segTableRaw <- data.frame(matrix(unlist(segtmp), ncol = 4, byrow = T), stringsAsFactors = F)
    colnames(segTableRaw) <- c("chromosome", "start", "end", "segVal")
    segTable <- rbind(segTable, segTableRaw)
  }  
    if (nrow(cn)>0)
    {	
    segment_tables[[i]] <- segTable
    }
  }
}
