extractCopynumberFeaturesDrews = function(CN_data, cores = 1, allowedError = 0.1, rmNorm = FALSE, build="hg19") {

    #get chromosome lengths
    chrlen<-read.table(paste(this_path,"data/hg19.chrom.sizes.txt",sep="/"),sep="\t",stringsAsFactors = F)[1:24,]
    
    #get centromere locations
    gaps<-read.table(paste(this_path,"data/gap_hg19.txt",sep="/"),sep="\t",header=F,stringsAsFactors = F)
    centromeres<-gaps[gaps[,8]=="centromere",]

     if(cores > 1) {
        require(foreach)
        doMC::registerDoMC(cores)
        
        temp_list = foreach::foreach(i=1:5) %dopar% {
            if(i == 1){
                list(segsize = getSegsizeDrews(CN_data, rmNorm = rmNorm) )
            } else if (i == 2) {
                list(bp10MB = getBPnumDrews(CN_data,chrlen) )
            } else if (i == 3) {
                list(osCN = getOscillationDrews(CN_data,chrlen) )
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCountsDrews(CN_data,centromeres,chrlen) )
            } else  {
                list(changepoint = getChangepointCNDrews(CN_data, allowedError, rmNorm = rmNorm) )
            } 
        }
        doParallel::stopImplicitCluster()

        # Another failsafe that the outcome is definitely numeric
        temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) {
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )

    } else {
        # Single core usage
        segsize<-getSegsizeDrews(CN_data,rmNorm = rmNorm)
        bp10MB<-getBPnumDrews(CN_data,chrlen)
        osCN<-getOscillationDrews(CN_data,chrlen)
        bpchrarm<-getCentromereDistCountsDrews(CN_data,centromeres,chrlen)
        changepoint<-getChangepointCNDrews(CN_data,rmNorm = rmNorm)
        copynumber<-getCNDrews(CN_data)

        temp_list = list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
        #temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) {
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )

    }

}

getSegsizeDrews = function(abs_profiles, rmNorm = FALSE) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)

    # Loop over samples
    for(i in samps) {

        # Get segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Make sure segment values are numeric
        segTab$segVal = as.numeric(segTab$segVal)

        # If wished, don't consider normal segments
        if(rmNorm) { segTab = segTab[ segTab$segVal != 2, ] }

        # Avoiding potential artefact
        segTab$segVal[segTab$segVal<0] = 0
        seglen = segTab$end-segTab$start
        seglen = seglen[seglen>0]

        # Double tap.
        out = rbind(out,cbind(ID=rep(i,length(seglen)),value=as.numeric(seglen)))
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))

}

getBPnumDrews = function(abs_profiles,chrlen, SIZE = 10000000) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)

    # Loop over samples
    for(i in samps) {
        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Loop over chromosomes and identify breaks
        chrs = unique(segTab$chromosome)
        allBPnum = c()
        for(c in chrs) {
            currseg = segTab[segTab$chromosome == c,]
            intervals = seq(1, chrlen[chrlen[,1] == paste0("chr",c),2]+SIZE, SIZE)
            res = graphics::hist(as.numeric(currseg$end[-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
            allBPnum = c(allBPnum,res)
        }
        # Make sure it's really numeric
        out = rbind(out, cbind(ID = rep(i,length(allBPnum)),
                               value = as.numeric(allBPnum)))
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}

getOscillationDrews = function(abs_profiles, chrlen) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)
    # Loop over samples
    for(i in samps) {

        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Loop over chromosomes to identify oscillation
        chrs = unique(segTab$chromosome)
        oscCounts = c()
        for(c in chrs) {

            currseg = as.numeric(segTab$segVal[segTab$chromosome == c])
            currseg = round(as.numeric(currseg))

            # Only take chains into consideration with a length of more than 3 elements
            if(length(currseg)>3) {
                prevval = currseg[1]
                count = 0
                for(j in 3:length(currseg)) {
                    if(currseg[j] == prevval & currseg[j] != currseg[j-1]) {
                        count = count+1
                    } else {
                        oscCounts = c(oscCounts,count)
                        count = 0
                    }
                    prevval = currseg[j-1]
                }
            }
        }
        # Make sure it's really numeric
        out = rbind(out, cbind(ID = rep(i,length(oscCounts)),
                               value = as.numeric(oscCounts)))
        if(length(oscCounts) == 0) {
            out = rbind(out,cbind(ID = i, value = 0))
        }
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}

getCentromereDistCountsDrews = function(abs_profiles,centromeres,chrlen) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)
    # Loop over samples
    for(i in samps) {

        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Loop over chromosomes
        chrs = unique(segTab$chromosome)
        all_dists = c()
        for(c in chrs) {
            if(nrow(segTab) > 1) {
                starts = as.numeric(segTab$start[segTab$chromosome==c])[-1]
                segstart = as.numeric(segTab$start[segTab$chromosome==c])[1]
                ends = as.numeric(segTab$end[segTab$chromosome==c])
                segend = ends[length(ends)]
                ends = ends[-length(ends)]
                centstart = as.numeric(centromeres[substr(centromeres[,2],4,5)==c,3])
                centend = as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
                chrend = chrlen[substr(chrlen[,1],4,5)==c,2]
                ndist = cbind(rep(NA,length(starts)),rep(NA,length(starts)))
                ndist[starts<=centstart,1] = (centstart-starts[starts<=centstart])/(centstart-segstart)*-1
                ndist[starts>=centend,1] = (starts[starts>=centend]-centend)/(segend-centend)
                ndist[ends<=centstart,2] = (centstart-ends[ends<=centstart])/(centstart-segstart)*-1
                ndist[ends>=centend,2] = (ends[ends>=centend]-centend)/(segend-centend)
		ndist = stats::na.omit(ndist)
		ndist = apply(ndist,1,min)

                all_dists = rbind(all_dists,sum(ndist>0))
                all_dists = rbind(all_dists,sum(ndist<=0))
            }
        }
        if(nrow(all_dists)>0) {
            # Make sure it's really numeric
            out = rbind(out,cbind(ID=i,ct1=as.numeric(all_dists[,1])))
        }
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}

getChangepointCNDrews = function(abs_profiles, allowedError = 0.1, rmNorm = FALSE) {

    # Prepare looping
    out = c()
    samps = names(abs_profiles)
    # Loop over samples
    for(i in samps) {

        # Retrieve segments
        segTab = abs_profiles[[i]]
        colnames(segTab)[4] = "segVal"

        # Initiate and prepare looping over chromosomes
        segTab$segVal = as.numeric(segTab$segVal)
        segTab$segVal[segTab$segVal<0] = 0
        chrs = unique(segTab$chromosome)
        allcp = c()

        # Loop over chromosomes
        for(c in chrs) {
            currseg = as.numeric(segTab$segVal[segTab$chromosome==c])
            firstSeg = abs(2 - currseg[1] )
            # As we look only at the left end of a CNA, we might miss a changepoint at the beginning of the p-arm
            # That's why we check manually but only regard this value if it is higher than an allowed error rate.
            if(firstSeg <= allowedError) {
                theseChanges = abs(currseg[-1]-currseg[-length(currseg)])
                if(rmNorm) { theseChanges = theseChanges[ currseg[-1] != 2 ] }
                allcp = c(allcp, theseChanges)
            } else {
                theseChanges = c( firstSeg, abs(currseg[-1]-currseg[-length(currseg)]) )
                if(rmNorm) { theseChanges = theseChanges[ currseg != 2 ] }
                allcp = c(allcp, theseChanges)
            }

        }
        if(length(allcp)==0) {
            allcp = 0 #if there are no changepoints
        }
        # Make sure it's really numeric
        out = rbind(out,cbind(ID=rep(i,length(allcp)),value=as.numeric(allcp)))
    }

    # Prepare return
    rownames(out) = NULL
    return(data.frame(out,stringsAsFactors = FALSE))
}