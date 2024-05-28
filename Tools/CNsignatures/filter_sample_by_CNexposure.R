CN_exposure <- readRDS("~/CerCNsig/Randomization/CN_exposure.rds")

CN_exposure <-column_to_rownames(CN_exposure, var = "sample")
tmp <- Filter(is.numeric, CN_exposure[1:7])

CN_exposure_over_3 <-rowSums(sapply(split.default(tmp, sub('_.*', '', names(tmp))), 
               function(x) rowSums(x) > 0)) %>% 
               as.data.frame() %>% 
               dplyr::rename(., "non-zero"=1) %>% 
               filter(`non-zero` %in% 4:7) %>%
               rownames_to_column(.) %>% 
               dplyr::rename(., "Library"=1)

saveRDS(CN_exposure_over_3, "CN_exposure_over_3.rds")
