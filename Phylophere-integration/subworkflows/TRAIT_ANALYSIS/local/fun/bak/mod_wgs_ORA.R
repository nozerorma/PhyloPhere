library(WebGestaltR)

# Fix zip function
zip <- function(zipfile, files = ".", flags = "-rq", extras = "") {
  # Find system zip
  zip_path <- Sys.which("zip")
  if (zip_path == "") stop("System 'zip' utility not found.")
  
  # Ensure zipfile is a single quoted string (handles spaces)
  zipfile <- shQuote(zipfile)
  
  # files can be a vector, quote each element
  files <- paste(shQuote(files), collapse = " ")
  
  # Construct command args for system2
  args <- c(flags, zipfile, files)
  
  # Call system zip command, pass extras if any
  res <- system2(zip_path, args = args, stdout = TRUE, stderr = TRUE)
  
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0) {
    stop(paste("zip command failed:", paste(res, collapse = "\n")))
  }
  
  invisible(TRUE)
}


mod.WebGestaltR <- function(enrichMethod = "ORA", organism = "hsapiens", enrichDatabase = NULL, 
          enrichDatabaseFile = NULL, enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL, 
          interestGeneFile = NULL, interestGene = NULL, interestGeneType = NULL, 
          interestGeneNames = NULL, collapseMethod = "mean", referenceGeneFile = NULL, 
          referenceGene = NULL, referenceGeneType = NULL, referenceSet = NULL, 
          minNum = 10, maxNum = 500, sigMethod = "fdr", fdrMethod = "BH", 
          fdrThr = 0.05, topThr = 10, reportNum = 20, perNum = 1000, 
          gseaP = 1, isOutput = TRUE, outputDirectory = getwd(), projectName = NULL, 
          dagColor = "continuous", saveRawGseaResult = FALSE, gseaPlotFormat = c("png", 
                                                                                 "svg"), setCoverNum = 10, networkConstructionMethod = NULL, 
          neighborNum = 10, highlightType = "Seeds", highlightSeedNum = 10, 
          nThreads = 1, cache = NULL, hostName = "https://www.webgestalt.org/", 
          useWeightedSetCover = FALSE, useAffinityPropagation = FALSE, 
          usekMedoid = TRUE, kMedoid_k = 25, listName = NULL, ...) 
{
  extraArgs <- list(...)
  if ("keepGSEAFolder" %in% names(extraArgs) | "keepGseaFolder" %in% 
      names(extraArgs)) {
    warning("Parameter keepGSEAFolder is obsolete.\n")
  }
  if ("is.output" %in% names(extraArgs)) {
    isOutput <- extraArgs$is.output
    warning("Parameter is.output is deprecated and changed to isOutput!\n")
    warning("Column names of the result data frame are modified.")
  }
  if ("methodType" %in% names(extraArgs)) {
    warning("Parameter methodType is obsolete.\n")
  }
  if ("lNum" %in% names(extraArgs)) {
    warning("Parameter lNum is obsolete.\n")
  }
  if ("dNum" %in% names(extraArgs)) {
    warning("Parameter dNum is deprecated and changed to reportNum.\n")
    reportNum <- extraArgs$dNum
  }
  if (!is.null(cache)) {
    cat("Use cache data if available.\n")
  }
  if (!is.null(enrichDatabase)) {
    if (length(enrichDatabase) > 1) {
      enrichDatabase <- unlist(sapply(enrichDatabase, 
                                      function(x) {
                                        return(get_gmt_file(hostName, interestGeneType, 
                                                            x, organism, cache))
                                      }))
    }
    else {
      enrichDatabase <- get_gmt_file(hostName, interestGeneType, 
                                     enrichDatabase, organism, cache)
    }
  }
  errorTest <- parameterErrorMessage(enrichMethod = enrichMethod, 
                                     organism = organism, collapseMethod = collapseMethod, 
                                     minNum = minNum, maxNum = maxNum, fdrMethod = fdrMethod, 
                                     sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr, 
                                     reportNum = reportNum, perNum = perNum, isOutput = isOutput, 
                                     outputDirectory = outputDirectory, dagColor = dagColor, 
                                     hostName = hostName, cache = cache)
  if (!is.null(errorTest)) {
    return(errorTest)
  }
  if (is.null(projectName)) {
    projectName <- as.character(as.integer(Sys.time()))
  }
  projectName <- sanitizeFileName(projectName)
  if (enrichMethod == "ORA") {
    enrichR <- mod.WebGestaltROra(organism = organism, enrichDatabase = enrichDatabase, 
                              enrichDatabaseFile = enrichDatabaseFile, enrichDatabaseType = enrichDatabaseType, 
                              enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile, 
                              interestGeneFile = interestGeneFile, interestGene = interestGene, 
                              interestGeneType = interestGeneType, collapseMethod = collapseMethod, 
                              referenceGeneFile = referenceGeneFile, referenceGene = referenceGene, 
                              referenceGeneType = referenceGeneType, referenceSet = referenceSet, 
                              minNum = minNum, maxNum = maxNum, fdrMethod = fdrMethod, 
                              sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr, 
                              reportNum = reportNum, setCoverNum = setCoverNum, 
                              isOutput = isOutput, outputDirectory = outputDirectory, 
                              projectName = projectName, dagColor = dagColor, 
                              nThreads = nThreads, cache = cache, hostName = hostName, 
                              useWeightedSetCover = useWeightedSetCover, useAffinityPropagation = useAffinityPropagation, 
                              usekMedoid = usekMedoid, kMedoid_k = kMedoid_k, 
                              listName = listName)
  }
  else if (enrichMethod == "GSEA") {
    enrichR <- WebGestaltR::WebGestaltRGsea(organism = organism, enrichDatabase = enrichDatabase, 
                               enrichDatabaseFile = enrichDatabaseFile, enrichDatabaseType = enrichDatabaseType, 
                               enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile, 
                               interestGeneFile = interestGeneFile, interestGene = interestGene, 
                               interestGeneType = interestGeneType, collapseMethod = collapseMethod, 
                               minNum = minNum, maxNum = maxNum, fdrMethod = fdrMethod, 
                               sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr, 
                               reportNum = reportNum, setCoverNum = setCoverNum, 
                               perNum = perNum, p = gseaP, isOutput = isOutput, 
                               outputDirectory = outputDirectory, projectName = projectName, 
                               dagColor = dagColor, saveRawGseaResult = saveRawGseaResult, 
                               plotFormat = gseaPlotFormat, nThreads = nThreads, 
                               cache = cache, hostName = hostName, useWeightedSetCover = useWeightedSetCover, 
                               useAffinityPropagation = useAffinityPropagation, 
                               usekMedoid = usekMedoid, kMedoid_k = kMedoid_k, 
                               listName = listName)
  }
  else if (enrichMethod == "NTA") {
    enrichR <- WebGestaltR::WebGestaltRNta(organism = organism, network = enrichDatabase, 
                              method = networkConstructionMethod, neighborNum = neighborNum, 
                              highlightSeedNum = highlightSeedNum, inputSeed = interestGene, 
                              inputSeedFile = interestGeneFile, interestGeneType = interestGeneType, 
                              sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr, 
                              outputDirectory = outputDirectory, projectName = projectName, 
                              highlightType = highlightType, cache = cache, hostName = hostName, 
                              listName = listName)
  }
  return(enrichR)
}

mod.WebGestaltRORA <- function (organism = "hsapiens", enrichDatabase = NULL, enrichDatabaseFile = NULL, 
          enrichDatabaseType = NULL, enrichDatabaseDescriptionFile = NULL, 
          interestGeneFile = NULL, interestGene = NULL, interestGeneType = NULL, 
          collapseMethod = "mean", referenceGeneFile = NULL, referenceGene = NULL, 
          referenceGeneType = NULL, referenceSet = NULL, minNum = 10, 
          maxNum = 500, fdrMethod = "BH", sigMethod = "fdr", fdrThr = 0.05, 
          topThr = 10, reportNum = 20, setCoverNum = 10, isOutput = TRUE, 
          outputDirectory = getwd(), projectName = NULL, dagColor = "binary", 
          nThreads = 1, cache = NULL, hostName = "https://www.webgestalt.org/", 
          useWeightedSetCover = TRUE, useAffinityPropagation = FALSE, 
          usekMedoid = FALSE, kMedoid_k = 10, listName = NULL) 
{
  enrichMethod <- "ORA"
  projectDir <- file.path(outputDirectory, paste0("Project_", 
                                                  projectName))
  enrichDatabase <- testNull(enrichDatabase)
  enrichDatabaseFile <- testNull(enrichDatabaseFile)
  enrichDatabaseType <- testNull(enrichDatabaseType)
  enrichDatabaseDescriptionFile <- testNull(enrichDatabaseDescriptionFile)
  interestGeneFile <- testNull(interestGeneFile)
  interestGene <- testNull(interestGene)
  interestGeneType <- testNull(interestGeneType)
  referenceGeneFile <- testNull(referenceGeneFile)
  referenceGene <- testNull(referenceGene)
  referenceGeneType <- testNull(referenceGeneType)
  referenceSet <- testNull(referenceSet)
  errorTest <- parameterErrorMessage(enrichMethod = enrichMethod, 
                                     organism = organism, collapseMethod = collapseMethod, 
                                     minNum = minNum, maxNum = maxNum, fdrMethod = fdrMethod, 
                                     sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr, 
                                     reportNum = reportNum, isOutput = isOutput, outputDirectory = outputDirectory, 
                                     dagColor = dagColor, hostName = hostName, cache = cache)
  if (!is.null(enrichDatabase)) {
    if (is.character(enrichDatabase) & length(enrichDatabase) == 
        1) {
      if (enrichDatabase == "all") {
        all_sets <- listGeneSet(organism = organism, 
                                hostName = hostName, cache = cache)
        all_sets <- all_sets[all_sets$idType == "entrezgene", 
        ]
        enrichDatabase <- all_sets$name
        enrichDatabaseType <- all_sets$idType
      }
    }
  }
  if (!is.null(errorTest)) {
    stop(errorTest)
  }
  cat("Loading the functional categories...\n")
  enrichD <- loadGeneSet(organism = organism, enrichDatabase = enrichDatabase, 
                         enrichDatabaseFile = enrichDatabaseFile, enrichDatabaseType = enrichDatabaseType, 
                         enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile, 
                         cache = cache, hostName = hostName)
  geneSet <- enrichD$geneSet
  geneSetDes <- enrichD$geneSetDes
  geneSetDag <- enrichD$geneSetDag
  geneSetNet <- enrichD$geneSetNet
  databaseStandardId <- enrichD$standardId
  rm(enrichD)
  cat("Loading the ID list...\n")
  interestingGeneMap <- loadInterestGene(organism = organism, 
                                         dataType = "list", inputGeneFile = interestGeneFile, 
                                         inputGene = interestGene, geneType = interestGeneType, 
                                         collapseMethod = collapseMethod, cache = cache, hostName = hostName, 
                                         geneSet = geneSet)
  if (organism == "others") {
    interestGeneList <- unique(interestingGeneMap)
  }
  else {
    interestStandardId <- interestingGeneMap$standardId
    interestGeneList <- unique(interestingGeneMap$mapped[[interestStandardId]])
  }
  cat("Loading the reference list...\n")
  referenceGeneList <- loadReferenceGene(organism = organism, 
                                         referenceGeneFile = referenceGeneFile, referenceGene = referenceGene, 
                                         referenceGeneType = referenceGeneType, referenceSet = referenceSet, 
                                         collapseMethod = collapseMethod, hostName = hostName, 
                                         geneSet = geneSet, interestGeneList = interestGeneList, 
                                         cache = cache)
  if (isOutput) {
    dir.create(projectDir)
    if (organism != "others") {
      if (databaseStandardId == "entrezgene") {
        cat("Summarizing the input ID list by GO Slim data...\n")
        goSlimOutput <- file.path(projectDir, paste0("goslim_summary_", 
                                                     projectName))
        re <- goSlimSummary(organism = organism, geneList = interestGeneList, 
                            outputFile = goSlimOutput, outputType = "png", 
                            isOutput = isOutput, cache = cache, hostName = hostName)
      }
      write_tsv(interestingGeneMap$mapped, file.path(projectDir, 
                                                     paste0("interestingID_mappingTable_", projectName, 
                                                            ".txt")))
      write(interestingGeneMap$unmapped, file.path(projectDir, 
                                                   paste0("interestingID_unmappedList_", projectName, 
                                                          ".txt")))
    }
    else {
      write(interestGeneList, file.path(projectDir, paste0("interestList_", 
                                                           projectName, ".txt")))
    }
  }
  cat("Performing the enrichment analysis...\n")
  oraRes <- oraEnrichment(interestGeneList, referenceGeneList, 
                          geneSet, minNum = minNum, maxNum = maxNum, fdrMethod = fdrMethod, 
                          sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr)
  if (is.null(oraRes)) {
    return(NULL)
  }
  enrichedSig <- oraRes$enriched
  insig <- oraRes$background
  clusters <- list()
  geneTables <- list()
  if (!is.null(enrichedSig)) {
    if (!is.null(geneSetDes)) {
      enrichedSig <- enrichedSig %>% left_join(geneSetDes, 
                                               by = "geneSet") %>% select(.data$geneSet, .data$description, 
                                                                          .data$link, .data$size, .data$overlap, .data$expect, 
                                                                          .data$enrichmentRatio, .data$pValue, .data$FDR, 
                                                                          .data$overlapId) %>% arrange(.data$FDR, .data$pValue, 
                                                                                                       desc(.data$size)) %>% mutate(description = ifelse(is.na(.data$description), 
                                                                                                                                                         "", .data$description))
    }
    else {
      enrichedSig <- enrichedSig %>% select(.data$geneSet, 
                                            .data$link, .data$size, .data$overlap, .data$expect, 
                                            .data$enrichmentRatio, .data$pValue, .data$FDR, 
                                            .data$overlapId) %>% arrange(.data$FDR, .data$pValue, 
                                                                         desc(.data$size))
    }
    geneTables <- getGeneTables(organism, enrichedSig, "overlapId", 
                                interestingGeneMap)
    if (organism != "others") {
      enrichedSig$link <- mapply(function(link, geneList) linkModification("ORA", 
                                                                           link, geneList, interestingGeneMap, hostName), 
                                 enrichedSig$link, enrichedSig$overlapId)
    }
    if ("database" %in% colnames(geneSet)) {
      enrichedSig <- enrichedSig %>% left_join(unique(geneSet[, 
                                                              c("geneSet", "database")]), by = "geneSet")
    }
    if (organism != "others" && interestGeneType != interestStandardId) {
      outputEnrichedSig <- mapUserId(enrichedSig, "overlapId", 
                                     interestingGeneMap)
    }
    else {
      outputEnrichedSig <- enrichedSig
    }
    if (isOutput) {
      write_tsv(outputEnrichedSig, file.path(projectDir, 
                                             paste0("enrichment_results_", projectName, ".txt")))
      idsInSet <- sapply(enrichedSig$overlapId, strsplit, 
                         split = ";")
      names(idsInSet) <- enrichedSig$geneSet
      minusLogP <- -log(enrichedSig$pValue)
      minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
      apRes <- NULL
      wscRes <- NULL
      kRes <- NULL
      if (useAffinityPropagation) {
        apRes <- affinityPropagation(idsInSet, minusLogP)
      }
      if (useWeightedSetCover) {
        wscRes <- weightedSetCover(idsInSet, 1/minusLogP, 
                                   setCoverNum, nThreads)
      }
      if (usekMedoid) {
        kRes <- kMedoid(idsInSet, minusLogP, maxK = kMedoid_k)
      }
      if (!is.null(apRes)) {
        writeLines(sapply(apRes$clusters, paste, collapse = "\t"), 
                   file.path(projectDir, paste0("enriched_geneset_ap_clusters_", 
                                                projectName, ".txt")))
      }
      else {
        apRes <- NULL
      }
      clusters$ap <- apRes
      if (!is.null(kRes)) {
        tryCatch({
          writeLines(sapply(kRes$clusters, paste, collapse = "\t"), 
                     file.path(projectDir, paste0("enriched_geneset_kmedoid_clusters_", 
                                                  projectName, ".txt")))
        }, error = function(e) {
          cat("Error in writing kmedoid clusters.\n")
        })
      }
      else {
        kRes <- NULL
      }
      clusters$km <- kRes
      if (!is.null(wscRes$topSets)) {
        writeLines(c(paste0("# Coverage: ", wscRes$coverage), 
                     wscRes$topSets), file.path(projectDir, paste0("enriched_geneset_wsc_topsets_", 
                                                                   projectName, ".txt")))
        clusters$wsc <- list(representatives = wscRes$topSets, 
                             coverage = wscRes$coverage)
      }
      else {
        clusters$wsc <- NULL
      }
    }
  }
  if (isOutput) {
    cat("Generate the final report...\n")
    createReport(hostName = hostName, outputDirectory = outputDirectory, 
                 organism = organism, projectName = projectName, 
                 enrichMethod = enrichMethod, geneSet = geneSet, 
                 geneSetDes = geneSetDes, geneSetDag = geneSetDag, 
                 geneSetNet = geneSetNet, interestingGeneMap = interestingGeneMap, 
                 referenceGeneList = referenceGeneList, enrichedSig = enrichedSig, 
                 background = insig, geneTables = geneTables, clusters = clusters, 
                 enrichDatabase = enrichDatabase, enrichDatabaseFile = enrichDatabaseFile, 
                 enrichDatabaseType = enrichDatabaseType, enrichDatabaseDescriptionFile = enrichDatabaseDescriptionFile, 
                 interestGeneFile = interestGeneFile, interestGene = interestGene, 
                 interestGeneType = interestGeneType, collapseMethod = collapseMethod, 
                 referenceGeneFile = referenceGeneFile, referenceGene = referenceGene, 
                 referenceGeneType = referenceGeneType, referenceSet = referenceSet, 
                 minNum = minNum, maxNum = maxNum, fdrMethod = fdrMethod, 
                 sigMethod = sigMethod, fdrThr = fdrThr, topThr = topThr, 
                 reportNum = reportNum, dagColor = dagColor, listName = listName)
    cwd <- getwd()
    setwd(projectDir)
    zip(paste0("Project_", projectName, ".zip"), paste0("Project_", projectName), flags = "-rq")
    setwd(cwd)
    cat("Results can be found in the ", projectDir, "!\n", 
        sep = "")
  }
  return(outputEnrichedSig)
}