

#' We use exclusive markers to define partitionning population centers to which we associate the closests cells (no combinatorial genes)

doCellPartition <- function(
                          known_template_file=paste0(dataset_path, "/markers3.csv"),
                          readcounts=m$getReadcounts(data_status='Raw'),
                          cell_level_min_step1 = 2,
                          cell_level_min_step2 = 1,
                          normalize=TRUE
                          ){

  cat(paste0("Load population definitions (", known_template_file,') \n'))

  # Load population definitions
  literature_markers.df = read.table(file=known_template_file, sep='\t', stringsAsFactors=F, header=TRUE, skipNul=T)

  # Make list of genes
  literature_markers.df <- literature_markers.df %>%
                    dplyr::mutate(
                        gmap_step1=strsplit(gsub(" ", "", Genes_map_step1, fixed = TRUE), split=","),
                        gmap_step2=strsplit(gsub(" ", "", Genes_map_step2, fixed = TRUE), split=","),
                        Neural_pop_unique=make.unique(Neural_pop), # New for MN
                        name=case_when(Type %in% c('Progenitor', 'Neuron') ~ Neural_pop, !Type %in% c('Progenitor', 'Neuron') ~ Type), # No need for rownames (?)
                        name_unique=make.unique(name) # No need for rownames (?)
                        )

  # ' `r literature_markers.df %>% dplyr::select(-c(gmap_step1, gmap_step2, DV, name, name_unique, Neural_pop_unique)) %>% knitr::kable(format="html", row.names=FALSE, caption="Published population markers used for data partition", escape=F) %>% kableExtra::kable_styling(bootstrap_options = "striped", font_size=12, full_width=F)`

  #' ### Step 1 / Progenitor vs Neuron vs extra populations

  cat(paste0("Step 1: Split macro populations (", paste0(unique(literature_markers.df$Type), collapse=", "), ") \n"))

  #' Create matrix of target population "centers"
  pop_def_mask_step1 = matrix(0,
                              nrow=length(unique(unlist(literature_markers.df$gmap_step1))), 
                              ncol=length(unique(literature_markers.df$Type)),
                              dimnames=list(
                                unique(unlist(literature_markers.df$gmap_step1)),
                                unique(literature_markers.df$Type)
                                )
                              )

  for(gl in unique(literature_markers.df$Type)) {
    pop_def_mask_step1[
          literature_markers.df %>% dplyr::filter(Type==gl) %>% dplyr::pull(gmap_step1) %>% unlist
          , gl] <- 1
  }

  #' Filter lower UMI count for distance calculation (less than 2)

  #' Associate cells to the closest pre-defined center (Euclidean distance)
  cell_levels_step1 = 1*(readcounts[unique(unlist(literature_markers.df$gmap_step1)),] >= cell_level_min_step1)
 
  if(normalize){
    # cell_levels_step1 <- apply(cell_levels_step1, 2, function(x) x/sum(x))
    cell_levels_step1 <- apply(cell_levels_step1, 2, function(x) x/sqrt(sum(x^2)))
    cell_levels_step1[is.na(cell_levels_step1)] <- 0
   
    # pop_def_mask_step1 <- apply(pop_def_mask_step1, 2, function(x) x/sum(x))
    pop_def_mask_step1 <- apply(pop_def_mask_step1, 2, function(x) x/sqrt(sum(x^2)))
    pop_def_mask_step1[is.na(pop_def_mask_step1)] <- 0
  }

  cell_centers_alldist_step1 = as.matrix(pdist::pdist(t(cell_levels_step1), t(pop_def_mask_step1))) # euclidean distance, no other option in pdist
  cell_center_id_step1 = apply(cell_centers_alldist_step1, 1, which.min)

  #' Store associated type
  Type_step1 = factor(colnames(pop_def_mask_step1)[cell_center_id_step1], levels=unique(literature_markers.df$Type))

  #' ### Step 2 / Independent partition of progenitors and neurons

  cat("Step 2: Split neural populations \n")

  #' Create matrices of target progenitor and neuron population centers

  Type_step2_unique = rep("dummy", ncol(readcounts))

  pop_def_mask_step2 = sapply(c('Progenitor', 'Neuron'), function(pop){
      
      pop_markers = literature_markers.df %>%
              dplyr::filter(Type==pop) %>%
              {'names<-'(.$gmap_step2, .$Neural_pop_unique)} %>%
              {lapply(., function(x) x[x %in% rownames(readcounts)])}

      mask = matrix(0,
                    nrow=length(unique(unlist(pop_markers))), 
                    ncol=length(pop_markers),
                    # ncol=length(pop_markers)+1,
                    dimnames=list(
                      unique(unlist(pop_markers)),
                      names(pop_markers)
                      # c(names(pop_markers), paste0(pop, '_NULL'))
                      )
                    )
      for(gl in names(pop_markers)) {
        mask[pop_markers[[gl]], gl] <- 1
      }
      mask
    })

  #' Associate cells to the closest pre-defined center (Euclidean distance)
  for(pop in c('Progenitor', 'Neuron')){
    pop_ids = which(Type_step1 == pop)
    cell_levels_step2 = 1*(readcounts[rownames(pop_def_mask_step2[[pop]]), pop_ids] >= cell_level_min_step2)
    
    if(normalize){
      # cell_levels_step2 <- apply(cell_levels_step2, 2, function(x) x/sum(x))
      cell_levels_step2 <- apply(cell_levels_step2, 2, function(x) x/sqrt(sum(x^2)))
      cell_levels_step2[is.na(cell_levels_step2)] <- 0

      # pop_def_mask_step2[[pop]] <- apply(pop_def_mask_step2[[pop]], 2, function(x) x/sum(x))
      pop_def_mask_step2[[pop]] <- apply(pop_def_mask_step2[[pop]], 2, function(x) x/sqrt(sum(x^2)))
      pop_def_mask_step2[[pop]][is.na(pop_def_mask_step2[[pop]])] <- 0
    }

    cell_centers_alldist_step2 = as.matrix(pdist::pdist(t(cell_levels_step2), t(pop_def_mask_step2[[pop]]))) # euclidean distance, no other option in pdist
    cell_center_id_step2 = apply(cell_centers_alldist_step2, 1, which.min)
    Type_step2_unique[pop_ids] <- colnames(pop_def_mask_step2[[pop]])[cell_center_id_step2]
  }

  # remaining "dummy" types are non-neural tissues
  Type_step2_unique[which(Type_step2_unique=="dummy")] <- as.character(Type_step1[which(Type_step2_unique=="dummy")])

  # Replace temparary neural pop names (for MN)
  popname_map = literature_markers.df %>% dplyr::filter(Type %in% c('Progenitor', 'Neuron')) %>% {"names<-"(.$Neural_pop, .$Neural_pop_unique)}
  Type_step2 <- popname_map[Type_step2_unique]

  Type_step2 <- factor(Type_step2, levels=unique(literature_markers.df$name))
  Type_step2_unique <- factor(Type_step2_unique, levels=unique(literature_markers.df$name_unique))

  # DV = (literature_markers.df %>% tibble::rownames_to_column("Pop") %>% {'names<-'(.$DV, .$Pop)})[as.character(Type_step2)]
  DV = (literature_markers.df %>% {'names<-'(.$DV, .$name)})[as.character(Type_step2)]
  DV[is.na(DV)] <- max(DV, na.rm=TRUE) + 1

  bothSteps_markers = setNames(
                          apply(literature_markers.df, 1, function(x) c(x[['gmap_step1']], x[["gmap_step2"]])) %>% {lapply(., function(x) x[x %in% rownames(readcounts)])},
                          literature_markers.df$name_unique)

  bothSteps_markers_neural_unique = bothSteps_markers[literature_markers.df %>% dplyr::filter(Type %in% c('Neuron', 'Progenitor')) %>% .$name_unique]

  bothSteps_markers_neural = literature_markers.df %>%
                                dplyr::mutate(allmarkers=bothSteps_markers) %>%           # reuse merged gene list
                                dplyr::filter(Type %in% c('Progenitor', 'Neuron')) %>%    # filter non neural
                                dplyr::mutate(Neural_pop=factor(Neural_pop, levels=unique(Neural_pop))) %>% # factor to keep order
                                dplyr::group_by(Neural_pop) %>%
                                dplyr::summarise(gl=list(unique(unlist(allmarkers))), step2_length=length(unique(unlist(gmap_step2)))) %>%
                                dplyr::filter(step2_length!=0) %>% # remove null neural populations
                                {"names<-"(.$gl, .$Neural_pop)}

  return(list(
          Type_step1=Type_step1,
          Type_step2=Type_step2,
          Type_step2_unique=Type_step2_unique,
          DV=DV,
          step1_markers=literature_markers.df %>% dplyr::select(Type, gmap_step1) %>% dplyr::distinct(Type, .keep_all=T) %>% {"names<-"(.$gmap_step1, .$Type)},
          step2_markers=literature_markers.df %>% dplyr::filter(Type %in% c('Neuron', 'Progenitor')) %>% {"names<-"(.$gmap_step2, .$name_unique)},
          bothSteps_markers=bothSteps_markers,
          bothSteps_markers_neural=bothSteps_markers_neural,
          bothSteps_markers_neural_unique=bothSteps_markers_neural_unique
    ))

}

getPopulationColors <- function(known_template_file=paste0(dataset_path, "/markers3.csv")){
  pop_colors = list()

  literature_markers.df = read.table(file=known_template_file, sep='\t', stringsAsFactors=F, header=TRUE, skipNul=T)

  pop_colors$Step1 = literature_markers.df %>% dplyr::select(Type, Step1_Color) %>% dplyr::distinct() %>% {"names<-"(.$Step1_Color, .$Type)}

  pop_colors$Step2 = literature_markers.df %>% dplyr::filter(Type %in% c('Neuron', 'Progenitor')) %>% dplyr::select(Neural_pop, Step2_Color) %>% dplyr::distinct() %>% {"names<-"(.$Step2_Color, .$Neural_pop)}

  pop_colors$DV = literature_markers.df %>% dplyr::filter(Type=='Progenitor' & !is.na(DV)) %>% dplyr::select(DV, DV_Color) %>% dplyr::distinct() %>% {"names<-"(.$DV_Color, .$DV)}

  pop_colors$timepoint = setNames(colorRampPalette(c("pink", "red"))(n=5), seq(9.5, 13.5, 1))

  return(pop_colors)
}

markersToMask <- function(markers){
  mask = matrix(0,
                nrow=length(unique(unlist(markers))), 
                ncol=length(markers),
                dimnames=list(
                  unique(unlist(markers)), 
                  names(markers)
                  )
                )
  for(pn in names(markers)) {
    mask[markers[[pn]], pn] <- 1
  }
  return(mask)
}

estimateDoublets <- function(){

  type_step1_map = list(
          "Neural"=c("Progenitor", "Neuron"),
          "Blood"=c("Blood", "Hematopoeitic", "Erythropoeitic", "Erythrocytes", "Erythrocytes II"),
          "Mesoderm"=c("Mesoderm I", "Mesoderm II", "Mesoderm III", "Mesoderm IV", "Mesoderm V", "Mesoderm VI"),
          "Skin"=c("Skin")
          )
  type_step1_genes = lapply(type_step1_map, function(st){
     unique(unlist(cell_partition$step1_markers[st]))
    })
  type_step1_genes <- Filter(function(x) length(x)>0, type_step1_genes)

  dcounts = matrix(0, nrow=length(type_step1_map), ncol=length(type_step1_genes), dimnames=list(names(type_step1_map), names(type_step1_genes)))

  for(st in names(type_step1_map)){
    print(type_step1_map[[st]])
    st_cells = which(pData(m$expressionSet)$Type_step1 %in% type_step1_map[[st]])
    print(length(st_cells))
    for(st2 in names(type_step1_genes)){
      subset = m$getReadcounts('Raw')[type_step1_genes[[st2]], st_cells, drop=F]
      dcounts[st, st2] = sum(subset>0) / (ncol(subset)*nrow(subset))
    }
  }
  diag(dcounts) <- NA

  return(dcounts)
}


resampleDataForCombinationDE <- function(
              orig_m_neural_path = paste0(m$plot_folder, '/m_neural.rds'),
              num_min_cells_per_type2 = 40,
              selected_type1 = c("Neuron", "Progenitor"), # "Neuron", "Progenitors"
              sampling_categories = c('Type_step2', 'timepoint'),
              num_cells_per_category = 20, 
              dispersion_prefiltering_zscore=0
              ){

  set.seed(1)

  m_combDE = readRDS(orig_m_neural_path)

  # Select either Progenitors or Neurons
  excluded_cell_ids1 = unique(c(
                        which(!pData(m_combDE$expressionSet)$Type_step1 %in% selected_type1)
                        ))
  m_combDE$excludeCellsFromIds(excluded_cell_ids1)

  # Store all Type_step2
  initial_type2 = levels(droplevels(pData(m_combDE$expressionSet)$Type_step2))

  # Eliminate small population types
  excluded_type2 = pData(m_combDE$expressionSet) %>% dplyr::count(Type_step2) %>% dplyr::filter(n < num_min_cells_per_type2) %>% .$Type_step2
  excluded_cell_ids2 = unique(c(
                        which(pData(m_combDE$expressionSet)$Type_step2 %in% excluded_type2)
                        ))

  m_combDE$excludeCellsFromIds(excluded_cell_ids2)

  m_combDE$excludeUnexpressedGenes(min.cells=1, min.level=1, verbose=T, data_status='Raw')

  # we keep the initial Type_step2 levels, even if empty
  pData(m_combDE$expressionSet)$Type_step2 <- factor(
                                                  pData(m_combDE$expressionSet)$Type_step2,
                                                  levels=initial_type2
                                                  )

  # Resample
  categories = expand.grid(
        lapply(sampling_categories, function(x){as.character(sort(unique(pData(m_combDE$expressionSet)[,x])))}),
        stringsAsFactors=F
        )
  colnames(categories) <- sampling_categories

  cell_ids_resampled = lapply(seq(nrow(categories)), function(i){
        
        x = categories[i, , drop=F]
        # print(x)

        pop_ids = Reduce(intersect, lapply(sampling_categories, function(cat){which(pData(m_combDE$expressionSet)[, cat] == x[[cat]])}))

        if(length(pop_ids)==0){
          print(paste0(paste0(x, collapse='/'), " is empty."))
          return(c())
        } else {
          sample(pop_ids, num_cells_per_category, replace=T)
        }
    })
  names(cell_ids_resampled) <- unlist(apply(categories, 1, paste0, collapse='_'))

  cell_ids_resampled <- cell_ids_resampled[lapply(cell_ids_resampled, length)>0] 

  # Combine subsamples into a single expressionSet
  # very slow with many subsets...
  t0=Sys.time()
  m_combDE$expressionSet = Reduce(BiocGenerics::combine,
                lapply(names(cell_ids_resampled), function(popn){

                    resample_eset = m_combDE$expressionSet[, cell_ids_resampled[[popn]]]
                    sampleNames(resample_eset) <- paste0(popn, "_", seq(length(cell_ids_resampled[[popn]])))
                    resample_eset
                  }))
  Sys.time()-t0

  m_combDE$readcounts_raw <- exprs(m_combDE$expressionSet)

  # Some genes may have disappeard with subsampling
  m_combDE$excludeUnexpressedGenes(min.cells=1, min.level=1, verbose=T, data_status='Raw')

  if(!is.null(dispersion_prefiltering_zscore)) {
    # Pre-filter genes with dispersion
    m_combDE$identifyGeneModules(
                  method="Dispersion_DR",
                  data_status='Raw',
                  zscore_threshold=dispersion_prefiltering_zscore,
                  mean_threshold=-Inf,
                  clustering_mode=F)

    length(unlist(m_combDE$dispersion_DR$genemodules))

    m_combDE$excludeGenesFromIds(which(!m_combDE$getGeneNames() %in% unlist(m_combDE$dispersion_DR$genemodules)))
  }

  return(m_combDE)
}


runSingleMonocleDE <- function(
  m_combDE,
  fullModelFormulaStr="~Type_step2", # "~type+timepoint"
  reducedModelFormulaStr = "~1"
  ){

  HSMM0 <- monocle::newCellDataSet(
              as.matrix(m_combDE$getReadcounts(data_status='Raw')), #+.1,
              phenoData = new("AnnotatedDataFrame", data = pData(m_combDE$expressionSet)),
              featureData = new('AnnotatedDataFrame', data = data.frame(row.names=m_combDE$getGeneNames(), "gene_short_name"=m_combDE$getGeneNames())),
              expressionFamily=VGAM::negbinomial.size() # appropriate for UMI count #http://cole-trapnell-lab.github.io/monocle-release/docs/
              )

  HSMM0 <- estimateSizeFactors(HSMM0) # only with VGAM::negbinomial.size()
  HSMM0 <- estimateDispersions(HSMM0) # only with VGAM::negbinomial.size()

  test0 = monocle::differentialGeneTest(HSMM0,
                          fullModelFormulaStr=fullModelFormulaStr,
                          reducedModelFormulaStr=reducedModelFormulaStr,
                          cores=1,
                          relative_expr=F, 
                          verbose=T)

  return(test0)
}

generateCombinatorialModels <- function(
                                  m_combDE,
                                  mode="exponential" # "exponential", "top_n_uplet", 
                                  ){

  if(mode=="exponential"){
  
    # gene_domain_comb_tbl = 
    n <- length(levels(pData(m_combDE$expressionSet)$Type_step2))
    l <- rep(list(0:1), n)
    allcomb = expand.grid(l)
    colnames(allcomb) <- levels(pData(m_combDE$expressionSet)$Type_step2)

    gene_domain_comb_tbl = as_tibble(allcomb) %>% 
                tidyr::unite(model, remove=F, sep="") %>% # add model name (concatenate all design)
                dplyr::mutate(model=paste0("m", model)) %>%
                dplyr::filter(!model %in% c(
                    paste0("m", paste0(rep(1, n), collapse="")),
                    paste0("m", paste0(rep(0, n), collapse=""))
                )) %>% # remove the no/all domain(s) models
                tidyr::crossing(tibble(genename=m_combDE$getGeneNames()), ., by=NULL)

  } else if (mode=="top_n_uplet"){

    # identify domains 

    levels = m_combDE$getReadcounts('Raw') #[selgenes_test, selcells_test]
    bin_levels = silent.binarize.array(levels)
    levels_tbl = as_tibble(as.data.frame(levels)) %>%
                      tibble::rownames_to_column("genename") %>%
                      tidyr::gather(cellname, level, -genename)
    bin_levels_tbl = as_tibble(as.data.frame(bin_levels)) %>%
                      tibble::rownames_to_column("genename") %>%
                      tidyr::gather(cellname, bin_level, -genename)

    # new strategy, we use a score mixing pos ratio and mean level, and do all decreasing n-uplet
    gene_domain_tbl = bin_levels_tbl %>%
              bind_cols(levels_tbl %>% dplyr::select(level)) %>%
              dplyr::left_join(pData(m_combDE$expressionSet) %>%
                                  tibble::rownames_to_column("cellname") %>%
                                  dplyr::select(cellname, Type_step2),
                               by="cellname") %>%
              # dplyr::filter(genename=="Adsl", Type_step2=="dp5") # test gene
              dplyr::group_by(genename, Type_step2) %>%
              dplyr::summarise(
                    domain_ratio_pos=mean(bin_level) + 1e-10, # 1e-10 -> we want to have a non-zero score for 0 bin domains
                    domain_ratio_expressed=mean(1*(level>0)),
                    domain_mean_level=mean(level),
                    domain_score=(domain_ratio_pos*domain_mean_level)/(domain_ratio_pos+domain_mean_level)
                    ) 

    gene_domain_score_wide =  gene_domain_tbl %>% 
                dplyr::filter(domain_score > 0)  %>%  # filter domains with null score
                dplyr::arrange(genename, desc(domain_score)) %>%
                # dplyr::top_n(4, domain_avg) %>% # select maximum number of domain considered per gene 
                dplyr::select(-domain_ratio_pos, -domain_mean_level, -domain_ratio_expressed) %>%
                tidyr::spread(Type_step2, domain_score, fill=0, drop=FALSE) %>%
                dplyr::ungroup()

    # For a gene having N positive domains (eg N=3, p3 p0 dp4 with avg level 1.5 3 2)
    # we create N partitions with 1 singleton, 1 pair, 1 triplet, ...etc
    # ordered by domain score (eg p3 p0 dp4, p0 dp4 and dp4)
    createNuplets <- function(x){
      values = data.frame(x[2:length(x)])
      valuesBin = values
      valuesBin[values>0] <- 1
      numpos = length(which(values > 0))
      proc_rows = rbind(
        valuesBin,
        do.call(rbind, lapply(seq(numpos-1), function(i){
            values2 = values
            values2[, order(as.matrix(values), decreasing=T)[numpos:(numpos-i+1)]] <- 0
            values2[values2>0] <- 1
            values2
        }))
        )
      res=cbind.data.frame(genename=rep(unlist(unname(x[1])), numpos), proc_rows)
      res
    }

    gene_domain_comb_tbl = gene_domain_score_wide %>%
                    dplyr::rowwise() %>%
                    dplyr::do(createNuplets(.)) %>%
                    dplyr::ungroup() %>%
                    tidyr::unite(model, -genename, remove=F, sep="") %>% # add model name (concatenate all design)
                    dplyr::mutate(model=paste0("m", model)) %>%
                    dplyr::filter(!model %in% c(
                        paste0("m", paste0(rep(1, length(levels(pData(m_combDE$expressionSet)$Type_step2))), collapse="")),
                        paste0("m", paste0(rep(0, length(levels(pData(m_combDE$expressionSet)$Type_step2))), collapse=""))
                    )) # remove the no/all domain(s) models

  }

  # remove models containing empty types
  empty_types = names(which(table(pData(m_combDE$expressionSet)$Type_step2)==0))
  selected_models = which(rowSums(gene_domain_comb_tbl[,empty_types, drop=F]) == 0)
  gene_domain_comb_tbl2 <- gene_domain_comb_tbl[selected_models, ]

  return(gene_domain_comb_tbl2)

}

runCombinatorialDETests <- function(m_combDE, combinatorialModels, numcores=4){

  all_design_domain_level = combinatorialModels %>% dplyr::select(-genename) %>% dplyr::distinct(model, .keep_all=T)

  domain_cell_tbl = pData(m_combDE$expressionSet) %>% dplyr::select(Type_step2) %>% tibble::rownames_to_column("cellname") %>% as_tibble

  all_design_cell_level = all_design_domain_level %>%
                  tidyr::gather(domain, part_of, -model) %>%
                  dplyr::right_join(domain_cell_tbl, by=c("domain"="Type_step2")) %>% # right_join because some type2 may be missing in domain_cell_tbl$Type_step2
                  dplyr::select(-domain) %>%
                  tidyr::spread(model, part_of) %>%
                  dplyr::arrange(match(cellname, m_combDE$getCellsNames())) %>% # order as in original dataset
                  tibble::column_to_rownames("cellname") %>%
                  as.data.frame()

  all_design_cell_level$timepoint = pData(m_combDE$expressionSet)$timepoint # pData(m_combDE$expressionSet)$Type_step2

  HSMM <- monocle::newCellDataSet(
              as.matrix(m_combDE$getReadcounts(data_status='Raw')), #+.1,
              phenoData = new("AnnotatedDataFrame", data = all_design_cell_level),
              featureData = new('AnnotatedDataFrame', data = data.frame(row.names=m_combDE$getGeneNames(), "gene_short_name"=m_combDE$getGeneNames())),
              lowerDetectionLimit = .5,
              # expressionFamily=VGAM::tobit()
              # expressionFamily=VGAM::negbinomial() # appropriate for UMI count, but slower than negbinomial.size #http://cole-trapnell-lab.github.io/monocle-release/docs/
              expressionFamily=VGAM::negbinomial.size() # appropriate for UMI count #http://cole-trapnell-lab.github.io/monocle-release/docs/
              )

  # HSMM <- estimateSizeFactors(HSMM) # only with VGAM::negbinomial.size()
  # HSMM <- estimateDispersions(HSMM) # only with VGAM::negbinomial.size()

  model_test <- function(x){
    cbind(
      my_differentialGeneTest(HSMM[unlist(x$genelist), ,drop=F],
                                # fullModelFormulaStr=paste0("~", x$model, "+timepoint"),
                                fullModelFormulaStr=paste0("~", x$model),
                                reducedModelFormulaStr = "~1",
                                cores=1,
                                relative_expr=F, 
                                verbose=T),
      model=x$model
      )
  }

  cl <- multidplyr::create_cluster(numcores) # too many cores trigger segmentation fault (lack of memory when HSMM is copied on each cluster)
  multidplyr::cluster_copy(cl, model_test)
  multidplyr::cluster_copy(cl, HSMM)

  multidplyr::cluster_assign_value(cl, "my_diff_test_helper", my_diff_test_helper)
  multidplyr::cluster_assign_value(cl, "my_differentialGeneTest", my_differentialGeneTest)
  multidplyr::cluster_assign_value(cl, "my_compareModels", my_compareModels)
  multidplyr::cluster_assign_value(cl, "my_smartEsApply", my_smartEsApply)

  cl %>% multidplyr::cluster_ls()

  # multidplyr::cluster_rm(cl, c('my_diff_test_helper', 'my_differentialGeneTest', 'my_compareModels', 'my_smartEsApply'))

  t0 = Sys.time()
  Monocle_run = combinatorialModels %>%
                      dplyr::group_by(model) %>%
                      dplyr::summarise(genelist = list(genename)) %>%
                      multidplyr::partition(model, cluster=cl) %>%
                      dplyr::do(model_test(.)) %>%
                      dplyr::collect()
  t1 = Sys.time()

  print(t1 - t0)

  parallel::stopCluster(cl)

  return(Monocle_run)
}

getTopModel <- function(currm, mrun, temporal_pattern_thres=.2){

  # # A few models may produce infinite likelihood (for genes expressed in very small population, eg without resampling)
  # Monocle_run %>% dplyr::filter(is.infinite(full_model_loglikelihood))

  # These models are good ones, we change loglikelihood value to 0 -> top loglik
  Monocle_run2 = mrun %>% dplyr::mutate(full_model_loglikelihood = replace(full_model_loglikelihood, is.infinite(full_model_loglikelihood), 0))
  # Monocle_run2 %>% dplyr::filter(full_model_loglikelihood==0) # check replacement

  # Some test may have failed (ill-defined models where one of the 2 samples is empty)
  Monocle_run2 <- Monocle_run2 %>% dplyr::filter(!is.na(full_model_loglikelihood))

  # pop_sign: 1 for positive, 0 for negative
  # type_step1: "All", "Progenitor", or "Neuron"
  mean_level_domains <- function(x, pop_sign=1, type_step1="All"){
    domain_list = levels(pData(currm$expressionSet)$Type_step2)[which(strsplit(substring(x$model, 2), split="")[[1]]==pop_sign)]
    cell_list = which(pData(currm$expressionSet)$Type_step2 %in% domain_list)
    if(type_step1 %in% c("Progenitor", "Neuron")){
      cell_list <- cell_list[which(pData(currm$expressionSet)$Type_step1[cell_list]==type_step1)]
    } 
    mean(currm$getReadcounts('Raw')[x$genename, cell_list, drop=F])
  }

  ratio_expr <- function(x, pop_sign=1, type_step1="All"){
    domain_list = levels(pData(currm$expressionSet)$Type_step2)[which(strsplit(substring(x$model, 2), split="")[[1]]==pop_sign)]
    cell_list = which(pData(currm$expressionSet)$Type_step2 %in% domain_list)
    if(type_step1 %in% c("Progenitor", "Neuron")){
      cell_list <- cell_list[which(pData(currm$expressionSet)$Type_step1[cell_list]==type_step1)]
    } 
    mean(1*(currm$getReadcounts('Raw')[x$genename, cell_list, drop=F]>0))
  }

  time_correlation_pos_domains <- function(x, type_step1="All"){

    domain_list = levels(pData(currm$expressionSet)$Type_step2)[which(strsplit(substring(x$model, 2), split="")[[1]]==1)]
    cell_list = which(pData(currm$expressionSet)$Type_step2 %in% domain_list)
    if(type_step1 %in% c("Progenitor", "Neuron")){
      cell_list <- cell_list[which(pData(currm$expressionSet)$Type_step1[cell_list]==type_step1)]
    }

    cor(
          currm$getReadcounts('Raw')[x$genename, cell_list, drop=T],
          pData(currm$expressionSet)$timepoint[cell_list],
          method="spearman"
          )
  }

  #' Select top model and calculate associated metrics (fold-change, mean level in both clusters...)
  Monocle_run2 %>%
        dplyr::ungroup() %>%
        # dplyr::mutate(num_domain=sum(as.integer(strsplit(model, split="_")[[1]]))) %>% # count number of domains in model (for filter below)
        dplyr::mutate(genename=as.character(gene_short_name),
                      gene_short_name=NULL) %>%
        dplyr::group_by(genename) %>%
        dplyr::filter(full_model_loglikelihood >= 1.00001 * max(full_model_loglikelihood)) %>% # we add this step before adding metrics to reduce computation time. Also full_model_loglikelihood is negative.
        dplyr::ungroup() %>%
        # dplyr::arrange(num_domain, model) %>%
        dplyr::rowwise() %>%
        dplyr::do({res=dplyr::as_data_frame(.)      # add time correlation
                   res$tc = time_correlation_pos_domains(res, type_step1="All")
                   res$mean_pos = mean_level_domains(res, pop_sign=1, type_step1="All")
                   res$mean_neg = mean_level_domains(res, pop_sign=0, type_step1="All")
                   res$ratio_expr_pos = ratio_expr(res, pop_sign=1, type_step1="All")
                   res$ratio_expr_neg = ratio_expr(res, pop_sign=0, type_step1="All")
                   # res$mean_prog_pos = mean_level_domains(res, pop_sign=1, type_step1="Progenitor")
                   # res$mean_prog_neg = mean_level_domains(res, pop_sign=0, type_step1="Progenitor")
                   # res$mean_neuron_pos = mean_level_domains(res, pop_sign=1, type_step1="Neuron")
                   # res$mean_neuron_neg = mean_level_domains(res, pop_sign=0, type_step1="Neuron")
                   res$log2fc = log2(res$mean_pos) - log2(res$mean_neg)
                   # res$log2fc_prog = log2(res$mean_prog_pos) - log2(res$mean_prog_neg)
                   # res$log2fc_neuron = log2(res$mean_neuron_pos) - log2(res$mean_neuron_neg)
                   res
                   }) %>%
        dplyr::mutate(temporal_pattern=factor(
                        ifelse(tc >= temporal_pattern_thres,
                               "up",
                               ifelse(tc <= -temporal_pattern_thres,
                                      "down",
                                      "Undetermined")),
                        levels=c("down", "up", "Undetermined"))
                        ) %>% 
        dplyr::ungroup() %>%
        dplyr::group_by(genename) %>%
        dplyr::filter(log2fc > 0) %>% # some model pairs have exactly the same pval/qval because they are mirrored. we keep the "positive" model
        dplyr::filter(full_model_loglikelihood == max(full_model_loglikelihood)) %>% # then we filter the top model 
        dplyr::ungroup()

}


getModelAnnotations <- function(currm, model_list=unique(Monocle_run_topModel$model)){

  as_tibble(do.call(rbind, lapply(model_list, function(model){

              data.frame(
                    model=model,
                    model_names=paste0(levels(pData(currm$expressionSet)$Type_step2)[which(strsplit(substring(model, 2), split="")[[1]]==1)], collapse='-'),
                    model_prog_names=paste0(intersect(levels(pData(currm$expressionSet)$Type_step2)[which(strsplit(substring(model, 2), split="")[[1]]==1)], pData(currm$expressionSet) %>% dplyr::filter(Type_step1=="Progenitor") %>% .$Type_step2 %>% unique), collapse='-'),
                    model_neuron_names=paste0(intersect(levels(pData(currm$expressionSet)$Type_step2)[which(strsplit(substring(model, 2), split="")[[1]]==1)], pData(currm$expressionSet) %>% dplyr::filter(Type_step1=="Neuron") %>% .$Type_step2 %>% unique), collapse='-'),
                    first_domain_factor=factor(levels(pData(currm$expressionSet)$Type_step2)[which(strsplit(substring(model, 2), split="")[[1]]==1)[1]], levels=levels(pData(currm$expressionSet)$Type_step2)),
                    first_neuron_domain_factor=factor(levels(
                        pData(currm$expressionSet)$Type_step2)[
                          which(strsplit(substring(model, 2), split="")[[1]]==1 &
                                (pData(currm$expressionSet) %>% as_tibble %>% dplyr::select(Type_step1, Type_step2) %>% dplyr::distinct() %>% .$Type_step1 == "Neuron")
                                )[1]],
                        levels=levels(pData(currm$expressionSet)$Type_step2)),
                    num_domain=sum(as.integer(strsplit(substring(model, 2), split="")[[1]])),
                    num_prog_domain=sum(as.integer(strsplit(substring(model, 2), split="")[[1]][seq(pData(currm$expressionSet) %>% dplyr::filter(Type_step1=="Progenitor") %>% .$Type_step2 %>% unique %>% length)])),
                    num_neuron_domain=sum(tail(as.integer(strsplit(substring(model, 2), split="")[[1]]), n=pData(currm$expressionSet) %>% dplyr::filter(Type_step1=="Neuron") %>% .$Type_step2 %>% unique %>% length)),
                    broad_strict_ventral_prog_domain=grepl("^m000000011111[0|1]", model),  # RP-dp1-dp2 or dp1-dp2
                    broad_strict_dorsal_prog_domain=grepl("^m[0|1]110000000000", model),  # RP-dp1-dp2 or dp1-dp2
                    broad_strict_intermediate_prog_domain=grepl("^m0001111000000", model),
                    broad_ventral_prog_domain=grepl("^m0000000[0|1][0|1][0|1][0|1][0|1][0|1]", model),  # RP-dp1-dp2 or dp1-dp2
                    broad_dorsal_prog_domain=grepl("^m[0|1][0|1][0|1]0000000000", model),  # RP-dp1-dp2 or dp1-dp2
                    broad_intermediate_prog_domain=grepl("^m000[0|1][0|1][0|1][0|1]000000", model),
                    broad_prog_domain=factor(
                                        mapply(
                                        function(x){
                                          if(grepl("^m000[0|1][0|1][0|1][0|1]000000", model)){
                                            "Intermediate"
                                          } else if(grepl("^m0000000[0|1][0|1][0|1][0|1][0|1][0|1]", model)) {
                                            "Ventral"
                                          } else if(grepl("^m[0|1][0|1][0|1]0000000000", model)) {
                                            "Dorsal"
                                          } else {
                                            "NA"
                                          }
                                        },
                                        model),
                                      levels=c("Dorsal", "Intermediate", "Ventral", "NA")
                                      )
                )
    })))
}

makeGeneModules <- function(mrt, models_annotations, 
                            pval_thres = 1e-9,
                            pos_pop_level_threshold = .5,
                            pos_pop_ratio_threshold=.15,
                            log2fc_thres = 2,
                            gm_grouping_vars = c("model"),
                            gm_ordering_vars = c("first_domain_factor", "num_domain")
                            ){

  level_filtered_genes = mrt %>% dplyr::filter(mean_pos >= pos_pop_level_threshold) %>% .$genename %>% unique%>% sort
  ratio_filtered_genes = mrt %>% dplyr::filter(ratio_expr_pos >= pos_pop_ratio_threshold) %>% .$genename %>% unique%>% sort
  fc_filtered_genes = mrt %>% dplyr::filter(log2fc > log2fc_thres | mean_neg == 0) %>% .$genename %>% unique%>% sort
  selected_genes = Reduce(intersect, list(level_filtered_genes, fc_filtered_genes, ratio_filtered_genes))

  Monocle_run_filtered = mrt %>%
                            dplyr::filter(genename %in% selected_genes) %>%
                            dplyr::filter(pval <= pval_thres)
  
  Monocle_run_genelist = Monocle_run_filtered %>%
          dplyr::group_by_(.dots=gm_grouping_vars) %>% 
          dplyr::summarise(genelist = list(genename)) %>%
          dplyr::ungroup() %>%
          dplyr::left_join(models_annotations, by="model")
  
  if("temporal_pattern" %in% gm_grouping_vars){
    Monocle_run_genelist <- Monocle_run_genelist %>% tidyr::unite(model_names, temporal_pattern, col="model_names", sep=" / ", remove=FALSE)
  }

  Monocle_run_genelist %>%
        dplyr::arrange_at(gm_ordering_vars) %>%
       {"names<-"(.$genelist, .$model_names)}
}

plotPredictedSelection <- function(
                            currm,
                            genelist,
                            selected_genes,
                            monocle_run_topModel,
                            pc = getPopulationColors(known_template_file=partition_path)$Step2,
                            max_cells=2000,
                            zscore=FALSE,
                            ...
                        ){

  counts = currm$getReadcounts('Raw')[unlist(genelist), ]

  if(zscore){
    counts <- counts[which(apply(counts, 1, sd) > 0), ]
    counts <- t(scale(t(counts), scale=T, center=T))
  }

  cat_level = tibble::as_tibble(counts) %>%
                    dplyr::mutate('genename'=rownames(counts)) %>%
                    tidyr::gather(cellname, count, -genename) %>%
                    dplyr::mutate(logcount=log(1+count)) %>%
                    dplyr::left_join(cbind.data.frame('cellname'=currm$getCellsNames(), Type_step2=pData(currm$expressionSet)$Type_step2, timepoint=pData(currm$expressionSet)$timepoint), by="cellname")

  # Subsample cells randomly for each gene (population proportion unchanged) to allow density visualization
  if(!is.na(max_cells)){

    if(max_cells > currm$getNumberOfCells()){
      stop("max_cells must be less than the number of cells...")
    }

    cat_level = cat_level %>%
                      dplyr::group_by(genename) %>%
                      dplyr::sample_n(max_cells, replace=T)
  }

  null_model = paste0("m", paste0(rep(0, length(levels(pData(currm$expressionSet)$Type_step2))), collapse=''))

  cat_level.3 = cat_level %>%
                    dplyr::left_join(
                        monocle_run_topModel %>%
                            dplyr::select(genename, model) %>%
                            dplyr::group_by(genename) %>%
                            dplyr::mutate(model_sel = ifelse(genename %in% selected_genes,
                                                             model,
                                                             null_model
                                                             )
                                        ),
                        by="genename"
                      ) %>%
                    dplyr::group_by(Type_step2, model_sel) %>% # speed up diff_expr calculation
                    # dplyr::rowwise() %>%
                    dplyr::mutate(
                        diff_expr = ifelse(Type_step2 %in% levels(pData(currm$expressionSet)$Type_step2)[which(strsplit(substring(model_sel, 2), split="")[[1]]==1)], as.character(Type_step2), "No")
                    ) %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(
                        diff_expr=factor(diff_expr, levels=c("No", levels(pData(currm$expressionSet)$Type_step2))),
                        genename=factor(genename, levels=unlist(genelist)),
                        timepoint=factor(timepoint, levels=seq(9.5, 13.5))
                        )

  for(cat_name in names(genelist)){
    
    print(cat_name)
    gl = genelist[[cat_name]]
    plot.df = cat_level.3[which(cat_level.3$genename %in% unlist(gl)),]

    # subplot
    p <- ggplot(plot.df) +
            # geom_violin(aes(x=timepoint, y=logcount, fill=Type_step2, color=factor(diff_expr=="No"))) +
            geom_point(aes(x=cellname, y=logcount, colour=diff_expr), size=.1) +
            # geom_jitter(aes(x=1, y=logcount, colour=diff_expr), size=.2) +
            coord_flip() +
            facet_grid("Type_step2 ~ genename", scales = "free_y", space = "free_y", shrink=FALSE) +
            scale_color_manual(
                values=c("gray70", unname(pc[levels(pData(currm$expressionSet)$Type_step2)])),
                breaks=c("No", levels(pData(currm$expressionSet)$Type_step2))
                ) +
            ylab("") + xlab('') + theme_void() +
            theme(
                  strip.text.x = element_text(angle = 90, vjust = 0),
                  strip.text.y = element_text(angle = 0, hjust = 0),
                  legend.position="none"
              )
    pdf(paste0(plot_path, '/', cat_name, '.pdf'), ...)
    print(p)
    graphics.off()
    
    # subplot
    p <- ggplot(plot.df) +
            # geom_violin(aes(x=timepoint, y=logcount, fill=Type_step2, color=factor(diff_expr=="No"))) +
            geom_point(aes(x=cellname, y=logcount, colour=diff_expr), size=.1) +
            # geom_jitter(aes(x=1, y=logcount, colour=diff_expr), size=.2) +
            coord_flip() +
            facet_grid("Type_step2 ~ genename") +
            scale_color_manual(
                values=c("gray70", unname(pc[levels(pData(currm$expressionSet)$Type_step2)])),
                breaks=c("No", levels(pData(currm$expressionSet)$Type_step2))
                ) +
            ylab("") + xlab('') + theme_void() +
            theme(
                  strip.text.x = element_text(angle = 90, vjust = 0),
                  strip.text.y = element_text(angle = 0, hjust = 0)
              )

    pdf(paste0(plot_path, '/', cat_name, '_shrink.pdf'), ...)
    print(p)
    graphics.off()

  }
}


plotAveragePerTypeTimeDV <- function(currm, gene, threshold = 10, basename=NULL) {

  order.cell.types <- c("FP", "p3", "pMN", "p2", "p1", "p0", "dp6", "dp5", "dp4", "dp3", "dp2", "dp1", "RP")

  dot.input = cbind.data.frame(pData(currm$expressionSet), glevel=currm$getReadcounts('Raw')[gene,]) %>%
              group_by(timepoint, DV, Type_step1) %>%
              dplyr::summarize(mean=mean(glevel), count=n()) %>%
              dplyr::mutate(DV_name=factor(order.cell.types[DV], levels=order.cell.types))

  dot.input[which(dot.input$count < threshold), "mean"] <- NA

  gg <- ggplot(data = dot.input, aes(x = factor(timepoint), y = DV_name)) +
      geom_count(aes(size = mean, color=factor(timepoint))) +
      scale_colour_manual(values=colorRampPalette(c("pink", "red"))(n=5)) +
      scale_size_area(max_size = 20) +
      theme_bw() +
      facet_wrap(~ unlist(dot.input$Type_step1)) +
      ylab("DV position") +
      xlab("embryonic days") +
      ggtitle(gene) +
      theme(axis.text.y = element_text(size = 24),
            axis.title.y = element_text(size = 28),
            axis.text.x = element_text(size = 24),
            axis.title.x = element_text(size = 28),
            legend.title = element_text(size = 20),
            legend.title.align=0.5,
            legend.text = element_text(size = 20),
            strip.text = element_text(size = 24, face = "bold"),
            plot.title = element_text(face = "bold", size = 32)) + 
      guides(color = 'none')
  
  pdf(paste0(plot_path, '/', basename, gene, ".pdf"), width =10, height = 10, useDingbats = FALSE)
  print(gg)
  graphics.off()

}



plotNeuroSpace <- function(
        cell_domains=2:12,
        genelist=c('Lin28a', 'Fabp7', 'Sox2', "Tubb3", "Elavl3"),
        basename=NULL){

  cell_ids = which(pData(m_neural$expressionSet)$DV %in% cell_domains)

  png(paste0(plot_path, '/', basename, 'PCA_Dotplot.pdf'), width=4000, height=4000)

  nplot = length(genelist) + 2
  ncol = as.integer(sqrt(nplot))+1
  nrow = as.integer(ceiling(nplot/ncol))

  par(mfrow = c(nrow, ncol),
      mar = c(3, 3, 3, 3)) 

  set.seed(1)
  random_order = sample(cell_ids)

  plot(pca_res_all_x[random_order, 1:2], col=pData(m_neural$expressionSet)[random_order, 'cells_colors'], pch=16, cex=3, main="Sample", xaxt='n', yaxt='n', asp=1, ann=T)

  pt_values = pca_res_all_x[, 1]

  # Pseudotime on PCA plot  
  plot(pca_res_all_x[random_order, 1:2],
        col=pt_values[random_order] %>% {1+as.integer(100*(.-min(.))/(max(.)-min(.)))} %>% colorRampPalette(c("gray80", "black"))(n = 101)[.],
        pch=16, cex=3, main="Monocole Pseudotime", xaxt='n', yaxt='n', asp=1, ann=T)

  for(gn in genelist){

    percentile_90 = quantile(log(1+m_neural$getReadcounts(data_status="Raw")[gn, random_order]), .9)
    if(percentile_90 == 0){
      percentile_90 <- 1
    }

    gene_cols= log(1+m_neural$getReadcounts(data_status="Raw")[gn, random_order]) %>% {1+as.integer(100*(./percentile_90))} %>% colorRampPalette(c("#0464DF", "#FFE800"))(n = 101)[.]

    plot(
          pca_res_all_x[random_order, 1:2], asp=1, col=gene_cols,
          pch=16, cex=3, xlab="", ylab="", main=gn, xaxt='n', yaxt='n')
  }

  dev.off()

  # Plot average level on hexagonal grid
  library(hexbin)
 
  h <- hexbin(x=pca_res_all_x[cell_ids,1], y=pca_res_all_x[cell_ids,2], xbins=20, IDs = TRUE)

  gg_input = cbind.data.frame(
            x=pca_res_all_x[cell_ids,1],
            y=pca_res_all_x[cell_ids,2],
            t(m_neural$getReadcounts(data_status="Raw")[genelist,cell_ids]),
            cID=h@cID,
            PT=pt_values[cell_ids],
            SampleTime=pData(m_neural$expressionSet)$timepoint[cell_ids]) %>%
        tidyr::gather(gene, val, -x, -y, -cID)

  rplots <- lapply(c("SampleTime", "PT", genelist), function(gn){
    print(gn)
    gg_input2 = gg_input %>%
                  dplyr::filter(gene==gn) %>%
                  group_by(cID) %>%
                  summarise(mean=mean(val)) %>%
                    dplyr::left_join(
                        cbind.data.frame(cID=h@cell, x=hcell2xy(h)$x, y=hcell2xy(h)$y),
                        by = "cID"
                        )

    if(gn %in% c("PT")) {
      lowcol = "gray80"
      highcol = "black"
    } else if (gn %in% c("SampleTime")) {
      lowcol = "pink"
      highcol = "red"
    } else {
      gg_input2$mean <- log(1+gg_input2$mean)
      lowcol = "#0464DF"
      highcol = "#FFE800"
    }

    ggplot(data=gg_input2) +
      geom_hex(aes(x=x, y=y, fill=mean), col="gray70", stat="identity", size=0) + 
      scale_fill_gradient(low=lowcol, high=highcol) +
      ggtitle(gn)  + xlab("") + ylab("") + 
      theme_void() + theme(legend.position="none", plot.title = element_text(size = 80))

  })

  pdf(paste0(plot_path, '/', basename, 'PCA_HEXplot.pdf'), height=40, width=40, useDingbats=FALSE)
  gridExtra::grid.arrange(grobs=rplots, layout_matrix=matrix(seq(3*(as.integer((length(genelist)+2)/3)+1)), ncol=3, byrow=T))
  graphics.off()

}



plotNeurogenesisHeatmap <- function(df.merge, num_pt, plotted_DV, pt_clusters_reordered, plotfullname, zscore_cap = 3, dv_colors) {

 
  col_ordered = paste0("1_", seq(num_pt))
  col_ordered_all <- unlist(lapply(plotted_DV, function(dv) paste0(dv, "_", col_ordered)))
  pt_cluster_id_from_gene_names = setNames(
                rep(seq(length(pt_clusters_reordered)), lapply(pt_clusters_reordered, length)),
                unlist(pt_clusters_reordered)
            )


  if(!is.null(zscore_cap)){

    # ! Zscore dataset with selected domains only
    bif_dataset_whole = t(apply(df.merge[unlist(pt_clusters_reordered), col_ordered_all], 1, function(x){x-mean(x, na.rm=T)}))
    bif_dataset_whole <- t(apply(bif_dataset_whole, 1, function(x){x/sd(x, na.rm=T)}))

     bif_dataset_whole <- t(apply(bif_dataset_whole, 1, function(x){
                x[which(x < -zscore_cap)] <- -zscore_cap
                x[which(x > zscore_cap)] <- zscore_cap
                x
    }))   
     col=colorRamps::blue2green2red(1000)
  } else {

    # ! Zscore dataset with selected domains only
    bif_dataset_whole = df.merge[unlist(pt_clusters_reordered), col_ordered_all]
    bif_dataset_whole <- t(apply(bif_dataset_whole, 1, function(x){
                              # q=quantile(x, .95)
                              # if(q==0)
                              #   q=1
                              # x/q
                              x/max(x, na.rm=T)
                          }))
    col=colorRampPalette(c("#0464DF", "#FFE800"))(n = 1000)

  }
  
  pdf(plotfullname, width=6, height=as.integer(0.01*length(unlist(pt_clusters_reordered)) + 6))
  heatmap.3(
    bif_dataset_whole,
    Colv = FALSE,
    Rowv = FALSE, 
    dendrogram="none",
    col=col,
    scale='none',
    margins=c(.5,.1),
    ColSideColors=cbind(
                      rep(colorRampPalette(c("gray80", "black"))(n=length(col_ordered)), length(plotted_DV)),
                      rep(rev(dv_colors[plotted_DV]), each=length(col_ordered))
                      ),
    RowSideColors=t(getClusterColors()[pt_cluster_id_from_gene_names]),
    labRow=unlist(pt_clusters_reordered),
    labCol=F,
    lmat = rbind(c(6, 5, 0), c(0,2,0), c(4, 3, 1)),
    lhei = c(.1, 1, 10),
    lwid = c(.1, 10, 1),
    key=F,
    rowsep=c(0, cumsum(unlist(lapply(pt_clusters_reordered, length)))),
    # colsep=NULL,
    sepcolor = "black",
    sepwidth = c(0.001, 0.001),
    cexRow = 3*(0 + .5/log10(nrow(bif_dataset_whole))),
    rect_overlay = c(0, cumsum(table(sort(rep(
                      # sort(unique(bif_dataset$DV)),
                      plotted_DV,
                      length.out=length(col_ordered_all)
                      ))))) %>%
                {lapply(seq(length(.)-1), function(i)
                      list(
                            xleft=.[[i]] + 0.5,
                            xright=.[[i+1]] + 0.5,
                            ytop=0.5,
                            ybottom=length(unlist(pt_clusters_reordered)) + 0.5,
                            width=.2
                            )
                  )}
    )
  graphics.off()
}

