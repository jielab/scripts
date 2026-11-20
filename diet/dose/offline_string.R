# Read only the pre-downloaded human STRING v12.0 reference database.
offline_string_network <- function(genes, directory, cutoff=.4, simulations=1000L, seed=12345L) {
    files <- setNames(file.path(directory, paste0('9606.protein.', c('info','aliases','links'), '.v12.0.txt.gz')),
                      c('info','aliases','links'))
    if (any(!file.exists(files))) stop('Missing local STRING files: ', paste(files[!file.exists(files)], collapse=', '))
    read <- function(path) as.data.frame(data.table::fread(path, showProgress=FALSE, nThread=1L))
    info <- read(files[['info']]); names(info)[1L] <- 'stringId'
    if (!all(c('stringId','preferred_name') %in% names(info)) || anyDuplicated(info$stringId)) stop('Invalid STRING info table')
    genes <- sort(unique(toupper(trimws(genes))))
    aliases <- read(files[['aliases']]); names(aliases)[1L] <- 'stringId'
    if (!all(c('stringId','alias') %in% names(aliases))) stop('Invalid STRING aliases table')
    aliases$alias <- toupper(aliases$alias)
    mapping <- lapply(seq_along(genes), function(i) {
        gene <- genes[i]
        # Logical NA subscripts introduce phantom IDs; only retain actual matches.
        ids <- unique(info$stringId[which(toupper(info$preferred_name) == gene)])
        if (!length(ids)) ids <- unique(aliases$stringId[which(aliases$alias == gene & aliases$stringId %in% info$stringId)])
        if (length(ids)>1L) stop('Ambiguous local STRING mapping for ',gene,': ',paste(ids,collapse=', '))
        if (!length(ids)) return(NULL)
        data.frame(queryIndex=i-1L, Gene=gene, stringId=ids, preferredName=info$preferred_name[match(ids,info$stringId)], ncbiTaxonId=9606L)
    })
    mapping <- do.call(rbind,mapping)
    if (is.null(mapping) || !nrow(mapping)) stop('No candidate genes map to local STRING reference')
    if (anyDuplicated(mapping$stringId)) stop('Multiple candidate genes map to the same STRING protein')
    rm(aliases)
    links <- data.table::fread(files[['links']], showProgress=FALSE, nThread=1L)
    if (!all(c('protein1','protein2','combined_score') %in% names(links))) stop('Invalid STRING links table')
    # Official files contain both AB and BA; retain each undirected edge once.
    links <- links[combined_score > cutoff*1000 & protein1 < protein2]
    if (any(!links$protein1 %in% info$stringId) || any(!links$protein2 %in% info$stringId)) stop('STRING link endpoints absent from info table')
    links <- unique(links, by=c('protein1','protein2'))
    ids <- mapping$stringId
    selected <- links[protein1 %in% ids & protein2 %in% ids]
    edges <- data.frame(stringId_A=selected$protein1, stringId_B=selected$protein2, score=selected$combined_score/1000)
    # Reproducible local random-set test, explicitly distinct from the online API test.
    graph <- igraph::graph_from_data_frame(as.data.frame(links[,c('protein1','protein2')]), directed=FALSE,
                                          vertices=data.frame(name=info$stringId))
    set.seed(seed)
    null <- replicate(simulations, igraph::ecount(igraph::induced_subgraph(graph, sample(info$stringId,length(ids)))))
    ppi <- data.frame(number_of_nodes=length(ids),number_of_edges=nrow(edges),
        expected_number_of_edges=mean(null),p_value=(1+sum(null>=nrow(edges)))/(simulations+1),
        method='Offline uniform random protein-set permutation; not STRING API PPI enrichment',
        background_proteins=nrow(info), simulations=simulations, seed=seed)
    list(mapping=mapping,edges=edges,ppi=ppi,
         manifest=data.frame(file=files,md5=unname(tools::md5sum(files))))
}
