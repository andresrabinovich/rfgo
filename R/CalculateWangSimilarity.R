#' Calculate pairwise Wang semantic similarity between a group of genes
#'
#' Calculate pairwise Wang semantic similarity between a group of genes as described in
#' James Z. Wang, Zhidian Du, Rapeeporn Payattakool, Philip S. Yu, Chin-Fu Chen,
#' A new method to measure the semantic similarity of GO terms,
#' Bioinformatics, Volume 23, Issue 10, 15 May 2007, Pages 1274–1281,
#' https://doi.org/10.1093/bioinformatics/btm087
#'
#' @param OrgDb OrgDb object
#' @param keytype keytype to select records from OrgDb. All possible keys are returned by using the keytypes method.
#' @param ont Type of ontology to be used. One of 'BP', 'MF', 'CC'. Defaults to BP.
#' @param genes List of genes for which to calculate pairwise Wang semantic similarity between them. Defaults to NULL (use all annotated genes)
#' @param dropCodes List of codes of annotation evidence types. Annotations with this codes will be dropped before calculating similiarity. Defaults to NULL, including all annotations.
#' @param minIC Minimum term information content 'IC' to take into account. Terms with less IC won't be used for calculations. Defaults to NULL
#' @keywords Wang similarity
#' @return A matrix of pairwise gene similarities
#' @export

CalculateWangSimilarity <- function(OrgDb, keytype = "ENTREZID", ont = "BP", genes = NULL, dropCodes = NULL, minIC = NULL){

  #Cargamos el objeto OrgDb
  OrgDb <- GOSemSim::load_OrgDb(OrgDb)

  #Sanity check
  if(!keytype %in% AnnotationDbi::keytypes(OrgDb)){
    stop(paste("keytype must be one of", paste(keytypes(OrgDb), collapse = ", ")))
  }
  if(is.null(dropCodes)){
    dropCodes <- c()
  }

  #Levantamos la data para calcular la similaridad
  semData <- GOSemSim::godata(OrgDb, keytype = keytype, ont = ont)

  #Si no pide algunos genes en particular, traemos todos los que estén anotados
  if(is.null(genes)){
    message("using all annotated genes")
    genes <- AnnotationDbi::keys(OrgDb, keytype = keytype)
  }

  message(paste("starting filtering with", length(genes), "genes"))

  #Levantamos los terminos en los que está anotado cada gen en una lista
  message("preparing annotations")
  gos        <- lapply(genes, GOSemSim:::gene2GO, godata = semData, dropCodes = dropCodes)
  names(gos) <- genes

  #Filtramos por minIC si es pertinente
  #minIC <- -log2(1/10)
  if(!is.null(minIC)){
    message(paste('filtering terms by IC larger than ', minIC))
    gos <- lapply(gos,function(x){
      y <- x[semData@IC[x]>minIC]
      y <- y[!is.na(y)]
      return(y)
    })
    message(paste(length(genes), "genes remaining after filtering by minIC"))
  }

  #Sacamos los genes que no esten anotados en ningun termino
  gos <- gos[unlist(lapply(gos, length)) > 0]

  if(length(gos) == 0){
    stop("No annotated gene remains after filtering. Try reducing minIC, drop less evidence or select more genes.")
  }

  message(paste(length(genes), "genes remaining after removing genes without annotations"))

  #Levantamos los scores semanticos
  message("preparing semantic scores, this might take a while, please wait...")
  uniqueGO                <- unique(unlist(gos))
  scores                  <- GOSemSim::termSim(uniqueGO, uniqueGO, semData = semData, method = "Wang")

  message(paste(nrow(scores), "genes remaining after filtering by genes whose scores can't be calculated"))

  similaridades           <- matrix(0, ncol = length(gos), nrow = length(gos))
  rownames(similaridades) <- colnames(similaridades) <- names(gos)

  #Calculamos la similaridad entre genes todos contra todos en el upper.tri + diag a partir de los scores
  message("Calculating semantic similarity, this might take a while, please wait...")
  pb <- txtProgressBar(min = 0, max = (length(gos)-1), style = 3)
  for(i in 1:length(gos)){ #hasta el max para capturar la diagonal
    setTxtProgressBar(pb, i)
    s_i <- scores[gos[[i]], , drop = F]
    for(j in i:length(gos)){ #desde i para capturar la diagonal
      s_ij <- s_i[, gos[[j]], drop = F]
      similaridades[i, j] <- (sum(apply(s_ij, 1, max)) +  sum(apply(s_ij, 2, max)))/(nrow(s_ij) + ncol(s_ij))
    }
  }
  #Terminamos de completar la matriz
  similaridades[lower.tri(similaridades)] <- t(similaridades)[lower.tri(similaridades)]

  #Devolvemos la matriz
  return(similaridades)
}


