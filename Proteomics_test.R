# Functions to Preprocess and Analyse Proteomics Experiments.

# Preprocess. (Tasks Zero Imputation and Log-Transformation)

impute_lognorm <- function(dataset, groupings, discard, discard_missing){
    # Args:
    #   dataset: A proteins x samples matrix with the dataset to preprocess
    #   groupings:  A samples x groups matrix with values of 1 if a sample belongs to a group (e.g. cases/controls)
    #   discard: minimum number of missing value a group can have not to be discarded.
    #   discard_missing: minimum number of missing value a protein can have not to be discarded.
    # First create new matrix for preprocessed results.
    prep_dataset = matrix(0,dim(dataset)[1],dim(dataset)[2])
    trash =matrix(FALSE, dim(dataset)[1],dim(dataset)[2]) # labels of proteins with too many missing values 
    # do a for loop across each group and analyse. 
    for (i in 1:(dim(groupings)[2])){
       # find all columns in groupings i.
       col_i = which(groupings[,i] == 1.)
       # copy the columns in col_i into a matrix.
       data_coli = dataset[,col_i]
       # now let's see which we need to discard.
       zeroes = rowSums(data_coli == 0.) >= discard
       trash[,col_i] = zeroes
       keepers = which(!zeroes)
       # for loop for each protein that is not to be discarded.
       for (j in 1:length(keepers)){
           prot_vals = data_coli[keepers[j],] # copy row.
           # check if there are missing values. 
           if (sum(prot_vals==0.)>0){
               nzs_j = which(prot_vals>0)
               zs_j = which(prot_vals==0.)
               prot_vals[nzs_j ] = log(prot_vals[nzs_j ]) # transform them.
               mu = mean(prot_vals[nzs_j ])# get mean. 
               std = sd(prot_vals[nzs_j ])# get std.
               sim_vals = rnorm(length(zs_j),mu,std)# simulate new data.
               prot_vals[zs_j] = sim_vals 
           } else {
               prot_vals = log(prot_vals)
           }
           prep_dataset[keepers[j],col_i] = prot_vals
       }

    } 
    # return the matrix with the values that are never discarded.
    tot_zeros = rowSums(dataset==0)
    trash[tot_zeros>= discard_missing,] = TRUE
    kept = which(rowSums(trash)==0)
    prep_dataset  = prep_dataset[kept,]
    return(list(prep_dataset =  prep_dataset, kept = kept))
}


# Differentially Expressed Proteins (Taks: t-tests and multiple comparison adjustment)

diff_exp_prots <- function(prep_dataset,groupings,vs){
    # Args:
    #    prep_dataset: A dataset that is log-transformed and imputed.
    #    groupings: A samples x groups matrix with values of 1 if a sample belongs to a group (e.g. cases/controls)
    #    vs: A list with the two groups to compare with each other. 
    # Returns:
    #    a matrix whose columns represent:
    #        logfc: log 2 fold change.
    #        pvals: P-value of the t test
    #        qvals: Multiple Comparisons Adjusted P-Value (FDR method)
    # pre-allocate memory for each result.
    results = matrix(0,dim(prep_dataset)[1],3)
    for (i in 1:dim(prep_dataset)[1]){
        # get each in the comparison
        prot_1 = prep_dataset[i,groupings[,vs[1]]==1]
        prot_2 = prep_dataset[i,groupings[,vs[2]]==1]
        # calculate fold change.
        m1 = log(mean(exp(prot_1)),2)
        m2 =  log(mean(exp(prot_2)),2)
        results[i,1] = m2-m1
        # Perform t-test
        results[i,2] = t.test(prot_1,prot_2,alternative = "less")$p.value
    }
    # perform multiple comparisons adjustements.
    results[,3] = p.adjust(results[,2],'fdr')
    return(results)
}


diff_exp_lr <- function(prep_dataset,des_mat){
    # Args:
    #    prep_dataset: A dataset that is log-transformed and imputed.
    #    groupings: A samples x groups matrix with values of 1 if a sample belongs to a group (e.g. cases/controls)
    #    vs: A list with the two groups to compare with each other. 
    # Returns:
    #    a matrix whose columns represent:
    #        coefs: coefficient of the regression.
    #        pvals: P-value of the t test
    #        qvals: Multiple Comparisons Adjusted P-Value (FDR method)
    # pre-allocate memory for each result.
    coefs =  matrix(0,dim(prep_dataset)[1],dim(des_mat)[2])
    pvals = matrix(0,dim(prep_dataset)[1],dim(des_mat)[2])
    qvals = matrix(0,dim(prep_dataset)[1],dim(des_mat)[2])
    for (i in 1:dim(prep_dataset)[1]){
        # do linear regression. 
        y = prep_data[i,]
        lin_model <- lm(y~des_mat)
        pvals[i,] = summary(lin_model)$coefficients[2:dim(des_mat)[2],4]
        coefs[i,] = summary(lin_model)$coefficients[2:dim(des_mat)[2],1]
    }
    # perform multiple comparisons adjustements.
    for (j in 1:dim(res_mat)[2]){
        qvals[,j] = p.adjust(pvals[,j],'fdr')
    }
    return(list(pvals = pvals, qvals = qvals, coefs = coefs))
}


read_lfqfile <-function(file_name, groups_id){
    # Read Proteomics Data.
    dataset <- read.delim(file_name)
    # read column names.
    samples <- colnames(dataset)
    gene_names_col <- which(regexpr('Gene',samples)>-1)
    samples <- samples[-gene_names_col]
    # get gene names
    # pre-allocate matris for data.
    # Remove first row.
    dataset <- dataset[-1,]
    genes <- dataset[,gene_names_col]
    datmat <- matrix(0,length(genes),length(samples))
    start_ = "y."
    end_ = "_"
    grouping = matrix(0,length(samples),length(groups_id))
    for (i in 1:length(samples)){
        datmat[,i] = as.numeric(dataset[,i])
        # match grouping.
        oo = regexpr(start_,samples[i])[1]+2
        ff = regexpr(end_,samples[i])[1]-1
        id = substr(samples[i], oo,ff)
        grouping[i,which(id==groups_id)] = 1
    }
    return(list(datmat = datmat, genes = genes, grouping = grouping))
}

read_analyse<- function(file_path, groups_id, discard,  discard_missing, vs, outliers){
    # first read file. 
    resfile <- read_lfqfile(file_path,groups_id)
    if (!missing(outliers)){
        resfile$datmat = resfile$datmat[,-outliers]
        resfile$grouping = resfile$grouping[-outliers,]
    }
    # impute and log normalize
    resimp <- impute_lognorm(resfile$datmat, resfile$grouping, discard, discard_missing)
    # keep the genes
    genes_kept <- resfile$genes[resimp$kept]
    allres = c()
    for (i in 1:dim(vs)[1]){
        print(paste(paste(paste("Comparison", as.character(vs[i,1])),"vs"),as.character(vs[i,2])))
        restest <- diff_exp_prots(resimp$prep_dataset,resfile$grouping,vs[i,])
        rownames(restest) = genes_kept
        colnames(restest) = c("lfc", "pvalue","qvalue")
        allres[[i]] <- restest
    }
    return(list(prepdata = resimp$prep_dataset, genes = genes_kept, t.test.results = allres))
}



get_interactome <- function(results,sig,lfc){
    overlap = matrix(0,length(results$t.test.result),length(results$t.test.result))
    prots = c()
    for (i in 1:length(results$t.test.result)){
        # access results
        resi = results$t.test.result[[i]]  
        if (!is.null(dim(resi))){
           id  = which((resi[,3]<sig)&(resi[,1]>lfc))
            prots[[i]] =  rownames(resi)[id] # statistically significant.
        }
    }
    i = 1
    j = 1 
    while (i <= length(prots)){
        while (j <= length(prots)){
              string_res = as.character(sum(is.element(prots[[i]], prots[[j]]))) 
             string_res = paste(string_res , as.character(length(prots[[i]])), sep = "/")
             overlap[i,j] = string_res
           j = j +1             
         }
        i = i +1
        j = 1
    }
    if (length(prots)>1) {
         each_net = c()
        each_net[[1]] = prots[[1]]
        all_network = prots[[1]]
        for (i in 2:length(prots)){
            all_network = append(all_network, prots[[i]])
            each_net[[i]] = prots[[i]]
        }
    } else if (length(prots)==1){
        all_network = prots[[1]]
        each_net = c()
    }   
    return(list(prots = prots, overlap = overlap, all_network = unique(all_network), each_net = each_net))
}

read_linear_regression <- function(file_path, groups_id, design_mat){
# first read file. 
    resfile <- read_lfqfile(file_path,groups_id)
    # impute and log normalize
    resimp <- impute_lognorm(resfile$datmat, resfile$grouping, discard)
    # keep the genes
    genes_kept <- resfile$genes[resimp$kept]
    # pre-allocate memory for results. 
    # coefficients of regression, P-values 
    for (i in 1:length(genes_kept)){
        print(genes_kept[i])
        y = (resimp$prep_dataset[i,])
        restest <- diff_exp_prots(resimp$prep_dataset,resfile$grouping,vs[i,])
        rownames(restest) = genes_kept
        colnames(restest) = c("lfc", "pvalue","qvalue")
        allres[[i]] <- restest
    }
    return(list(prepdata = resimp$prep_dataset, genes = genes_kept, t.test.results = allres))
}

create_designmat <- function(tbl,grouping){
    # find unique variables.
    # fill in first column
    uni_tb = unique(tbl[,1])
    # make a design matrix
    des_j = model.matrix(~tbl[,1])
    desmat = des_j
    if (dim(tbl)[2]>1){
       for (j in 2:dim(tbl)[2]){
          des_j = model.matrix(~tbl[,2])
          desmat = cbind(desmat, des_j)
       }
    }
    # pre-allocate memory for the deisgnmat
    designmat = matrix(0,dim(grouping)[1],dim(desmat)[2])
    for (i in 1:dim(grouping)[1]){
       idx = which(grouping[i,]==1)
       designmat[i,] = desmat[idx,]
    }
    return(designmat)
}


overlap <- function(inter, nodes){
     pres_mat = matrix(0,length(nodes),length(inter$each_net))
     for (i in 1:length(inter$each_net)){
        w = which(is.element(nodes,inter$each_net[[i]]))
        pres_mat[w,i] =1
     }
     return(pres_mat)
}

table_2_adj <- function(interactome){
    #get unique set of nodes.
    node_1 = interactome[,1]
    node_2 = interactome[,2]
    nodes = unique(c(node_1,node_2))
    # pre-allocate memory for the adjancency matrix.
    adj = matrix(0,length(nodes),length(nodes))
    for (i in 1:length(nodes)){
       # find nodes that belong to the nodes[i]
       wj = which(interactome[,1]==nodes[i])
       wj2 = which(interactome[,2] == nodes[i])
       node_js = unique(c(node_2[wj],node_1[wj2]))
       j = which(is.element(nodes, node_js))
       adj[i,j] = 1 
       adj[j,i] = 1
    }
    return(list(adj =adj,nodes =  nodes))
}

laplacian <- function(adj){
      D = diag(rowSums(adj))
      L = D - adj
      return(L)
}
