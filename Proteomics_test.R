# Proteomics tests

source("F:/RAF_dimers/Proteomics.R", chdir = TRUE)

library(testthat)

# Test impute_lognorm

test_that("impute_lognorm_kept",{
    # simulate some toy data 
    # with 10 proteins and two groups
    # of 3 biological reps with 2 technical reps.
    # Pre-allocate data.
    toy_data = matrix(0,10,12)
    # create a groupings matrix.   
    groupings = matrix(0,12,2)
    groupings[1:6,1] = 1
    groupings[7:12,2] = 1
    # The data will be assumed to be log normal
    # mean and standard deviation group 1. 
    mu1 = 25
    sd1 = 1
    data_group1 = round(exp(rnorm(60,mu1,sd1)))
    # mean and standard deviation group 2.
    mu2 = 50
    sd2 = 1.5
    data_group2 = round(exp(rnorm(60,mu2,sd2)))
    # allocate the data. 
    toy_data[,groupings[,1]==1] = data_group1
    toy_data[,groupings[,2]==1] = data_group2
    # delete some observations (like misssing values). 
    toy_data[1,1] = 0
    toy_data[3,c(1,3)] = 0
    toy_data[5,c(1,4,5,10,12)] = 0
    toy_data[7,c(1,2,3,4,5,6,7)] = 0
    toy_data[9,c(2,10)] = 0
    # impute and log-transform.
    res_test = impute_lognorm(toy_data,groupings,4,4)
    # the ouput has two objects: 
    # kept: a list with the indexes of the proteins kept.
    # prep_dataset: a matrix with the imputed log normalized data.
    #### Test that we kept the right proteins 
    # Only the proteins 5 and 7 should have been deleted.
    expect_equal(6,res_test$kept[5])
    expect_equal(8,res_test$kept[6])
})

test_that("impute_lognorm_norm",{
    # simulate some toy data 
    # with 10 proteins and two groups
    # of 3 biological reps with 2 technical reps.
    # Pre-allocate data.
    toy_data = matrix(0,10,12)
    # create a groupings matrix.   
    groupings = matrix(0,12,2)
    groupings[1:6,1] = 1
    groupings[7:12,2] = 1
    # The data will be assumed to be log normal
    # mean and standard deviation group 1. 
    mu1 = 25
    sd1 = 1
    data_group1 = round(exp(rnorm(60,mu1,sd1)))
    # mean and standard deviation group 2.
    mu2 = 50
    sd2 = 1.5
    data_group2 = round(exp(rnorm(60,mu2,sd2)))
    # allocate the data. 
    toy_data[,groupings[,1]==1] = data_group1
    toy_data[,groupings[,2]==1] = data_group2
    # delete some observations (like misssing values). 
    toy_data[1,1] = 0
    toy_data[3,c(1,3)] = 0
    toy_data[5,c(1,4,5,10,12)] = 0
    toy_data[7,c(1,2,3,4,5,6,7)] = 0
    toy_data[9,c(2,10)] = 0
    # impute and log-transform.
    res_test = impute_lognorm(toy_data,groupings,4,4)
    # Now let's test the log normalization. 
    # the mean of the log res_test with no missing should be the same.
    actual_mp2 = mean(res_test$prep_dataset[2,1:6])
    expected_mp2 = mean(log(toy_data[2,1:6]))
    expect_equal(actual_mp2,expected_mp2)
    actual_mp3_n13 = mean(res_test$prep_dataset[3,c(2,4,5,6)])
    expected_mp3_n13 = mean(log(toy_data[3,c(2,4,5,6)]))
    expect_equal(actual_mp3_n13,expected_mp3_n13)
})

test_that("impute_lognorm_imputation",{
# simulate some toy data 
    # with 10 proteins and two groups
    # of 3 biological reps with 2 technical reps.
    # Pre-allocate data.
    toy_data = matrix(0,10,12)
    # create a groupings matrix.   
    groupings = matrix(0,12,2)
    groupings[1:6,1] = 1
    groupings[7:12,2] = 1
    # The data will be assumed to be log normal
    # mean and standard deviation group 1. 
    mu1 = 25
    sd1 = 1
    data_group1 = round(exp(rnorm(60,mu1,sd1)))
    # mean and standard deviation group 2.
    mu2 = 50
    sd2 = 1.5
    data_group2 = round(exp(rnorm(60,mu2,sd2)))
    # allocate the data. 
    toy_data[,groupings[,1]==1] = data_group1
    toy_data[,groupings[,2]==1] = data_group2
    # delete some observations (like misssing values). 
    toy_data[1,1] = 0
    toy_data[3,c(1,3)] = 0
    toy_data[5,c(1,4,5,10,12)] = 0
    toy_data[7,c(1,2,3,4,5,6,7)] = 0
    toy_data[9,c(2,10)] = 0
    # impute and log-transform.
    res_test = impute_lognorm(toy_data,groupings,4,4)
    # finally, let's test the imputed values. 
    # they should have similar mean to the mean of the non-missing ones.
    actual_mp3_13 = mean(res_test$prep_dataset[3,c(2,4,5,6)])
    expected_mp3_13 = mean(log(toy_data[3,c(2,4,5,6)]))
    dist = (actual_mp3_13 - expected_mp3_13)**2
    expect_lt(dist,2)
})

test_that("impute_lognorm_reordering",{
# simulate some toy data 
    # with 10 proteins and two groups
    # of 3 biological reps with 2 technical reps.
    # Pre-allocate data.
    toy_data = matrix(0,10,12)
    # create a groupings matrix.   
    groupings = matrix(0,12,2)
    groupings[1:6,1] = 1
    groupings[7:12,2] = 1
    # The data will be assumed to be log normal
    # mean and standard deviation group 1. 
    mu1 = 25
    sd1 = 1
    data_group1 = round(exp(rnorm(60,mu1,sd1)))
    # mean and standard deviation group 2.
    mu2 = 50
    sd2 = 1.5
    data_group2 = round(exp(rnorm(60,mu2,sd2)))
    # allocate the data. 
    toy_data[,groupings[,1]==1] = data_group1
    toy_data[,groupings[,2]==1] = data_group2
    # delete some observations (like misssing values). 
    toy_data[1,1] = 0
    toy_data[3,c(1,3)] = 0
    toy_data[5,c(1,4,5,10,12)] = 0
    toy_data[7,c(1,2,3,4,5,6,7)] = 0
    toy_data[9,c(2,10)] = 0
    # impute and log-transform.
    res_test = impute_lognorm(toy_data,groupings,4,4)
    # check no mistakes ocurred when re-ordering. 
    expect_equal(res_test$prep_dataset[8,1],log(toy_data[10,1]))
    expect_equal(res_test$prep_dataset[5,10],log(toy_data[6,10]))
})

test_that("impute_lognorm_dimensions",{
    # simulate some toy data 
    # with 10 proteins and two groups
    # of 3 biological reps with 2 technical reps.
    # Pre-allocate data.
    toy_data = matrix(0,10,12)
    # create a groupings matrix.   
    groupings = matrix(0,12,2)
    groupings[1:6,1] = 1
    groupings[7:12,2] = 1
    # The data will be assumed to be log normal
    # mean and standard deviation group 1. 
    mu1 = 25
    sd1 = 1
    data_group1 = round(exp(rnorm(60,mu1,sd1)))
    # mean and standard deviation group 2.
    mu2 = 50
    sd2 = 1.5
    data_group2 = round(exp(rnorm(60,mu2,sd2)))
    # allocate the data. 
    toy_data[,groupings[,1]==1] = data_group1
    toy_data[,groupings[,2]==1] = data_group2
    # delete some observations (like misssing values). 
    toy_data[1,1] = 0
    toy_data[3,c(1,3)] = 0
    toy_data[5,c(1,4,5,10,12)] = 0
    toy_data[7,c(1,2,3,4,5,6,7)] = 0
    toy_data[9,c(2,10)] = 0
    # impute and log-transform.
    res_test = impute_lognorm(toy_data,groupings,4,4)
    # check number of columns in the prep_dataset is the same. 
    expect_equal(dim(res_test$prep_dataset)[2],12)  
    expect_equal(dim(res_test$prep_dataset)[1],8)
})


# Test diff_exp_prots
test_that("diff_exp_prots_t_test",{
    # simulate some toy data 
    # with 10 proteins and two groups
    # of 3 biological reps with 2 technical reps.
    # Pre-allocate data.
    toy_data = matrix(0,2,12)
    # create a groupings matrix.   
    groupings = matrix(0,12,2)
    groupings[1:6,1] = 1
    groupings[7:12,2] = 1
    # The data will be assumed to be log normal
    # mean and standard deviation group 1. 
    mu1 = 25
    sd1 = 1
    data_group1 = round(exp(rnorm(6,mu1,sd1)))
    # mean and standard deviation group 2.
    mu2 = 50
    sd2 = 1.5
    data_group2 = round(exp(rnorm(6,mu2,sd2)))
    # allocate the data. 
    toy_data[1,groupings[,1]==1] = data_group1
    toy_data[1,groupings[,2]==1] = data_group2
    data_group3 = round(exp(rnorm(12,mu2,sd2)))
    toy_data[2,] = data_group3
    # impute and log-transform.
    res_test = impute_lognorm(toy_data,groupings,4,4)
    res = diff_exp_prots(res_test$prep_dataset,groupings,c(1,2))
    # expect to have P-values lower than 0.05 row 1.
    expect_lt(res[1,3],0.05)
    # expect to have P-values grater than 0.05 row 2.
    expect_gt(res[2,3],0.05)
    # check difference in means in Ok.
    expected_lfc = mean(log(data_group1,2)-log(data_group2,2))
    expect_equal(expected_lfc,res[1])
})
    
# test the function to read lfq file. 


test_that("read_lfqfile_dimensions",{
    # create grups ids, which are samples identifiers
    groups_id = c("1")
    for (i in 2:24){
        groups_id = append(groups_id,as.character(i))
    }
    resfile = read_lfqfile('F:/RAF dimers/jens_RAF.txt',groups_id)
    dataset <- read.delim('F:/RAF dimers/jens_RAF.txt')
    # check that the dimensions are correct. 
    expect_equal(dim(resfile$datmat)[1],length(resfile$genes))
})

test_that("read_lfqfile_values",{
    # create grups ids, which are samples identifiers
    groups_id = c("1")
    for (i in 2:24){
        groups_id = append(groups_id,as.character(i))
    }
    resfile = read_lfqfile('F:/RAF dimers/jens_RAF.txt',groups_id)
    dataset <- read.delim('F:/RAF dimers/jens_RAF.txt')
    # check that the same genes have the same values. 
    expect_equal(as.numeric(dataset[2,1]),resfile$datmat[1,1])
    # sample from the end.
    expect_equal(as.numeric(dataset[dim(dataset)[1],10]),resfile$datmat[dim(resfile$datmat)[1],10])
})

test_that("read_lfqfile_groups",{
    # create grups ids, which are samples identifiers
    groups_id = c("1")
    for (i in 2:24){
        groups_id = append(groups_id,as.character(i))
    }
    resfile = read_lfqfile('F:/RAF dimers/jens_RAF.txt',groups_id)
    dataset <- read.delim('F:/RAF dimers/jens_RAF.txt')
    # There's 6 replicates of group 1 so expect the right assignment.
    expect_equal(1,which(resfile$groupings[1,]==1))
    expect_equal(2,which(resfile$groupings[7,]==1))

})

test_that("read_analyse",{
    # simulate some toy data 
    # with 10 proteins and two groups
    # of 3 biological reps with 2 technical reps.
    # Pre-allocate data.
    toy_data = matrix(0,2,12)
    # create a groupings matrix.   
    groupings = matrix(0,12,2)
    groupings[1:6,1] = 1
    groupings[7:12,2] = 1
    # The data will be assumed to be log normal
    # mean and standard deviation group 1. 
    mu1 = 25
    sd1 = 1
    data_group1 = round(exp(rnorm(6,mu1,sd1)))
    # mean and standard deviation group 2.
    mu2 = 50
    sd2 = 1.5
    data_group2 = round(exp(rnorm(6,mu2,sd2)))
    # allocate the data. 
    toy_data[1,groupings[,1]==1] = data_group1
    toy_data[1,groupings[,2]==1] = data_group2
    data_group3 = round(exp(rnorm(12,mu2,sd2)))
    cols_1 = c("intensity 1_a","intensity 1_b","intensity 1_c","intensity 1_d","intensity 1_e","intensity 1_f")
    cols_2 = c("intensity 2_a","intensity 2_b","intensity 2_c","intensity 2_d","intensity 2_e","intensity 2_f")
    cols = c(cols_1,cols_2)
    cols = c(cols,"Gene")
    dataset_fake = matrix("E",3,13)
    dataset_fake[2,1:12] = as.character(toy_data[1,])
    dataset_fake[3,1:12] =as.character( toy_data[2,])
    dataset_fake[1,13] = "E"
    dataset_fake[2,13] = "BRAF"
    dataset_fake[3,13] = "TP53"
    colnames(dataset_fake) = cols
    # write the dataset
    write.table(dataset_fake,"dataset_test.txt", sep = "\t",quote = FALSE, row.names = FALSE)
    # read and analyse the dataset. 
    vs = matrix(0,1,2)
    vs[1] = 1
    vs[2] = 2
    res = read_analyse("dataset_test.txt",c("1","2"),4,4,vs)
    expect_lt(res[[3]][[1]][1,3],0.05)
    expect_gt(res[[3]][[1]][2,3],0.05)
    expect_equal(res[[2]][1],"BRAF")
})


    


