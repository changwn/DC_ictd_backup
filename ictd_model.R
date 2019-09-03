

print("test docker with Hello word!")

# Then, test the ICTD
library('ICTD')
library('readr')
#data_bulk = GSE72056_diri_example[[1]]
#print(data_bulk[1:5,1:6])
# ictd_result <- ICTD(data_bulk)
# ictd_result[[2]]
#print(sessionInfo())


#----------------function part---------------
ictd_2_output <- function(ictd_prop, dataset.name)
{
  #2
  col_name <- matrix(colnames(ictd_prop),1,length(colnames(ictd_prop)))
  vv <- c()
  for(i in 1:length(col_name))
  {
    vv <- c(vv, rep(col_name[i], 8))
  }
  sample_col <- vv
  
  #3  
  celltype_col <- rep(c('CD4.T.cells','CD8.T.cells','NK.cells','B.cells','monocytic.lineage','neutrophils','endothelial.cells','fibroblasts'), ncol(ictd_prop))
  #4
  prediction_col <- rep(runif(8,min=0.01,max=0.8), ncol(ictd_prop))  
  
  #bind
  output_df <- data.frame(dataset.name=dataset.name,sample.id=sample_col,cell.type=celltype_col,prediction=prediction_col)
  
  return(output_df)
}

ICTD_round1 <- function(data_bulk)
{
  data.matrix = data_bulk
  if (length(colnames(data.matrix)) == 0) {
    warning("input data do NOT have colnames")
    colnames(data.matrix) <- paste("Setsample", 1:ncol(data.matrix),sep = "")
  }
  data.matrix <- rm_zero_row(data.matrix)
  if (max(data.matrix) > 20) {
    d.matrix <- log(data.matrix + 1)
    d.matrix <- as.matrix(d.matrix)
  }else{
    d.matrix <- as.matrix(data.matrix)
  }
  data0 <- d.matrix
  data2 <- data0
  data21 <- data2[order(-apply(data2, 1, mean)), ]
  data22 <- data21[unique(rownames(data21)), ]
  data22_code <- data22[intersect(rownames(data22), TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3] == "protein_coding" & TCGA_ensem_annotation[,4] != ""), 4]), ]
  data23 <- data22[intersect(rownames(data22), TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3] == "protein_coding" & TCGA_ensem_annotation[,4] != ""), 4]), ]
  data23 <- normalize_data2(data23)
  data_CORS_cancer <- data23
  data_ccc <- data23
  list_c1 <- MRHCA_IM_compute_MR(data_CORS_cancer = data_ccc, IM_id_list, immune_cell_uni_table = immune_cell_uni_table0_GS)
  MR_IM_result_new_c <- MRHCA_IM_compute_full_pub_new(data_CORS_cancer = data_ccc, list_c = list_c1, IM_id_list, immune_cell_uni_table = immune_cell_uni_table0_GS)
  tg_key = "nonono"
  list_new_c2 <- Process_MR_IM_result_new(MR_IM_result_new_c, tg_key_c = tg_key, cell_type_enrich_cut = 0.4, cor_cut0 = 0.7, num_cut = 7, num_cut2 = 5, IM_id_list, immune_cell_uni_table = immune_cell_uni_table0_GS)
  R1_filter_step1_results_new <- R1_list_filtering_step1_new(list_new_c2, data_CORS_cancer = data_ccc, max_cut = 7, cutn0 = 7, cut10 = 0.7, IM_id_list, immune_cell_uni_table = immune_cell_uni_table0_GS)
  tg_R1_lists <- R1_filter_step1_results_new[[4]]
  print("R4 length :")
  print(length(tg_R1_lists))
  
  
}

#-------------------------------


print("test ictd_round1!")


print(list.files())
print('current dir:')
print(getwd())
print("what in the input folder:")
print(list.files('input/'))
input_df <- read.csv('input/input.csv')

print('read <input.csv> done!')

# Extract the names of each dataset
dataset_name <- as.character(input_df$dataset.name)

# Extract the names of the expression files that use gene name
expression_files <- as.character(input_df$native.expr.file)


output_all_ds <- c()
for(i in 1:length(expression_files))
{
  ff_tmp <- paste('input/', expression_files[i],sep='')
  print(ff_tmp)
  data_tmp <- read.csv(ff_tmp)
  rownames(data_tmp) <- data_tmp[,1]
  data_tmp <- data_tmp[,-1]
  data_bulk <- data_tmp
  ictd_result <- ICTD_round1(data_tmp)
  
  
  # fake output!!!
  ictd_prop <- data_tmp[1:12,]
  dn_tmp <- unlist(strsplit(expression_files[i],split='.',fixed=T))[[1]]
  output_tmp <- ictd_2_output(ictd_prop, dn_tmp)
  #combine prediction into big dataframe
  output_all_ds <- rbind(output_all_ds, output_tmp)
  
}

# Create the directory the output will go into
dir.create("output")

# Write the result into output directory
readr::write_csv(output_all_ds, "output/predictions.csv")
print("files exist in the output folder:")
print(list.files("output"))
print("output file dim:")
print(dim(output_all_ds))
print("output :::")
print(output_all_ds)



