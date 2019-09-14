

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
# ictd_2_output <- function(ictd_prop, dataset.name)
# {
#   #2
#   col_name <- matrix(colnames(ictd_prop),1,length(colnames(ictd_prop)))
#   vv <- c()
#   for(i in 1:length(col_name))
#   {
#     vv <- c(vv, rep(col_name[i], 8))
#   }
#   sample_col <- vv
#   
#   #3  
#   celltype_col <- rep(c('CD4.T.cells','CD8.T.cells','NK.cells','B.cells','monocytic.lineage','neutrophils','endothelial.cells','fibroblasts'), ncol(ictd_prop))
#   #4
#   prediction_col <- rep(runif(8,min=0.01,max=0.8), ncol(ictd_prop))  
#   
#   #bind
#   output_df <- data.frame(dataset.name=dataset.name,sample.id=sample_col,cell.type=celltype_col,prediction=prediction_col)
#   
#   return(output_df)
# }

ictd_2_output_real <- function(ictd_prop, dataset.name)
{
  #2
  col_name <- matrix(colnames(ictd_prop),1,length(colnames(ictd_prop)))
  vv <- c()
  for(i in 1:length(col_name))
  {
    vv <- c(vv, rep(col_name[i], 8))  #8 cell types
  }
  sample_col <- vv

  #3
  celltype_col <- rep(c('CD4.T.cells','CD8.T.cells','NK.cells','B.cells','monocytic.lineage','neutrophils','endothelial.cells','fibroblasts'), ncol(ictd_prop))
  #4
  #prediction_col <- rep(runif(8,min=0.01,max=0.8), ncol(ictd_prop))
  prediction_col <- c()
  cell_type_order <- c('CD4.T.cells','CD8.T.cells','NK.cells','B.cells','monocytic.lineage','neutrophils','endothelial.cells','fibroblasts')
  for(i in 1:length(col_name))
  {
    prop_tmp <- ictd_prop[cell_type_order,i]
    prediction_col <- c(prediction_col, prop_tmp)
  }


  #bind
  output_df <- data.frame(dataset.name=dataset.name,sample.id=sample_col,cell.type=celltype_col,prediction=prediction_col)

  return(output_df)
}


#-------------------------------


print("test load RData!!!")


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
  #ictd_result <- ICTD(data_tmp)
  #ictd_prop <- data_tmp[1:12,]
  ictd_prop <- matrix(runif(8*ncol(data_tmp), 0.1,0.4), 8, ncol(data_tmp))
  colnames(ictd_prop) <- paste("Setsample", 1:ncol(ictd_prop),sep = "")
  rownames(ictd_prop) <- c('CD4.T.cells','CD8.T.cells','NK.cells','B.cells','monocytic.lineage','neutrophils','endothelial.cells','fibroblasts')
  
  dn_tmp <- unlist(strsplit(expression_files[i],split='.',fixed=T))[[1]]
  
  #output_tmp <- ictd_2_output(ictd_prop, dn_tmp)
  output_tmp <- ictd_2_output_real(ictd_prop, dn_tmp)
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



