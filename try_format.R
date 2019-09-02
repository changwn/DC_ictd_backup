

setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_06_DC/Tumor-Deconvolution-Challenge-Workflow-master/example_files/fast_lane_dir')



library('ICTD')


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
  celltype_col <- rep(c('CD4.T.cell','CD8.T.cell','NK.cells','B.cells','monocytic.lineage','neutrophils','endothelial.cells','fibroblats'), ncol(ictd_prop))
  #4
  prediction_col <- rep(c(0.1:0.8), ncol(ictd_prop))  
  
  #bind
  output_df <- data.frame(dataset.name=dataset.name,sample.id=sample_col,cell.type=celltype_col,prediction=prediction_col)

  return(output_df)
}




#-------------------------------




print(list.files())
print(getwd())
input_df <- readr::read_csv('input.csv')

# Extract the names of each dataset
dataset_name <- input_df$dataset.name

# Extract the names of the expression files that use gene name
expression_files <- input_df$native.expr.file


output_all_ds <- c()
for(i in 1:length(expression_files))
{
  data_tmp <- read.csv(expression_files[i])
  rownames(data_tmp) <- data_tmp[,1]
  data_tmp <- data_tmp[,-1]
  #ictd_result <- ICTD(data_tmp)
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

print(list.files("output"))


