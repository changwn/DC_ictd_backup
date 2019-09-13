## Start from this Docker image
FROM ubuntu:16.04

#FROM rocker/r-ver:3.5.1

FROM zamora/r-devtools

## Install R packages in Docker image
##RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("plyr")'
##RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('readr')"
#RUN Rscript -e "sessionInfo()"
RUN Rscript -e "library('devtools')"
RUN Rscript -e "devtools::install_github('GfellerLab/EPIC')"

#install dependent pkg
#RUN Rscript -e 'print(sessionInfo() )'
RUN Rscript -e 'install.packages("BiocManager") '
RUN Rscript -e 'BiocManager::install("impute") '
RUN Rscript -e 'BiocManager::install("GO.db")'
RUN Rscript -e 'BiocManager::install("sva")'
RUN Rscript -e 'BiocManager::install("preprocessCore")'
RUN Rscript -e 'install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE) '


#install ICTD
#ARG CACHEBUST=1	#don't use old image
RUN Rscript -e 'devtools::install_github("zy26/ICTD") '


## Copy your files into Docker images
##COPY my_model.R /usr/local/bin/
##COPY model.rds /usr/local/bin/
##RUN chmod a+x /usr/local/bin/my_model.R
#COPY ictd_model.R /usr/local/bin/
#RUN chmod a+x /usr/local/bin/ictd_model.R
COPY TCGA_IM_core_markers.RData ./
COPY ictd_model.R ./
RUN chmod a+x ictd_model.R



## Make Docker container executable
#ENTRYPOINT ["Rscript", "/usr/local/bin/my_model.R"]
#ENTRYPOINT ["Rscript", "/usr/local/bin/ictd_model.R"]
ENTRYPOINT ["Rscript", "ictd_model.R"]

