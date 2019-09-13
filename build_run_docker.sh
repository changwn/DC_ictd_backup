# shell for docker

docker build -t ictdtest:v1.0 .

docker image ls

docker rmi -f `docker images  | grep "<none>" | awk '{print $3}'`

docker run ictdtest:v1.0



# test docker model locally
#docker run --network none -v /path/to/input_files/:/input:ro -v /path/to/output_directory/:/output:rw your_docker_image

#docker run --network none -v /home/aelamb/repos/Tumor-Deconvolution-Challenge-Workflow/example_files/fast_lane_dir/:/input:ro -v /home/aelamb/output/:/output:rw docker.synapse.org/syn11738103/td_cibersort_course

docker run --network none -v /Users/chang/Documents/dc_ictd/fast_lane_dir/input/:/input:ro -v /Users/chang/Documents/dc_ictd/fast_lane_dir/output/:/output:rw ictdtest:v1.0


#--------------------------
#
#push image to synapse
#
#-------------------------

docker login docker.synapse.org
chang91
cwn3821365

docker logout docker.synapse.org

#docker image

#push ICTD account
#syntax docker push docker.synapse.org/<Your project ID>/<Repo name>:<Tag>
docker tag ictdtest:v1.0 docker.synapse.org/syn20551227/ictdtest:v1.0
docker image ls
docker push docker.synapse.org/syn20551227/ictdtest

#view docker
https://www.synapse.org/#!Synapse:syn20551227/docker

#-----------------------
#push to test account
#-----------------------
docker login docker.synapse.org --username=hearthewind 
hearthewind
cwn3821365
docker images ls

docker tag ictdtest:v1.0 docker.synapse.org/syn20757641/test:v1.0



docker image ls
docker push docker.synapse.org/syn20757641/ictdtest

docker push docker.synapse.org/syn20757641/test:v1.0





