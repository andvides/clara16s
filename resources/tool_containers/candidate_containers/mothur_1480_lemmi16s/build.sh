#docker build --rm -t ezlab/mothur_148_lemmi16s:v1.0_cv1 .
#docker login quay.io -u ezlab
#docker tag ezlab/mothur_148_lemmi16s:v1.0_cv1 quay.io/ezlab/mothur_148_lemmi16s:v1.0_cv1

#To build the local image
docker build --rm . -t $USER/mothur
#docker run --rm -it -v $(pwd):/app $USER/mothur
