#docker build --rm -t ezlab/qiime2_20228_lemmi16s:v1.0_cv1 .
#docker login quay.io -u ezlab
#docker tag ezlab/qiime2_20228_lemmi16s:v1.0_cv1 quay.io/ezlab/qiime2_20228_lemmi16s:v1.0_cv1

#To build the local image
docker build --rm . -t $USER/qiime2
#docker run --rm -it -v $(pwd):/app $USER/qiime2
