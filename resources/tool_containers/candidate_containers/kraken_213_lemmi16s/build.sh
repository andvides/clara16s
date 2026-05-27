#docker build --rm -t ezlab/kraken2_213_lemmi16s:v1.0_cv1 .
#docker login quay.io -u ezlab
#docker tag ezlab/kraken2_213_lemmi16s:v1.0_cv1 quay.io/ezlab/kraken2_213_lemmi16s:v1.0_cv1

#To build the local image
docker build --rm . -t $USER/kraken2_16s
#docker run --rm -it -v $(pwd):/app $USER/kraken2_16s

