#!/usr/bin/env bash
set -o xtrace
set -e
ref_seq_qza=$1
ref_taxo_qza=$2
class_name=$3
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $ref_seq_qza \
  --i-reference-taxonomy $ref_taxo_qza \
  --o-classifier $class_name

