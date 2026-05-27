#!/usr/bin/env bash
set -o xtrace
set -e
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path $1 \
  --output-path $2
