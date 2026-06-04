#!/usr/bin/env bash

set -o xtrace
set -e

source activate clara_report

QMD=$(realpath "$1")
OTU=$(realpath "$2")
TAXONOMY=$(realpath "$3")
METADATA=$(realpath "$4")

REPORT_NAME="$5"

WORKDIR=$(dirname "${OTU}")

pushd "${WORKDIR}" >/dev/null

TMP_PARAMS="${WORKDIR}/${REPORT_NAME}_quarto_params.yml"

cat > "${TMP_PARAMS}" << EOF
otu: "${OTU}"
taxonomy: "${TAXONOMY}"
metadata: "${METADATA}"
EOF

echo "PWD=$(pwd)"
echo "QMD=${QMD}"

ls -lh "${QMD}"
ls -lh "${TMP_PARAMS}"

quarto render "${QMD}"  \
--execute-params "${TMP_PARAMS}" \
--output "${REPORT_NAME}" 

rm -f "${TMP_PARAMS}"

popd >/dev/null

source deactivate


