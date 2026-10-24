# CLARA: A Containerized Lightweight Automated Reproducible Analysis for 16S rRNA Amplicon Sequencing

In short:

```bash
git clone https://github.com/andvides/clara16s.git
cd clara16s
git checkout tags/v1.0.1

export CLARA16s_ROOT=$(pwd)
export PATH=$CLARA16s_ROOT/workflow/scripts:$PATH

conda install -n base -c conda-forge mamba
conda create -n clara16s
conda env update -n clara16s --file $CLARA16s_ROOT/workflow/envs/clara16s.yaml
conda activate clara16s


clara16s --cores 8 > out 2>&1 & # running locally on Docker
clara16s --cores 8 --use-singularity # running locally on Singularity

clara16s_web_docker​
```

