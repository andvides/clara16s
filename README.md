# CLARA16s

**A Containerized Lightweight Automated Reproducible Analysis for 16S rRNA Amplicon Sequencing**

CLARA16s is a containerized workflow designed to automate the analysis of 16S rRNA amplicon sequencing data. It integrates workflow orchestration and software containerization to support reproducible execution across computing environments.

The project aims to simplify pipeline execution, manage software dependencies, and provide a consistent framework for processing 16S rRNA sequencing datasets.

## Installation and Setup

### 1. Clone the Repository

```bash
git clone https://github.com/andvides/clara16s.git
cd clara16s
git checkout tags/v1.0.1
```

### 2. Configure Environment Variables

```bash
export CLARA16s_ROOT=$(pwd)
export PATH=$CLARA16s_ROOT/workflow/scripts:$PATH
```

### 3. Create the Conda Environment

```bash
conda install -n base -c conda-forge mamba
conda create -n clara16s
conda env update -n clara16s --file $CLARA16s_ROOT/workflow/envs/clara16s.yaml
conda activate clara16s
```

## Configuration

### 1. Edit the Main Configuration File

Edit `config/config.yaml` to configure the analysis.

### 2. Select Tools to Execute

Visit the [CLARA16s Container Registry](https://quay.io/organization/clara16s) to explore the available container images and identify the tools to execute.

Include the corresponding tool configuration files in `/data/tools/<tool>.yaml`.

## Execution

### 1. Run CLARA16s Locally

```bash
clara16s --cores 8 > out 2>&1 &
```

### 2. Run CLARA16s Using Singularity

```bash
clara16s --cores 8 --use-singularity
```

## Results Visualization

### Run CLARA16s Web Interface

```bash
clara16s_web_docker
```
