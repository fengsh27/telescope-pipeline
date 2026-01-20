#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=256G
#SBATCH --job-name=immuno_retrovirus
#SBATCH --account=PDE0005
#SBATCH --error=immuno_retrovirus-%j.err
#SBATCH --output=immuno_retrovirus-%j.out

set -euo pipefail
set -x

cd "${SLURM_SUBMIT_DIR}"

module load miniconda3/24.1.2-py310
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate immuno_telescope

export PROJECT_PATH="/fs/scratch/PDE0005/projects/telescope_TCC00573"
export SAMPLES_PATH="/fs/ess/PDE0054/TCC00573-ROLF/RNAseq/fastq"

rm -rf "${PROJECT_PATH}/tmp_0"
mkdir -p "${PROJECT_PATH}/tmp_0"

# copy code + manifests into node-local space
cp -R   "/fs/scratch/PDE0005/telescope_lungsamples/code/lib"   "/fs/scratch/PDE0005/telescope_lungsamples/code/run.py"   "${PROJECT_PATH}/manifests/manifest_0.tsv"   "${TMPDIR}/"

cp -R   "/fs/scratch/PDE0005/telescope_lungsamples/data/Indexes"   "/fs/scratch/PDE0005/telescope_lungsamples/data/REF"   "${TMPDIR}/"

cd "${TMPDIR}"
mkdir -p output

python run.py   --manifest manifest_0.tsv   --gtf REF/gencode.v39.annotation.gtf   --genome REF/GRCh38.p13.genome.fa   --transcript REF/gencode.v39.transcripts.fa   --herv_gtf REF/HG38_HERV_LINE_all_families_telescope_ann.gtf   --bowtiew2_idx Indexes/gencode.v39_bowtie2/human   --workflows BOWTIE   --out_dir output   --samples_dir "${SAMPLES_PATH}"   --n_cores "${SLURM_CPUS_PER_TASK}"   --trimgalore_n_cores 6   --telescope_n_cores 1   --seed 123456   --scratch_dir "${PROJECT_PATH}/tmp_0"   --fastq_mode   --sample_check_n_cores 6

LOGS="${TMPDIR}/output/BOWTIE/logs.tsv"
if [ -f "${LOGS}" ]; then
  cp "${LOGS}" "${PROJECT_PATH}/output/BOWTIE/logs_0.tsv"
fi

BOWTIE_DIR="${TMPDIR}/output/BOWTIE"
if [ -d "${BOWTIE_DIR}" ]; then
  rsync -a "${BOWTIE_DIR}/" "${PROJECT_PATH}/output/BOWTIE/"
fi

