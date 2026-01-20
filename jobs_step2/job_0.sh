#!/bin/bash
#SBATCH --time=144:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=256GB
#SBATCH --job-name=immuno_retrovirus
#SBATCH --account=PDE0005
#SBATCH -e immuno_retrovirus-%j.err
#SBATCH -o immuno_retrovirus-%j.out

set echo

# slurm starts job in working DIR
cd $SLURM_SUBMIT_DIR

# set up software environment
module load intel/2021.10.0
module load miniconda3/24.1.2-py310

source activate immuno_telescope

export CMPLR_ROOT="/apps/spack/0.21/ascend/linux-rhel9-zen2/intel-oneapi-compilers/gcc/11.4.1/2023.2.3-nkzjnam/compiler/2023.2.3"
export INTEL64="$CMPLR_ROOT/linux/compiler/lib/intel64_lin"
export LD_LIBRARY_PATH="$INTEL64:$CMPLR_ROOT/linux/lib:$CMPLR_ROOT/linux/lib/x64:$CONDA_PREFIX/lib:${LD_LIBRARY_PATH}"
# Optional (only if you still see unresolved Intel symbols):
export LD_PRELOAD="$INTEL64/libimf.so:$INTEL64/libsvml.so:$INTEL64/libirc.so"

# 4) (Optional) sanity check that the Telescope C extension loads
python - <<'PY'
import ctypes, glob, os, sys
pats = glob.glob(os.path.expanduser("~/.conda/envs/immuno_telescope/lib/python3.6/site-packages/telescope/utils/calignment*.so"))
if not pats:
    sys.exit("Telescope extension not found")
ctypes.CDLL(pats[0])
print("Loaded Telescope extension OK:", pats[0])
PY

export PROJECT_PATH="/fs/scratch/PDE0005/projects/telescope_TCC00573"
export SAMPLES_PATH="/fs/ess/PDE0054/TCC00573-ROLF/RNAseq/fastq"

# copy data to compute node local space $TMPDIR
cp -R     /fs/scratch/PDE0005/telescope_lungsamples/code/lib     /fs/scratch/PDE0005/telescope_lungsamples/code/run.py     $TMPDIR/
    
cp -R     /fs/scratch/PDE0005/telescope_lungsamples/data/Indexes     /fs/scratch/PDE0005/telescope_lungsamples/data/REF     $TMPDIR/

cd $TMPDIR

mkdir output

# $PFSDIR: /fs/scratch/PAS0438/osu9787/18091549.owens

    
python run.py     --gtf REF/gencode.v39.annotation.gtf     --genome REF/GRCh38.p13.genome.fa     --transcript REF/gencode.v39.transcripts.fa     --herv_gtf REF/HG38_HERV_LINE_all_families_telescope_ann.gtf     --bowtiew2_idx Indexes/gencode.v39_bowtie2/human     --workflows TELESCOPE     --out_dir output     --scratch_dir $PROJECT_PATH/tmp_0     --n_cores 24     --trimgalore_n_cores 6     --telescope_n_cores 1     --seed 123456     --fastq_mode


LOGS="$TMPDIR/output/TELESCOPE/logs.tsv"
if [ -f "$LOGS" ]; then
  cp $LOGS $PROJECT_PATH/output/TELESCOPE/logs_0.tsv
fi

BOWTIE_DIR="$TMPDIR/output/TELESCOPE"
if [ -d "$BOWTIE_DIR" ]; then
  cp -R $BOWTIE_DIR/* $PROJECT_PATH/output/TELESCOPE/
fi
