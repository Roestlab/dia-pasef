#!/bin/bash
#SBATCH --time=07:00:00
#SBATCH --account=def-hroest
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu 1024M
#SBATCH --array=4-27

# Python venv
source /project/6011811/bin/pyenv_27/bin/activate

# Setup environment
source_dir=/project/6011811/frankmak/code/tims/src/
bruker_sdk=/project/6011811/frankmak/code/tdf-sdk-2.3.0/linux64/libtimsdata.so
python_bindings=/project/6011811/frankmak/code/tdf-sdk-2.3.0/linux64/timsdata.py
ln -sf  $bruker_sdk libtimsdata.so # Caution: overwrites existing symlinks
ln -sf  $python_bindings timsdata.py # Caution: overwrites existing symlinks
cp ${source_dir}/convert_single.py . # Use most current version of the conversion script

# Declare filepaths
mzML_dir=/scratch/frankmak/tims/dia/20180327_firstRuns_data/mzML
# This can be run with a single file
# analysis_dir=20180320_AnBr_SA_diaPASEF_200ng_HeLa_Rost_Method_4_b_1800V_01_A1_01_2151.d
# outfile=$(echo $analysis_dir | sed 's/\(.*\.\)d/\1mzML/' )
# Or multiple files
analysis_dir=$( ls | grep ".d" | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' )
outfile=$(echo $analysis_dir | sed 's/\(.*\.\)d/\1mzML/' )

# Convert to mzML
mkdir -p $mzML_dir
echo "Converting $analysis_dir to $outfile ..."
python convert_single.py $analysis_dir ${mzML_dir}/$outfile
echo "Successfully completed"
echo "Storing output in $analysis_dir  ${mzML_dir}/$outfile"


