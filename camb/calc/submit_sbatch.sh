#! /bin/bash

#=====================================================
#===== Modify the following options for your job =====
#=====    DON'T remove the #! /bin/bash lines    =====
#=====      DON'T comment #SBATCH lines          =====
#=====        of partition,account and           =====
#=====                qos                        =====
#=====================================================

# Specify the partition name from which resources will be allocated  
#SBATCH --partition=ali

# Specify which expriment group you belong to.
# This is for the accounting, so if you belong to many experiments,
# write the experiment which will pay for your resource consumption
#SBATCH --account=alicpt

# Specify which qos(job queue) the job is submitted to.
#SBATCH --qos=regular


# ====================================
#SBATCH --job-name=wym

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=1000
# SBATCH --exclude=aliws007
#SBATCH --nodelist=aliws005

#SBATCH -o out.log
#SBATCH --error=err.log
#SBATCH --mail-type=END
# SBATCH --mail-user=wangyiming@ihep.ac.cn

# this setting is used when you use openmp
# if [ -n "$SLURM_CPUS_PER_TASK" ]; then
#   omp_threads=$SLURM_CPUS_PER_TASK
# else
#   omp_threads=1
# fi
# export OMP_NUM_THREADS=$omp_threads

# or use relative path
# mpiexec python -u tod_gen4cc.py
# mpirun -np 7 ./cosmomc test.ini
# python as.py
# mpiexec python -u tmp.py
# mpiexec python -u dm_fsr_7params.py
mpiexec python -u getpd.py
