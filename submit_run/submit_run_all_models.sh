#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=l40s
#SBATCH --gpus=1
#SBATCH --time=24:00:00
#SBATCH --job-name=test_short

source /home/tatzber/miniconda3/bin/activate syntheseus-full

bash /home/tatzber/syntheseus_benchmark/run_all_models.sh \
  /home/tatzber/syntheseus_benchmark/test_run.yml \
  /home/tatzber/syntheseus_benchmark/pesticide_all_models_centroids
