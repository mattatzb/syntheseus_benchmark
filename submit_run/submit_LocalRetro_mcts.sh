#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=l40s
#SBATCH --gpus=1
#SBATCH --time=24:00:00
#SBATCH --job-name=LocalRetro_mcts

source /home/tatzber/miniconda3/bin/activate syntheseus-full

syntheseus search \
  --config /home/tatzber/syntheseus_benchmark/results_6_runs/configs/LocalRetro_mcts.yml
