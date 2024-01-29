#!/bin/bash
#SBATCH -J cellpose_gpu
#SBATCH -A gts-wratcliff3
#SBATCH -q inferno
#SBATCH -N 1 --gres=gpu:V100:1
#SBATCH --mem-per-gpu=16G
#SBATCH -t 24:00:00
#SBATCH -o cellpose_gpu.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ktong34@gatech.edu
    
# Run cellpose segmentation given directory containing input images
source activate cellpose_gpu_v2.2.2
python -m cellpose \
--dir $1 \
--pretrained_model $2 --chan 0 --chan2 0 --diameter $3 \
--save_png --no_npy \
--verbose --use_gpu

