#$ -q all.q
#$ -j y                      # Merge stdout and stderr
#$ -cwd                      # Use the current working directory
#$ -S /bin/bash              # Use Bash shell
#$ -terse                    # Output only job ID
#$ -pe smp 8                # Request 24 cores
#$ -N M2M88       # Name the job
#$ -l h_vmem=150G              # Request 8GB RAM
#$ -l h_rt=24:00:00          # Request 24 hour runtime
#$ -m beas                   # Email at the beginning and end of the job
#$ -M june.yoon@nationwidechildrens.org
#$ -o /igm/home/hxy008/SoloTE/M2M88.log           # Save combined stdout/stderr to this log file


module load bedtools_2.31.0
module load samtools_1.19

python3 /igm/home/hxy008/SoloTE/SoloTE_pipeline.py \
  --bam /igm/home/afh005/Bedrosian/Collaborations/switzerland_071025/P2_250303G/P2_250303G/outs/per_sample_outs/M_2M_88/count/sample_alignments.bam \
  --teannotation /igm/home/hxy008/SoloTE/te7.txt \
  --outputprefix M2M88 \
  --outputdir /igm/home/hxy008/SoloTE/M2M88 \
  -t 4
