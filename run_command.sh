cd /home/pr422/RDS/live/Users/Parisa/Demultiplexing

# (Optional) test only D2A row in examples/samples.csv first
nextflow run main.nf -profile docker \
  --manifest /home/pr422/RDS/live/Users/Parisa/Demultiplexing/examples/samples.csv \
  --outdir   /home/pr422/RDS/live/Users/Parisa/demux_out
