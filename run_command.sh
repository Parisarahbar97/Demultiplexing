cd /home/pr422/RDS/live/Users/Parisa/Demultiplexing

# Single sample (e.g. D10A) from your CSV:
nextflow run main.nf -profile dsi \
  --samples   "/home/pr422/RDS/live/Users/Parisa/Demultiplexing/examples/samples.csv" \
  --outdir    "/home/pr422/RDS/live/Users/Parisa/demux_out_nf/D10A" \
  -with-report -with-timeline -with-trace