QC commands for S2A 1000G demuxlet sweep
========================================

All commands run from `/home/pr422/RDS/live/Users/Parisa`.

1. Singlet/doublet counts per geno-error setting

```
DIR=/home/pr422/RDS/live/Users/Parisa/demux_manual/S2A/1000G_test/demuxlet_GT_1000G
for f in demuxlet_err0.01.best demuxlet_err0.02.best demuxlet_err0.03.best demuxlet_err0.05.best; do
  echo "=== $f ==="
  awk 'NR>1{c[$5]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' "$DIR/$f" | sort
  echo
  echo "Donor counts (SNG only):"
  awk -F'\t' 'NR>1 && $5=="SNG"{split($6,a,","); s[a[1]]++} END{for(k in s) printf "%s\t%d\n",k,s[k]}' "$DIR/$f" | sort
  echo
done
```

2. Median NUM.SNPS / NUM.READS among singlets

```
python3 - <<'PY'
import csv,statistics
DIR='/home/pr422/RDS/live/Users/Parisa/demux_manual/S2A/1000G_test/demuxlet_GT_1000G'
files=['demuxlet_err0.01.best','demuxlet_err0.02.best','demuxlet_err0.03.best','demuxlet_err0.05.best']
for fname in files:
    ns=[]; rd=[]
    with open(f"{DIR}/{fname}") as fh:
        rdr=csv.DictReader(fh, delimiter='\t')
        for row in rdr:
            if row['DROPLET.TYPE']=='SNG':
                ns.append(int(row['NUM.SNPS']))
                rd.append(int(row['NUM.READS']))
    med_ns = statistics.median(ns) if ns else 0
    med_rd = statistics.median(rd) if rd else 0
    print(f"{fname}: SNG={len(ns)}, median NUM.SNPS={med_ns}, median NUM.READS={med_rd}")
PY
```
