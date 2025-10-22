process qc_demuxlet {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/${sample_id}_demuxlet", mode: 'copy'

  input:
  tuple val(sample_id), path(best), path(barcodes), path(sm_list)

  output:
  path("QC_summary.txt")
  path("donor_QC.tsv")

  script:
  """
  set -euo pipefail
  BEST="$best"
  BARCODES="$barcodes"
  SM_LIST="$sm_list"

  BC_TMP=$(mktemp)
  case "$BARCODES" in
    *.gz)  gzip -cd "$BARCODES" > "$BC_TMP.raw" ;;
    *)     cat  "$BARCODES" > "$BC_TMP.raw" ;;
  esac
  first_line=$(head -n1 "$BC_TMP.raw")
  if echo "$first_line" | grep -qi 'barcode'; then
    awk -F'[,\t]' 'NR>1{print $1}' "$BC_TMP.raw" > "$BC_TMP"
  else
    mv "$BC_TMP.raw" "$BC_TMP"
  fi

  {
    echo -e "donor\tsinglets\tdoublets_as_d1\tdoublets_as_d2\tdoublets_total"
    tail -n +2 "$BEST" | awk -F'\\t' '
      \$5=="SNG"{split(\$6,a,\",\"); s[a[1]]++; seen[a[1]]=1}
      \$5=="DBL"{split(\$6,a,\",\"); d1[a[1]]++; d2[a[2]]++; seen[a[1]]=1; seen[a[2]]=1}
      END{for(d in seen){sd=s[d]+0; x=d1[d]+0; y=d2[d]+0; printf "%s\\t%d\\t%d\\t%d\\t%d\\n", d, sd, x, y, x+y}}' \
    | sort -k1,1
  } > donor_QC.tsv

  {
    echo "rows_vs_barcodes"
    echo -n "barcodes "; wc -l < "$BC_TMP"
    echo -n "best_rows "; awk 'END{print NR-1}' "$BEST"
    echo -n "overlap   "; comm -12 <(tail -n +2 "$BEST" | cut -f2 | sort) <(sort "$BC_TMP") | wc -l
    echo
    echo "status_counts"
    tail -n +2 "$BEST" | awk -F'\\t' '{c[\$5]++} END{for(k in c) printf "%s\\t%d\\n",k,c[k]}' | sort
    echo
    echo "sng_posteriors"
    tail -n +2 "$BEST" | awk -F'\\t' '\$5=="SNG"{t++; if(\$12>=0.99) p99++; else if(\$12>=0.9) p90++}
      END{printf "SNG=%d  >=0.99=%d (%.1f%%)  0.9-<0.99=%d (%.1f%%)\\n", t,p99,100*p99/t,p90,100*p90/t+0}'
    echo
    echo "top_dbl_pairs"
    tail -n +2 "$BEST" | awk -F'\\t' '\$5=="DBL"{split(\$6,a,\",\"); d1=a[1]; d2=a[2]; pair=(d1<d2?d1\"+\"d2:d2\"+\"d1); p[pair]++} END{for(k in p) print p[k],k}' | sort -k1,1nr | head
    echo
    echo "donors_seen_vs_expected"
    tail -n +2 "$BEST" | awk -F'\\t' '\$5=="SNG"{print \$13} \$5=="DBL"{split(\$6,a,\",\"); print a[1]; print a[2]}' | sort -u > donors_found.txt
    echo "expected_but_NOT_seen:"; comm -13 donors_found.txt <(sort "$SM_LIST") || true
    echo "seen_but_NOT_expected:"; comm -23 donors_found.txt <(sort "$SM_LIST") || true
  } > QC_summary.txt

  rm -f "$BC_TMP" "$BC_TMP.raw" donors_found.txt
  """
}
