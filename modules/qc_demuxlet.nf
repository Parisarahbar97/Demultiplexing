process QC {
  tag "${sample_id}"
  label 'QC'

  input:
  tuple val(sample_id), path(best), path(barcodes), path(sm_list)

  output:
  path("${sample_id}.demuxlet.qc.txt")

  script:
  """
  set -euo pipefail

  N_BEST=\$(( \$(wc -l < "${best}") - 1 ))
  N_10X=\$( ( [[ "${barcodes}" == *.gz ]] && gzip -cd "${barcodes}" || cat "${barcodes}" ) | wc -l )
  OVERLAP=\$(comm -12 <(tail -n +2 "${best}" | cut -f2 | sort) <( ( [[ "${barcodes}" == *.gz ]] && gzip -cd "${barcodes}" || cat "${barcodes}" ) | sort ) | wc -l)

  {
    echo "sample\t${sample_id}"
    echo "best_rows\t\$N_BEST"
    echo "barcodes_10x\t\$N_10X"
    echo "barcode_overlap\t\$OVERLAP"
    echo ""

    echo "[status_counts]"
    tail -n +2 "${best}" | awk -F'\\t' '{c[\$5]++} END{for(k in c) printf "%s\\t%d\\n",k,c[k]}' | sort

    echo ""
    echo "[snp_depth]"
    tail -n +2 "${best}" | awk -F'\\t' '{print \$3}' \\
    | awk '{s+=\$1;n++; if(\$1<20) l20++} END{printf "mean_SNPS\\t%.2f\\n<20_SNPS\\t%d (%.2f%%)\\n", s/n, l20, 100*l20/n}'

    echo ""
    echo "[posterior_on_SNG]"
    tail -n +2 "${best}" \\
    | awk -F'\\t' '\$5=="SNG"{t++; if(\$12>=0.99) p99++; else if(\$12>=0.9) p90++} END{printf "SNG\\t%d\\n>=0.99\\t%d (%.1f%%)\\n0.9-<0.99\\t%d (%.1f%%)\\n", t,p99,100*p99/t,p90,100*p90/t+0}'

    echo ""
    echo "[donors_seen_vs_expected]"
    tail -n +2 "${best}" | awk -F'\\t' '\$5=="SNG"{print \$13} \$5=="DBL"{split(\$6,a, ","); print a[1]; print a[2]}' | sort -u > donors_found.txt
    echo "expected_but_NOT_seen:"
    comm -13 <(sort donors_found.txt) <(sort "${sm_list}")
    echo "seen_but_NOT_expected:"
    comm -23 <(sort donors_found.txt) <(sort "${sm_list}")
  } > ${sample_id}.demuxlet.qc.txt
  """
}
