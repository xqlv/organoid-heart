meta=~/cellphonedb/D8/meta.txt
counts=~/cellphonedb/D8/counts.txt
output=~/cellphonedb/D8/output

meanspath=$output/means.txt
pvaluespath=$output/pvalues.txt

mkdir $output
/opt/cellphonedb/3.1.0/bin/cellphonedb method statistical_analysis $meta $counts \
--output-path $output \
--counts-data hgnc_symbol \
--threads 10 \
--iterations 1000 \
--threshold 0.1 \
--pvalue 0.05

/opt/cellphonedb/3.1.0/bin/cellphonedb plot dot_plot --means-path $meanspath --pvalues-path $pvaluespath --output-path $output
/opt/cellphonedb/3.1.0/bin/cellphonedb plot heatmap_plot $meta --pvalues-path $pvaluespath --output-path $output