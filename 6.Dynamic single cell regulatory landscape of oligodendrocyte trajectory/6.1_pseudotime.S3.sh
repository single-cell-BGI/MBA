# run "Rscript /hwfssz1/ST_PRECISION/USER/zhongyu/2.Scripts/PseudotemporalAnalysis/Monocle_v1.1/Monocle_pseudotime.R -h" for help
/hwfssz1/ST_PRECISION/PUB/softwares/R-3.5.1/bin/Rscript 3.5.1Monocle_pseudotime.R \
	--exp /hwfssz1/ST_PRECISION/USER/lizihao/download/OPC-OLI/C4/OPC_OLI_matrix.txt \
	--cell_sheet /hwfssz1/ST_PRECISION/USER/lizihao/download/OPC-OLI/C4/info.list \
	--gene_anno /hwfssz1/ST_PRECISION/USER/lizihao/download/OPC-OLI/C4/OPC-OLI.pseudotiomemetadata \
        --data_type UMI \
	--cell_group seurat_clusters \
	--batch NULL \
	--min_expr 0.5 \
	--min_pct 0.05 \
	--qval_threshold 1e-10 \
	--mean_cutoff 1 \
	--dispersion_cutoff 1 \
	--thread 4 \
	--reverse FALSE \
	--n_cluster 6 \
	--out /hwfssz1/ST_PRECISION/USER/lizihao/download/OPC-OLI/C4/out/c4

