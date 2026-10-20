###expr.csv is the cell * gene expression
###chaneg expr to loom save name expr.csv
####we need create the envs of pyscenic : https://github.com/aertslab/pySCENIC
conda activate pyscenic
python /00data_Processing/1.7change.py

dir=/DATA/User/hemingmin/SCENIC_DB/mouse #
tfs=$dir/allTFs_mm.txt
feather=$dir/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v9-nr.mgi-m0.001-o0.0.tbl
echo $dir
echo $tfs
echo $feather
echo $tbl
# ensure the data 
input_loom=./sample.loom
ls $tfs  $feather  $tbl 

##grn
pyscenic grn --num_workers 20 --output adj.sample.tsv --method grnboost2 sample.loom $tfs

##cistarget
pyscenic ctx adj.sample.tsv $feather --annotations_fname $tbl --expression_mtx_fname $input_loom  --mode "dask_multiprocessing" --output reg.csv --num_workers 20  --mask_dropouts

##AUCell
pyscenic aucell $input_loom reg.csv --output out_SCENIC.loom --num_workers 3

