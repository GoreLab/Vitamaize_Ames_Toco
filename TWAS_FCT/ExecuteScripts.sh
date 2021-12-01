cd /workdir/rt475/TWAS_2020-07
mkdir LOGFILE

# Merge overlapped genes: do this before FC test
cd /workdir/rt475/TWAS_2020-07/RAWDATA/Annotation
nohup R --vanilla --slave < 1-MergeGenes_B73.R > log_1.txt

# Run 1.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.1-MakeDataset.R --args v1.1 B73 > LOGFILE/log_1.1_v1.1_B73.txt

# Run 1.2
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.2-AttachGDD.R --args v1.1 B73 > LOGFILE/log_1.2_v1.1_B73.txt

# Run 1.3
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.3-GeneFiltering.R --args v1.1 B73 0.5 > LOGFILE/log_1.3_v1.1_B73.txt

# Run 999.001
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 999.001-Diagnosis_for_gene_filtering.R --args B73 > LOGFILE/log_999.001_B73.txt

# Run 999.002
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 999.002-Diagnosis_for_Outlier_Removal.R --args v1.1 B73 > LOGFILE/log_999.002_v1.1_B73.txt

# Run 999.003
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 999.003-Diagnosis_for_Outlier_Removal_2.R --args v1.1 B73 > LOGFILE/log_999.003_v1.1_B73.txt

# Run 1.4
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.4-OutlierRemoval_and_Imputation.R --args v1.1 B73 100 > LOGFILE/log_1.4_v1.1_B73_100.txt

# Run 1.5.1 -- this part is done by other shell scripts
cd /workdir/rt475/TWAS_2020-07
nohup sh 1.5-RunMultipleR_v1.1_B73.sh > LOGFILE/log_1.5_v1.1_B73.txt &
# You must wait until the BLUE calculation ends

# Run 1.5.2
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.5.2-BLUE_merge.R --args v1.1 B73 > LOGFILE/log_1.5.2_v1.1_B73.txt

# Run 1.5.3
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.5.3-BLUE_RmErr.R --args v1.1 B73 > LOGFILE/log_1.5.3_v1.1_B73.txt

# Run 1.5.4
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.5.4-BLUE_SampleFiltering.R --args v1.1 B73 > LOGFILE/log_1.5.4_v1.1_B73.txt

# Run 2.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 2.1-Peer_Use25Fact.R --args v1.1 B73 > LOGFILE/log_2.1_v1.1_B73.txt

# Run 2.2
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 2.2-MakeScreePlot.R --args v1.1 B73 > LOGFILE/log_2.2_v1.1_B73.txt

# Run 2.3
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 2.3-Peer_UseOptFact.R --args v1.1 B73 11 > LOGFILE/log_2.3_v1.1_B73.txt

# Run 3.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 3.1-OutlierRemoval.R --args v1.1 B73 11 > LOGFILE/log_3.1_v1.1_B73.txt

# Run 4.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 4.1-BIC_for_TWAS.R > LOGFILE/log_4.1.txt

# Run 4.2
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 4.2-TWAS.R --args v1.1 B73 > LOGFILE/log_4.2_v1.1_B73.txt

# Run 5.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 5.1-CombinedTest.R --args v1.1 B73 > LOGFILE/log_5.1_v1.1_B73.txt

# Run VTE7
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 999.999.001-AnalyseVte7_01.R > LOGFILE/log_vte7_01.txt
nohup R --vanilla --slave < 999.999.002-AnalyseVte7_02.R > LOGFILE/log_vte7_02.txt
