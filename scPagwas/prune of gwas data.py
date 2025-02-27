####----- step1: Extract 503 EUR-Sample's 1000G.EUR.QC.chr------######
### 1.1 1kg-plink.txt prepare
cd ./scPagwas/gwas_data_process
ssh s001
 
ls ./reference/1000Genomes/raw_vcf/plink/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.bim|cut -f1-6 -d '.' |head
ls ./reference/1000Genomes/raw_vcf/plink/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.bim|cut -f1-6 -d '.' > ./scPagwas/gwas_data_process/aa.txt
 
ls ./reference/1000Genomes/raw_vcf/plink/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.bim|cut -f2 -d '.' |head
ls ./reference/1000Genomes/raw_vcf/plink/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.bim|cut -f2 -d '.' > ./scPagwas/gwas_data_process/bb.txt
 
# paste  -d\\t aa.txt bb.txt > 1kg-plink.txt 
paste aa.txt bb.txt  > 1kg-plink.txt  
 
 
### 1.1 for-subset-sample.txt prepare
awk '$3=="EUR" {print $1,$1}' integrated_call_samples_v3.20130502.ALL.panel |head
awk '$3=="EUR" {print $1,$1}' integrated_call_samples_v3.20130502.ALL.panel > for-subset-sample.txt

### 1.2 plink---get_extract_sample.py
cat get_extract_sample.py
import sys
import os
 
if __name__ == '__main__':
    keep_sample_file = sys.argv[1]
    plink_list_file = sys.argv[2]
    
    with open(plink_list_file) as f:
        for line in f:
            b_file, chrome, = line.split()
            cmd = rf'''/data/home/lingw/bin/plink --allow-no-sex --bfile {b_file} \
                --keep {keep_sample_file} --make-bed --out 1000G.EUR.QC.{chrome}'''
            print(cmd)
            print('#' * 20)
            print()
            
python get_extract_sample.py for-subset-sample.txt 1kg-plink.txt > run.sh
nohup sh run.sh  
## We get 1000G.EUR.QC.chr*

# step2: Extract 1000G.EUR.QC.chr-sub-chr of SNP in my-MDD
### 2.1 sub-sample-plink-list.txt prepare
ls 1000G.EUR.QC.chr*.bim|cut -f1-4 -d '.' > aa.txt
ls 1000G.EUR.QC.chr*.bim|cut -f4 -d '.' > bb.txt
#paste  -d\\t aa.txt bb.txt 
paste  -d aa.txt bb.txt  > sub-sample-plink-list.txt
 
### 2.1 extract_snp.txt prepare
awk -F' ' '{print $3}' MDD_gwas_data.txt > extract_snp.txt 
touch 
ll -rt 

### 2.2 plink---get_extract_snp.py
cat > get_extract_snp.py
import sys
import os
 
if __name__ == '__main__':
    keep_snp_file = sys.argv[1]
    plink_list_file = sys.argv[2]
    
    with open(plink_list_file) as f:
        for line in f:
            b_file, chrome, = line.split()
            cmd = rf'''/data/home/lingw/bin/plink --allow-no-sex --bfile {b_file} \
                --extract {keep_snp_file} --make-bed --out 1000G.EUR.QC.{chrome}-sub-SNP'''
            print(cmd)
            print('#' * 20)
            print()


python get_extract_snp.py extract_snp.txt sub-sample-plink-list.txt > run-extract-snp.sh
nohup sh ./run-extract-snp.sh &  
## We get 1000G.EUR.QC.chr-sub-chr*

#check data
wc -l extract_snp.txt  #8481297 extract_snp.txt
wc -l 1000G.EUR.QC.chr*-sub-SNP.bim  #7829695 total

# step3: Pruning and get MDD.chr*_plink_prune_EUR_filtered_LD0.8.prune.in
### 3.1 sub-sample-SNP-list.txt  prepare
ls 1000G.EUR.QC.chr*-sub-SNP.bim|cut -f1-3 -d '-' > aa.txt
ls 1000G.EUR.QC.chr*-sub-SNP.bim|cut -f4 -d '.' |cut -f1 -d '-' >bb.txt
paste aa.txt bb.txt  > sub-sample-SNP-list.txt

### 3.2 plink----get_pruning.py
cat > get_pruning.py
import sys
import os
 
if __name__ == '__main__':
    # keep_snp_file = sys.argv[1]
    plink_list_file = sys.argv[1]
    
    with open(plink_list_file) as f:
        for line in f:
            b_file, chrome, = line.split()
            cmd = rf'''/data/home/lingw/bin/plink --allow-no-sex --bfile {b_file} \
                --indep-pairwise 50 5 0.8  --out MDD.{chrome}_plink_prune_EUR_filtered_LD0.8'''
            print(cmd)
            print('#' * 20)
            print()
            
python get_pruning.py sub-sample-SNP-list.txt > run-pruning.sh
nohup sh run-pruning.sh
jobs

### 3.3 Combine all chr and get MDD_EUR_LD0.8.prune
cat MDD.chr*_plink_prune_EUR_filtered_LD0.8.prune.in > MDD_EUR_LD0.8.prune
 
head MDD_EUR_LD0.8.prune
wc -l MDD_EUR_LD0.8.prune  #2716381
