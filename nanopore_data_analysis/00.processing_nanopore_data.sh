
### Written by Yiming Huang, Wang Lab at Columbia University, yiminghuang0405@gmail.com ###

### a demo of compressed raw nanopore reads in FASTQ format are in ./rawdata/merge.fastq.gz ###
### the information of barcodes usage for each sample pool are in ./sample.bc.list ###
### the list of samples are in ./sample.list ###
### MUSCLE need to be in the system PATH for this script ###



### enable other script used in the analysis
chmod +x ./bin/*

### make output folder ###
mkdir -p output_00_readsProcess
mkdir -p output_01_readsDemultiplex

### uncompressed the input reads ###
gunzip ./rawdata/merge.fastq.gz

### the pool to be analyzed ###
targetPool=rsgApool


### remove reads that are too short ###
./bin/cutoff_length_fq.py ./rawdata/merge.fastq ./rawdata/merge.pass.fastq 1600

### rename reads and output to FASTA format ###
./bin/convertFastq2Fasta.py ./rawdata/merge.pass.fastq \
        ./output_00_readsProcess/all.merge.fasta

### rename reads and output to FASTQ format ###
./bin/convertFastq2Fastq.py ./rawdata/merge.pass.fastq \
        ./output_00_readsProcess/all.merge.fastq

### split reads into 96 batches so we can run the barcode identification in parallel
./bin/splitFasta.py ./output_00_readsProcess/all.merge.fasta \
	./output_00_readsProcess/merge.split 96

### run barocde identification based on archor sequences at both end using 96 threads
for i in `seq 0 95`
do
    ./bin/statFastaQC_multiple.py ./output_00_readsProcess/merge.split.$i\.fasta \
            ./sample.bc.tsv ./output_00_readsProcess/$targetPool.metadata.split.$i\.tsv \
            $targetPool /dev/shm/$targetPool\.batch_$i &
done
wait

### merge the barcode identification result from 96 threads
./bin/mergeMetadata.py ./output_00_readsProcess/$targetPool.metadata.split \
        ./output_00_readsProcess/$targetPool.metadata.tsv 96


### filter out reads with poor quality, >3bp single-base difference or >1bp gap
./bin/statFastaBarcode_multiple.py output_00_readsProcess/$targetPool.metadata.tsv \
        ./sample.bc.tsv \
        ./output_00_readsProcess/all.merge.fastq \
        ./output_00_readsProcess/$targetPool.metadata.pass.tsv \
        ./output_00_readsProcess/$targetPool.pass.fastq \
        $targetPool 3 1 200

### assign reads to each sample based on their barcodes
./bin/assignSample_multiple.py ./output_00_readsProcess/$targetPool.metadata.pass.tsv \
    ./sample.bc.tsv \
    ./output_00_readsProcess/$targetPool.metadata.pass.demultiplex.tsv \
    2 $targetPool

### generate reads for each sample in FASTQ format for mutation identification
./bin/splitSeq.py ./output_00_readsProcess/$targetPool.metadata.pass.demultiplex.tsv \
    ./output_00_readsProcess/$targetPool.pass.fastq \
    ./output_01_readsDemultiplex

### identify mutation in nanopore reads based on reference sequences provided
for eachS in `cat ./sample.list`
do
    echo $eachS

	### convert reads in FASTQ format to FASTA format
    ./bin/fq2fa.py ./output_01_readsDemultiplex/$eachS.fastq \
            ./output_01_readsDemultiplex/$eachS.corrected.fasta

	### identify mutation sites in each reads using muscle
    ./bin/statMutants_ind_multiple.py ./output_01_readsDemultiplex/$eachS.corrected.fasta \
        ./sample.bc.tsv \
        ./output_01_readsDemultiplex/$eachS.MSA.stat $eachS 

	### process output and merge indels if possible
    ./bin/processMutants.py ./output_01_readsDemultiplex/$eachS.MSA.stat \
        ./output_01_readsDemultiplex/$eachS.MSA.parsed.stat

	### generate statistics of all potential mutations
    ./bin/statPosition.py ./output_01_readsDemultiplex/$eachS.MSA.parsed.stat \
        ./output_01_readsDemultiplex/$eachS.mutation.stat

done
wait


