conda activate qiime2

qiime tools import
–type ‘SampleData[PairedEndSequencesWithQuality]’
–input-path biofilm
–input-format CasavaOneEightSingleLanePerSampleDirFmt
–output-path demux-paired-end.qza

qiime cutadapt trim-paired
–i-demultiplexed-sequences demux-paired-end.qza
–p-cores 24
–p-front-f ACTCCTACGGGAGGCAGCAG
–p-front-r GGACTACHVGGGTWTCTAAT
–o-trimmed-sequences primer-trimmed-demux-paired-end.qza
–verbose

qiime demux summarize
–i-data demux-paired-end.qza
–o-visualization demux-paired-end.qzv

qiime dada2 denoise-paired
–i-demultiplexed-seqs demux-paired-end.qza
–p-trunc-len-f 285
–p-trim-left-f 26
–p-trunc-len-r 220
–p-trim-left-r 26
–p-n-threads 36
–o-representative-sequences rep-seqs.qza
–o-table table.qza
–o-denoising-stats stats.qza
–verbose

qiime feature-table summarize
–i-table table.qza
–o-visualization table.qzv
–m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs
–i-data rep-seqs.qza
–o-visualization rep-seqs.qzv


qiime feature-classifier extract-reads
–i-sequences silva-138-99-seqs.qza
–p-f-primer ACTCCTACGGGAGGCAGCAG
–p-r-primer GGACTACHVGGGTWTCTAAT
–p-n-jobs 24
–o-reads 338f806r-ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes
–i-reference-reads 338f806r-ref-seqs.qza
–i-reference-taxonomy silva-138-99-tax.qza
–o-classifier 338f806R-classifier.qza

qiime feature-classifier classify-sklearn
–i-classifier 338f806R-classifier.qza
–i-reads rep-seqs.qza
–o-classification taxonomy.qza

qiime metadata tabulate
–m-input-file taxonomy.qza
–o-visualization taxonomy.qzv

mkdir filtered qiime feature-table filter-samples
–i-table table.qza
–p-min-frequency 15010
–o-filtered-table filtered/filtered-table.qza

qiime feature-table filter-features
–i-table filtered/filtered-table.qza
–p-min-frequency 2
–o-filtered-table filtered/feature-frequency-filtered-table.qza

qiime taxa filter-table
–i-table filtered/feature-frequency-filtered-table.qza
–i-taxonomy taxonomy.qza
–p-include p__
–p-exclude mitochondria,chloroplast
–o-filtered-table filtered/table-with-phyla-no-mitochondria-chloroplast.qza


qiime taxa filter-table
–i-table filtered/table-with-phyla-no-mitochondria-chloroplast.qza
–i-taxonomy taxonomy.qza
–p-exclude "k__Archaea"
–o-filtered-table filtered/table-with-phyla-no-mitochondria-chloroplasts-archaea.qza

qiime taxa filter-table
–i-table filtered/table-with-phyla-no-mitochondria-chloroplasts-archaea.qza
–i-taxonomy taxonomy.qza
–p-exclude "k__Eukaryota"
–o-filtered-table filtered/table-with-phyla-no-mitochondria-chloroplasts-archaea-eukaryota.qza


qiime taxa filter-seqs
–i-sequences rep-seqs.qza
–i-taxonomy taxonomy.qza
–p-include p__
–p-exclude mitochondria,chloroplast
–o-filtered-sequences filtered/rep-seqs-with-phyla-no-mitochondria-chloroplast.qza


qiime taxa filter-seqs
–i-sequences filtered/rep-seqs-with-phyla-no-mitochondria-chloroplast.qza
–i-taxonomy taxonomy.qza
–p-exclude "k__Archaea"
–o-filtered-sequences filtered/rep-seqs-with-phyla-no-mitochondria-chloroplasts-archaea.qza


qiime taxa filter-seqs
–i-sequences filtered/rep-seqs-with-phyla-no-mitochondria-chloroplasts-archaea.qza
–i-taxonomy taxonomy.qza
–p-exclude "k__Eukaryota"
–o-filtered-sequences filtered/rep-seqs-with-phyla-no-mitochondria-chloroplasts-archaea-eukaryota.qza


mv filtered/table-with-phyla-no-mitochondria-chloroplasts-archaea-eukaryota.qza filtered/filtered-15010-table.qza

mv filtered/rep-seqs-with-phyla-no-mitochondria-chloroplasts-archaea-eukaryota.qza filtered/filtered-15010-rep-seqs.qza

qiime phylogeny align-to-tree-mafft-fasttree
–i-sequences filtered/filtered-15000-rep-seqs.qza
–o-alignment aligned-rep-seqs.qza
–o-masked-alignment masked-aligned-rep-seqs.qza
–o-tree unrooted-tree.qza
–o-rooted-tree rooted-tree.qza

mkdir tree

qiime tools export
–input-path unrooted-tree.qza
–output-path exported-tree

qiime tools export
–input-path filtered/filtered-15000-table.qza
–output-path filtered

biom convert -i filtered/feature-table.biom -o filtered/table.txt –to-tsv