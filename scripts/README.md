# Tool wrapper scripts

Run these from the project root:

```bash
cd PhagePipeline
```

Use a valid FASTA file as input. The current `data/genomes/test/T4.fasta` is not valid FASTA; known valid T4 inputs include:

```bash
results/pipeline_test/bacphlip/T4.fasta
results/smoke_tests/input/T4.fasta
results/proteins/T4.faa
```

## Genome tools

```bash
scripts/checkv/run_checkv.sh results/pipeline_test/bacphlip/T4.fasta T4
scripts/phanotate/run_phanotate.sh results/pipeline_test/bacphlip/T4.fasta T4
scripts/pharokka/run_pharokka.sh results/pipeline_test/bacphlip/T4.fasta T4
scripts/bacphlip/run_bacphlip.sh results/pipeline_test/bacphlip/T4.fasta T4
scripts/rgi/run_rgi.sh results/pipeline_test/bacphlip/T4.fasta T4 contig
scripts/abricate/run_abricate.sh results/pipeline_test/bacphlip/T4.fasta T4 card
```

## Protein tools

```bash
scripts/acrdb/run_acrdb_blast.sh results/proteins/T4.faa T4
scripts/phagetailfinder/run_phagetailfinder.sh results/proteins/T4.faa T4
scripts/phagerbpdetect/run_phagerbpdetect_v4.sh results/proteins/T4.faa T4
```

## CDS-folder tool

PhageDPO expects a folder containing nucleotide CDS FASTA files:

```bash
scripts/phagedpo/run_phagedpo.sh results/smoke_tests/phagedpo/cds_input
```

## Current known blockers

- PhageTailFinder is missing `tools/PhageTailFinder/dbs/tail_pfam` and `tools/PhageTailFinder/dbs/nontail_pfam`.
- PhageRBPdetect v4 is missing `tools/PhageRBPdetection/data/RBPdetect_v4_ESMfine`.
- PhageDPO currently reaches the model but can fail with a sequence parsing/counting error on the available T4 files.
