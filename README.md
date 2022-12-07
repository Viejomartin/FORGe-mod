# FORGe (modified) tool for ranking variants and building an optimal graph genome

There is a Jupyter Notebook in the folder /src/ where you can interact with the scripts without the use of a console.

This modified script contains 2 scripts -- `rank_mod.py` and `build.py` -- as well as a helper script `vcf_to_1ksnp.py` for generating input files in the proper format.

## vcf_to_1ksnp.py ##

FORGe takes information about genetic variants in a `1ksnp` and an optional phasing file. The `1ksnp` format is less complex than VCF, but you can convert VCF to using the included `vcf_to_1ksnp.py`.

In a `1ksnp` file, each alternate allele is stored as a separate line with the following columns:

```
Chromosome	Position	Reference Allele	Alternate Allele	Population Frequency	??	# Alternates	Variant Name
```

The phasing file is a comma-separated table with a column for each individual and row for each variant, in the same order they appear in the 1ksnp file. An element contains the allele for the corresponding individual (column) and variant (row), with `0` indicating the reference allele and `k` indicating the kth alternate allele, in order, appearing in the 1ksnp file.

FORGe includes a helper script `vcf_to_1ksnp.py` to facilitate convertsion from a VCF file and set of ingroup individuals to FORGe variant and phasing files.  Ingroup individuals can be specified as a list of names to either include (`--ingroup`) or exclude (`--outgroup`).

Example usage:

```
./vcf_to_1ksnp.py --reference ref.fa --vcf variants.vcf --ingroup names.txt --out variants.1ksnp --individuals phasing.txt
```

## rank_mod.py ##

`rank_mod.py` takes a linear reference genome fasta, a `1ksnp` variant file and an optional file containing phasing information. The user also specifies a model using ` --method` (options are: `popcov` or `popcov-blowup`) and a window size using `--window-size`.

The user can indicate a specific chromosome to process using `--chrom`. This modified script does not include the hybrid ranking since the packages used in this algorithm are deprecated.

Example usage:

```
./rank_mod.py --method popcov --reference ref.fa --vars variants.1ksnp --window-size 100 --phasing phasing.txt --output ordered.txt
```

## build.py ##

`build.py` takes as input a set of ranked variants (output by `rank.py`) and a percentage of variants to include in the graph. It produces the necessary input files to build an index with HISAT2 or with Bowtie (ERG).

Example usage:

```
./build.py --reference ref.fa --vars variants.1ksnp --window-size 100 --hisat variants.snp --sorted ordered.txt --pct 10
```

## Full pipeline

From beginning to end, running the FORGe pipeline with HISAT2 might look like this:

```
./vcf_to_1ksnp.py --reference ref.fa --vcf variants.vcf --ingroup names.txt --out variants.1ksnp --individuals phasing.txt
./rank_mod.py --method popcov --reference ref.fa --vars variants.1ksnp --window-size 100 --phasing phasing.txt --output ordered.txt
./build.py --reference ref.fa --vars variants.1ksnp --window-size 100 --hisat variants.snp --sorted ordered.txt --pct 10
$HISAT_HOME/hisat2-build --snp variants.snp ref.fa index
```
