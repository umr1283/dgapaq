% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_genotype.R
\name{create_genotype_matrix}
\alias{create_genotype_matrix}
\title{create_genotype_matrix}
\usage{
create_genotype_matrix(
  sample_sheet,
  output_name = "raw_merged",
  output_directory = NULL,
  remove_homozygous_mutation = TRUE,
  split_multiallelic = FALSE,
  split_type = "both",
  ref_fasta = NULL,
  bin_path = list(vcftools = "/usr/bin/vcftools", bcftools = "/usr/bin/bcftools", tabix
    = "/usr/bin/tabix", bgzip = "/usr/bin/bgzip"),
  nb_cores = 1L
)
}
\arguments{
\item{sample_sheet}{A \code{data.frame}. A data frame with two mandatory columns \code{id} and \code{vcf_file}.}

\item{output_name}{A \code{character}. The prefix for output genotype data. Default is "raw_merged"}

\item{output_directory}{A \code{character}. The path to the output directory. Default is \code{NULL}.}

\item{remove_homozygous_mutation}{A \code{logical}. Does the homozygous mutation (i.e.: genotype composed
of different alternative alleles) should be removed. Default is \code{TRUE}.}

\item{split_multiallelic}{A \code{logical}. Does the multi-allelic variants should be splitted into
bi-allelic variants. Default is \code{FALSE}.}

\item{split_type}{A \code{character}. The type of variants should be splitted by \verb{bcftools norm}. Default is \code{both}.}

\item{ref_fasta}{A \code{character}. A path to the reference fasta file, it should be provided
if multiallelic INDELs need left-alignment after splitting. Default is \code{NULL} which lead to no left-alignment.}

\item{bin_path}{A \code{list(character)}. A list giving the binary path of \code{vcftools}, \code{bcftools}, \code{tabix} and \code{bgzip}.}

\item{nb_cores}{An \code{integer}. The number of CPUs to use. Default is \code{1}.}
}
\value{
A genotype matrix coded NA, 1, or 2 (additive model), with variants in rows and individuals in
columns. The first column \code{var_id} is composed of variant's chromosome, position, reference allele
and alternative allele, separated by "_".
}
\description{
Creates a raw genotype matrix with provided VCF files.
}