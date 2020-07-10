#' Compute quality-control of genotyping array (PLINK format) using a rmarkdown template.
#'
#' @param input_files A `character`. The path to the plink files.
#'     The path should contains the file name without the extension, i.e., without `*.bed`, `*.bim` or `*.fam`.
#' @param output_directory A `character`. The path to the output directory.
#' @param cohort_name A `character`. The name of the studied cohort / population.
#' @param output_file A `character`. The name of the html file produced.
#' @param array A `character`. The array name, e.g., "Illumina Omni2.5".
#' @param callrate_samples A `numeric`. The call rate threshold for samples, under which samples are excluded.
#'     Default is `0.95`.
#' @param callrate_snps A `numeric`. The call rate threshold for probes, under which probes are excluded.
#'     Default is `0.95`.
#' @param sex_threshold  A `vector(numeric)`. The sex threshold to use to define the genotypic sex in PLINK.
#' @param heterozygosity_threshold A `numeric`. The heterozygosity threshold for samples
#'    (number of standard deviation from the mean), under/above which samples are excluded.
#'     Default is `4`.
#' @param maf_threshold A `numeric`. The minor allele frequency under which variants are considered "rare".
#'     Default is `0.01`.
#' @param hwe_pvalue A `numeric`. The p-value threshold for Hardy-Weinberg equilibrium test.
#'     Default is `0.0001`.
#' @param includes_relatives A `logical`. Does the data contain related samples?
#'     Default is `FALSE`.
#' @param mendelian_samples A `numeric`. The Mendel error rate threshold above which samples are excluded.
#'     Default is `0.05`.
#' @param mendelian_snps A `numeric`. The Mendel error rate threshold above which variants are excluded.
#'     Default is `0.1`.
#' @param ibd_threshold A `numeric`. The threshold for IBD (identical by descent) above which
#'     samples are characterised as relatives.
#'     Default is `0.2`.
#' @param population A `character`. The ethnicity of the studied population if known, e.g., `"EUR"`.
#'     Default is `NULL`.
#' @param pca_components A `numeric`. The number of principal components to be computed.
#'     Default is `10`.
#' @param pca_threshold A `numeric`. The threshold to define outliers on the principal component analysis,
#'     the as number of standard deviation from the cohort centroid.
#'     Default is `3`.
#' @param max_labels A `numeric`. The maximum number of labels to show on plots. Default is `15`.
#' @param check_bim_script A `character`. The PERL script to use to check PLINK files to allow later imputation.
#'     Default is `system.file("perl", "HRC-1000G-check-bim.pl", package = "dgapaq")`.
#'     Script from https://www.well.ox.ac.uk/~wrayner/tools/.
#' @param imputation_ref A `character` giving the imputation panel to check against, *e.g.*, `1KG` (default) or `HRC`.
#' @param imputation_panel A `character`. The `*.legend` file from 1,000 Genome project or `*.tab` file from HRC.
#'     Default is `"1000GP_Phase3_combined.legend"`.
#'     `*.legend` file can be downloaded from https://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz.
#' @param ref1kg_panel A `character`. The `*.panel` file from 1,000 Genome project.
#'     Default is `"integrated_call_samples_v3.20130502.ALL.panel"`.
#'     `integrated_call_samples_v3.20130502.ALL.panel` can be downloaded from
#'     ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/.
#' @param ref1kg_fasta A `character`. The `*.fasta` file from 1,000 Genome project.
#'     Default is `NULL`.
#' @param ref1kg_population A `character`. The `*.tsv` file from 1,000 Genome project describing samples and ethnicity.
#'     Default is `"20131219.populations.tsv"`.
#'     The `*.tsv` file can be downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv.
#' @param ref1kg_genotypes A `character`. The PLINK files from 1,000 Genome project.
#'     Default is `NULL`.
#'     `*vcf` files can be downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/.
#' @param bin_path A `list(character)`. A list giving the binary path of `bcftools`, `bgzip`, `gcta` and `plink1.9`.
#' @param title A `character`. The report's title. Default is `paste(array, "Array Quality-Control")`.
#' @param author_name A `character`. The author's name to be printed in the report.
#'     Default is `Unknown`.
#' @param author_affiliation A `character`. The affiliation to be printed in the report.
#'     Default is `NULL`.
#' @param author_email A `character`. The email to be printed in the report.
#'     Default is `NULL`.
#' @param encoding A `character`. The encoding to be used for the html report.
#'     Default is `"UTF-8"`.
#' @param ... Parameters to pass to `rmarkdown::render()`.
#'
#' @return NULL
#'
#' @export
qc_plink <- function(
  input_files = NULL,
  output_directory = NULL,
  cohort_name = "COHORT",
  output_file = paste0(cohort_name, "_QC.html"),
  array = NULL,
  callrate_samples = 0.95,
  callrate_snps = 0.95,
  sex_threshold = c(0.2, 0.8),
  heterozygosity_threshold = 4,
  maf_threshold = 0.01,
  hwe_pvalue = 0.0001,
  includes_relatives = FALSE,
  mendelian_samples = 0.05,
  mendelian_snps = 0.1,
  ibd_threshold = 0.2,
  population = NULL,
  pca_components = 10,
  pca_threshold = 3,
  max_labels = 15,
  check_bim_script = system.file("perl", "HRC-1000G-check-bim.pl", package = "dgapaq"),
  imputation_ref = "1KG",
  imputation_panel = "1000GP_Phase3_combined.legend",
  ref1kg_panel = "integrated_call_samples_v3.20130502.ALL.panel",
  ref1kg_fasta = NULL,
  ref1kg_population = "20131219.populations.tsv",
  ref1kg_genotypes = NULL,
  bin_path = list(bcftools = "/usr/bin/bcftools", bgzip = "/usr/bin/bgzip", plink = "/usr/bin/plink1.9", gcta = "/usr/bin/gcta64"),
  title = paste(array, "Array Quality-Control"),
  author_name = "Unknown",
  author_affiliation = NULL,
  author_email = NULL,
  encoding = "UTF-8",
  ...
) {
  message_prefix <- "[dgapaq] "

  message(message_prefix, "Quality-Control started ...")

  invisible(file.copy(
    from = system.file("rmarkdown", "templates", "qc_plink", "skeleton.Rmd", package = "dgapaq"),
    to = file.path(tempdir(), "qc_plink.Rmd"),
    overwrite = TRUE
  ))

  if (!all(sapply(
    X = strsplit(gsub("\\([^()]+\\)|\n| ", "", utils::packageDescription("dgapaq", fields = "Suggests")), ",")[[1]],
    FUN = function(x) nchar(system.file(package = x)) > 0
  ))) {
    stop(message_prefix, "Packages listed in 'Suggests' field are needed to run this function.")
  }

  rmarkdown::render(
    input = file.path(tempdir(), "qc_plink.Rmd"),
    output_file = output_file,
    output_dir = output_directory,
    encoding = encoding,
    params = list(
      title = title,
      author_name = author_name,
      author_affiliation = author_affiliation,
      author_email = author_email,
      input_files = input_files,
      output_directory = output_directory,
      cohort_name = cohort_name,
      array = array,
      callrate_samples = callrate_samples,
      callrate_snps = callrate_snps,
      sex_threshold = sex_threshold,
      heterozygosity_threshold = heterozygosity_threshold,
      maf_threshold = maf_threshold,
      hwe_pvalue = hwe_pvalue,
      includes_relatives = includes_relatives,
      mendelian_samples = mendelian_samples,
      mendelian_snps = mendelian_snps,
      ibd_threshold = ibd_threshold,
      population = population,
      pca_components = pca_components,
      pca_threshold = pca_threshold,
      max_labels = max_labels,
      check_bim_script = check_bim_script,
      imputation_ref = imputation_ref,
      imputation_panel = imputation_panel,
      ref1kg_panel = ref1kg_panel,
      ref1kg_fasta = ref1kg_fasta,
      ref1kg_population = ref1kg_population,
      ref1kg_genotypes = ref1kg_genotypes,
      bin_path = bin_path
    ),
    ...
  )

  message(message_prefix, "Quality-Control ended!")

  invisible()
}
