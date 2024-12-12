# **Nextflow Ribo-seq Pipeline**

This pipeline is designed for the analysis and interpretation of Ribosome profiling data. It is currently still being re-implemented into the Nextflow workflow language. 

## **Overview**

1. **Quality Control, Trimming and Filtering**: Trim reads with FastP and use a custom Bowtie2 index to remove unwanted RNA species, keeping mRNA ribosome protected fragments in the library.
2. **Alignment**: Use STAR to align the RPFs to the chosen reference genome and transcriptome.
3. **QC stats**: RiboseQC and multiQC for the output of relevant RPF statistics.
4. **ORF calling**: PRICE and ORFquant algorithms for calling new ORFs.
5. *ORF annotation*: Currently not implemented yet
6. *ORF expression*: Currently not implemented yet

## **Usage**

1. **Configure Parameters**:

    Configure the parameters located in `nextflow.config` to suit your analysis.

    - **Inputs**: Adjust settings to toggle pipeline components on and off.

    - **References and other files**: Define (file) paths to reference files, packages and other settings used in the pipeline.

    - **Containers**: Currently, we have loose containers for the entire pipeline that are described in `config/base.config`. We would like to provide a single for the pipeline in the future.


## Inputs

### Samplesheet

- The samplesheet is a critical input for the ribosome profiling pipeline.  It is a structured file,in CSV format, that lists all the samples to be processed. Each row represents a sample with detailed information required by the pipeline. Below is a breakdown of the expected columns and their constraints:

| **Column Name**  | **Description**                                                                     | **Required** | **Constraints**                                                                 |
|------------------|---------------------------------------------------------------------------------|--------------|---------------------------------------------------------------------------------|
| `subject_id`     | Unique identifier for the subject (e.g., patient).                              | Yes          | Must be a string without spaces.                                                |
| `sample_id`      | Unique identifier for the sample (e.g., sample barcode).                        | Yes          | Must be a string without spaces.                                                |
| `group_id`       | (Optional) Identifier for the sample group (e.g., cohort).                      | No           | Must be a string without spaces. If left empty, it must be omitted entirely.    |
| `sample_type`    | Type of sample: either `tumor` or `normal`.                                     | Yes          | Must be either `tumor` or `normal`.                                             |
| `sequence_type`  | Specifies the type of sequencing data: `rna` or `dna` or `ribo`.                | Yes          | Must be `ribo`, `rna` or `dna`.                                                  |
| `file_type`      | Format of the input files, e.g., `fastq`, `bam`, `cram`, `vcf`, `csv`, `tsv`.   | Yes          | Must be one of the supported formats: `fastq`, `bam`, `cram`, `vcf`, `csv`, etc.|
| `filename_1`     | Path to the first file (e.g., R1 FASTQ file for paired-end or single-end data).  | Yes          | Must be a valid file path with no spaces. File extension must match `file_type`. |
| `filename_2`     | (Optional) Path to the second file (e.g., R2 FASTQ file for paired-end data).    | No           | Must be a valid file path with no spaces. If not applicable, leave empty.       |

**Example samplesheet**:
```
subject_id,sample_id,group_id,sample_type,sequence_type,file_type,filename_1,filename_2
subject1,sampleA,cohort1,tumor,rna,fastq,/path/to/sampleA_R1.fastq.gz,/path/to/sampleA_R2.fastq.gz
subject2,sampleB,,normal,rna,bam,/path/to/sampleB.bam,
subject3,sampleC,cohort2,tumor,dna,vcf,/path/to/sampleC.vcf,
subject4,sampleA,cohort1,tumor,ribo,fastq,/path/to/sampleA.fastq.gz
```

Notes:

- All file paths (filename_1 and filename_2) must not contain spaces and should have extensions that match the declared file_type.

- For single-end data or non-FASTQ files (e.g., BAM, VCF), filename_2 can be omitted or left empty.

### Parameter specification


## Outputs

The pipeline generates output files including quality reports, trimmed reads, alignment results and ORF calls. The output directory is specified by the parameter `outdir` and has the following structure:

```
{outdir}
├── orfquant/
├── price/
├── fastp/
├── bowtie2/
├── fastqc/
└── star/
```

Additionally it produces Nextflow execution reports in `{project_folder}/log`

## Support and Contributions

- For issues or questions, [create an issue](https://github.com/VanHeeschTools/nf_ribo_pipeline/issues).
- Contributions to enhance this pipeline are welcomed through pull requests.
