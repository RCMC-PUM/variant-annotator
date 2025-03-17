# Variant Annotation and Interpretation Pipeline

## Overview
This repository contains a Nextflow pipeline for annotating, filtering, and interpreting genetic variants from VCF files. The pipeline integrates ClinVar annotations, filters variants based on user-defined pathogenicity levels, and utilizes OpenAI's API for variant interpretation.

## Features
- **VCF Annotation:** Uses `bcftools` to annotate VCF files with ClinVar data.
- **Filtering:** Filters variants based on pathogenicity.
- **Variant Table Export:** Extracts relevant variant information into tab-separated format.
- **AI-based Interpretation:** Leverages OpenAI API for prioritization and analysis of variants.

## Requirements
- [Nextflow](https://www.nextflow.io/)
- OpenAI API key

## Installation
Ensure Nextflow and docker are installed

## Usage
### Parameters
The pipeline requires the following parameters:
```yaml
params.sampleSheet = ''  # Path to the sample sheet (CSV format)
params.clinvarAnnotations = '' # Path to ClinVar VCF file
params.pathogenicityLevel = 'Pathogenic'  # Pathogenicity filter level
params.outputDir = '' # Output directory path
params.chatgptApiKey = '' # OpenAI API key (temporary solution, move to secrets)
```

### Running the Pipeline
To execute the pipeline, use:
```bash
nextflow run main.nf --sampleSheet path/to/sample_sheet.csv --clinvarAnnotations path/to/clinvar.vcf --pathogenicityLevel Pathogenic --outputDir path/to/output --chatgptApiKey YOUR_OPENAI_KEY
```

### Expected Output
The pipeline generates the following outputs for each sample:
- Annotated VCF files (`*.annot.vcf.gz`)
- Filtered VCF files (`*.annot.filtered.vcf.gz`)
- Variant tables (`*.annot.filtered.tsv`)
- AI-interpreted variant reports (`*.interpreted.txt`)

## Workflow Structure
1. **ANNOTATE_VCF**: Adds ClinVar annotations to VCF files.
2. **FILTER_VCF**: Filters variants based on the specified pathogenicity level.
3. **EXPORT_VARIANT_TABLE**: Converts filtered VCF files to tab-separated variant tables.
4. **INTERPRET_VARIANTS**: Uses OpenAI API to analyze and prioritize variants.

## Example Sample Sheet
The input sample sheet should be a CSV file structured as follows:
```csv
Sample_Name,Path
Sample1,sample1.vcf.gz
Sample2,sample2.vcf.gz
```

## Notes
- Ensure the OpenAI API key is handled securely.
- Consider moving API key management to environment variables or secrets.
- Future improvements include better handling of missing or malformed inputs.

## License
This project is licensed under the Apache2.0 License.

## Contact
For issues or contributions, please open an issue or submit a pull request.
