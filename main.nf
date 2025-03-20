import groovy.json.JsonOutput
include { validateParameters } from 'plugin/nf-schema'

params.sampleSheet = ''  // Path to the sample sheet
params.clinvarAnnotations = '' // Path to ClinVar VCF file
params.pathogenicityLevel = 'Pathogenic'  // User-defined pathogenicity filter
params.outputDir = '' // Path to output directory 

// NDC params
params.ndc_max = 1.1
params.ndc_min = 0.9

// HTML templates
params.template = 'templates/report.html'

process SAVE_PARAMS {
    output:
        publishDir "$params.outputDir", mode: 'copy', overwrite: true, pattern: 'params.json'
        path 'params.json'

    script:
      "echo '${JsonOutput.toJson(params)}' > params.json"
}


process INDEX_ANNOTS {
    input:
    path annots

    output:
    tuple file(annots), file("${annots.name}.tbi")

    script:
    """
    bcftools index -t $annots
    """
}


process SAVE_ANNOTS_PARAMS {
    input:
    tuple file(annots), file(annots_index)

    output:
    publishDir "$params.outputDir", mode: 'copy', overwrite: true, pattern: 'annots_params.json'
    path 'annots_params.json'

    script:
    """
    #!/usr/bin/python
    import vcfpy
    import json

    reader = vcfpy.Reader.from_path("${annots}")
    params = {line.key: line.value for line in reader.header.lines if line.key in ["fileformat", "fileDate", "source", "reference"]}

    with open("annots_params.json", "w") as handle:
        json.dump(params, handle, indent=2)
    """
}


process FILTER_PASSING_VARIANTS {
    input:
    tuple val(sample), file(vcf), val(caller), val(clinical_description)

    output:
    publishDir "$params.outputDir/samples/$sample/raw", mode: 'copy', overwrite: true, pattern: '*.pass.vcf.gz'
    tuple val(sample), file("${vcf.getSimpleName()}.${caller}.pass.vcf.gz"), val(caller), val(clinical_description)

    script:
    """
    if [[ "$caller" == "ploidy" ]]; then
        bcftools view -O z -o "${vcf.getSimpleName()}.${caller}.pass.vcf.gz" $vcf 
    else
        bcftools view -i 'FILTER=="PASS"' -O z -o "${vcf.getSimpleName()}.${caller}.pass.vcf.gz" $vcf 
    fi
    """


}


process ANNOTATE_VCF {
    input:
    tuple val(sample), file(vcf), val(caller), val(clinical_description)
    tuple file(annots), file(annots_index)

    output:
    publishDir "$params.outputDir/samples/$sample/annotated_raw", mode: 'copy', overwrite: true, pattern: '*.pass.annot.vcf.gz'
    tuple val(sample),file("${vcf.getSimpleName()}.${caller}.pass.annot.vcf.gz"), val(caller), val(clinical_description)

    script:
    """
    bcftools index -t $vcf
    bcftools annotate -a $annots -c CHROM,POS,REF,ALT,INFO -O z -o "${vcf.getSimpleName()}.${caller}.pass.annot.vcf.gz" $vcf
    """
}


process FILTER_VCF {
    input:
    tuple val(sample), file(vcf), val(caller), val(clinical_description)

    output:
    publishDir "$params.outputDir/samples/$sample/annotated_filtered", mode: 'copy', overwrite: true, pattern: '*.filtered.pass.annot.vcf.gz'
    tuple val(sample), file("${vcf.getSimpleName()}.${caller}.filtered.pass.annot.vcf.gz"), val(caller), val(clinical_description)

    script:
    """
    if [[ "$caller" == "ploidy" ]]; then
        bcftools filter -i 'NDC>=${params.ndc_max} || NDC<=${params.ndc_min}' -O z -o "${vcf.getSimpleName()}.${caller}.filtered.pass.annot.vcf.gz" $vcf
    else
        bcftools filter -i 'INFO/CLNSIG="${params.pathogenicityLevel}"' -O z -o "${vcf.getSimpleName()}.${caller}.filtered.pass.annot.vcf.gz" $vcf
    fi
    """

}


process CREATE_REPORTS {
    input:
    tuple val(sample), val(callers), file(vcf_files)
    file workflow_params
    file annots_params
    file template

    output:
    publishDir "$params.outputDir/samples/$sample", mode: 'copy', overwrite: true, pattern: '*.html'
    file "*.html"

    script:
    """
    reports.py $sample ${callers.join(',')} ${vcf_files.join(',')} $workflow_params $annots_params $template
    """
}

workflow {
    validateParameters()
    params_exported = SAVE_PARAMS()

    sampleSheet = file(params.sampleSheet)
    annots = file(params.clinvarAnnotations)

    // Create the channel for samples from the CSV sample sheet file
    samples = Channel.fromPath(sampleSheet).splitCsv(header: true, sep: ",").map { row -> tuple(row.Sample_Name, file(row.Path), row.Caller, row.Clinical_description ) }

    annots = INDEX_ANNOTS(annots)
    annots_params_exported = SAVE_ANNOTS_PARAMS(annots)

    // Run the annotation process
    samples = FILTER_PASSING_VARIANTS(samples)
    annotated_vcf = ANNOTATE_VCF(samples, annots)

    // Filter annotated VCF files according to provided pathogenicity level
    filtered_vcf = FILTER_VCF(annotated_vcf)

    // Export reports
    results = channel
    .fromPath("$params.outputDir/samples/*/annotated_filtered/*.vcf.gz", type: 'file')
    .map { file -> tuple(file.parent.parent.name, file.name.split("\\.")[1] ,file) }
    .groupTuple()

    template = file(params.template)
    CREATE_REPORTS(results, params_exported, annots_params_exported, template)
}
