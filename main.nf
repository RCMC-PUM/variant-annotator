
include { validateParameters } from 'plugin/nf-schema'

params.sampleSheet = ''  // Path to the sample sheet
params.clinvarAnnotations = '' // Path to ClinVar VCF file
params.pathogenicityLevel = 'Pathogenic'  // User-defined pathogenicity filter
params.outputDir = '' // Path to output directory 
params.chatgptApiKey = ''// API key [temporary solution, move to secrets]


process ANNOTATE_VCF {
    input:
    tuple val(sample), file(vcf)
    path annots

    output:
    publishDir "$params.outputDir/samples/$sample/annotated_raw", mode: 'copy', overwrite: true, pattern: '*annot.vcf.gz'
    publishDir "$params.outputDir/samples/$sample/annotated_raw", mode: 'copy', overwrite: true, pattern: '*annot.vcf.gz.tbi'
    tuple val(sample), file("${vcf.name.replace('.vcf.gz', '.annot.vcf.gz')}")

    script:
    """
    bcftools index -t $vcf 
    bcftools index -t $annots

    bcftools annotate -a $annots -c CHROM,POS,REF,ALT,INFO -O z -o "${vcf.name.replace('.vcf.gz', '.annot.vcf.gz')}" $vcf
    """
}


process FILTER_VCF {

    input:
    tuple val(sample), file(annotated_vcf)

    output:
    publishDir "$params.outputDir/samples/$sample/annotated_filtered", mode: 'copy', overwrite: true, pattern: '*annot.filtered.vcf.gz'
    publishDir "$params.outputDir/samples/$sample/annotated_filtered", mode: 'copy', overwrite: true, pattern: '*annot.filtered.vcf.gz.tbi'
    tuple val(sample), file("${annotated_vcf.name.replace('.annot.vcf.gz', '.annot.filtered.vcf.gz')}")

    script:
    """
    bcftools view -i 'INFO/CLNSIG="${params.pathogenicityLevel}"' -O z -o "${annotated_vcf.name.replace('.annot.vcf.gz', '.annot.filtered.vcf.gz')}" $annotated_vcf
    """

}


process EXPORT_VARIANT_TABLE {
    input:
    tuple val(sample), file(filtered_vcf)

    output:
    publishDir "$params.outputDir/samples/$sample", mode: 'copy', overwrite: true, pattern: '*tsv'
    tuple val(sample), file("${filtered_vcf.name.replace('.annot.filtered.vcf.gz', '.annot.filtered.tsv')}")

    script:
    """

    if [[ "$filtered_vcf" == *sv* ]]; then
        f_query="%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n"

    elif [[ "$filtered_vcf" == *cnv* ]]; then
        f_query="%CHROM\t%POS\t%ID\t%INFO\t%FORMAT/CN\n"

    elif [[ "$filtered_vcf" == *repeats* ]]; then
        f_query="%CHROM\t%POS\t%INFO\n"

    elif [[ "$filtered_vcf" == *ploidy* ]]; then
        f_query="%CHROM\t%POS\t%INFO/END\t%FORMAT/DC\t%FORMAT/NDC\n"

    else
        f_query="%CHROM\t%POS\t%INFO/ALLELEID\t[%GT]\t%INFO\n"

    fi

    bcftools query -f "\$f_query" $filtered_vcf > "${filtered_vcf.name.replace('.annot.filtered.vcf.gz', '.annot.filtered.tsv')}"

    """
}


process INTERPRET_VARIANTS {
    input:
    tuple val(sample), file(filtered_tsv)
    val api_key

    output:
    publishDir "$params.outputDir/samples/$sample", mode: 'copy', overwrite: true, pattern: '*interpreted.txt'
    tuple val(sample), file("${filtered_tsv.baseName}.interpreted.txt"), optional: true

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    from openai import OpenAI
    import json

    client = OpenAI(api_key="${api_key}")

    # Load the TSV file
    try:
        df = pd.read_csv("${filtered_tsv}", sep="\t", header=None)

        # Prepare variant descriptions for analysis
        variants = df.to_dict(orient='records')
        
        prompt = f"Given the following variant data, prioritize the variants based on pathogenicity. If genotype information is available, include it in the analysis. Order the variants from most to least significant and suggest validation methods if necessary: {json.dumps(variants, indent=2)}"
        response = client.responses.create(
            model="gpt-4o-mini",
            instructions="You are a geneticist specializing in rare diseases.",
            input=prompt,
            )
        
        interpretation = response.output_text
        with open("${filtered_tsv.baseName}.interpreted.txt", "w") as handle:
            handle.write(interpretation)
        
    except pd.errors.EmptyDataError:
        print("No variants")
    
    """
}

workflow {

    validateParameters()

    sampleSheet = file(params.sampleSheet)
    annots = file(params.clinvarAnnotations)

    // Create the channel for samples from the CSV sample sheet file
    samples = Channel.fromPath(sampleSheet).splitCsv(header: true, sep: ",").map { row -> tuple(row.Sample_Name, file(row.Path) ) }
    
    // Run the annotation process and wait until all annotated VCF files are available
    annotated_vcf = ANNOTATE_VCF(samples, annots)

    // Filter annotated VCF files according to provided pathogenicity level
    filtered_vcf = FILTER_VCF(annotated_vcf)

    variant_tables = EXPORT_VARIANT_TABLE(filtered_vcf)
    INTERPRET_VARIANTS(variant_tables, params.chatgptApiKey)
}
