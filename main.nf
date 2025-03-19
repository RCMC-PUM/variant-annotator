include { validateParameters } from 'plugin/nf-schema'

params.sampleSheet = ''  // Path to the sample sheet
params.clinvarAnnotations = '' // Path to ClinVar VCF file
params.pathogenicityLevel = 'Pathogenic'  // User-defined pathogenicity filter
params.outputDir = '' // Path to output directory 
params.model_configuration = ''// Path to model configuration file


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
        bcftools filter -i 'NDC>1.1 || NDC<0.9' -O z -o "${vcf.getSimpleName()}.${caller}.filtered.pass.annot.vcf.gz" $vcf
    else
        bcftools filter -i 'INFO/CLNSIG="${params.pathogenicityLevel}"' -O z -o "${vcf.getSimpleName()}.${caller}.filtered.pass.annot.vcf.gz" $vcf
    fi
    """

}


process EXPORT_VARIANT_TABLE {
    input:
    tuple val(sample), file(vcf), val(caller), val(clinical_description)

    output:
    publishDir "$params.outputDir/samples/$sample/final_variant_tables", mode: 'copy', overwrite: true, pattern: '*tsv'
    tuple val(sample), file("${vcf.getSimpleName()}.${caller}.tsv"), val(caller), val(clinical_description)
    
    script:
    """
    if [[ "$caller" == "sv" ]]; then
        f_query="%CHROM\t%POS\t%REF\t%ALT\t[%SVTYPE]\t[%BND_DEPTH]\t[%GQ]\t[%GT]\t[%ALLELEID]\t[%CLNDISDB]\t[%CLNREVSTAT]\t[%CLNSIG]\t[%CLNVC]\t[%CLNDN]\n"
        f_cols="CHROM\tPOS\tREF\tALT\tSVTYPE\tBND_DEPTH\tGQ\tGT\tALLELEID\tCLNDISDB\tCLNREVSTAT\tCLNSIG\tCLNVC\tCLNDN"

    elif [[ "$caller" == "cnv" ]]; then
        f_query="%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%SM]\t[%CN]\t[%ALLELEID]\t[%CLNDISDB]\t[%CLNREVSTAT]\t[%CLNSIG]\t[%CLNVC]\t[%CLNDN]\n"
        f_cols="CHROM\tPOS\tREF\tALT\tGT\tSM\tCN\tALLELEID\tCLNDISDB\tCLNREVSTAT\tCLNSIG\tCLNVC\tCLNDN"

    elif [[ "$caller" == "repeats" ]]; then
        f_query="%CHROM\t%POS\t%REF\t%ALT\t[%VARID]\t[%GT]\t[%SO]\t[%REPCN]\t[%REPCI]\t[%ADSP]\t[%ADFL]\t[%ADIR]\t[%ALLELEID]\t[%CLNDISDB]\t[%CLNREVSTAT]\t[%CLNSIG]\t[%CLNVC]\t[%CLNDN]\n"
        f_cols="CHROM\tPOS\tREF\tALT\tVARID\tGT\tSO\tREPCN\tREPCI\tADSP\tADFL\tADIR\tALLELEID\tCLNDISDB\tCLNREVSTAT\tCLNSIG\tCLNVC\tCLNDN"

    elif [[ "$caller" == "ploidy" ]]; then
        f_query="%CHROM\t%QUAL\t%FILTER\t[%DC]\t[%NDC]\n"
        f_cols="CHROM\tQUAL\tFILTER\tDC\tNDC"

    else
        f_query="%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%AC]\t[%AF]\t[%AN]\t[%DP]\t[%FractionInformativeReads]\t[%GQ]\t[%GT]\t[%ALLELEID]\t[%CLNDISDB]\t[%CLNREVSTAT]\t[%CLNSIG]\t[%CLNVC]\t[%CLNDN]\n"
        f_cols="CHROM\tPOS\tREF\tALT\tQUAL\tAC\tAF\tAN\tDP\tFractionInformativeReads\tGQ\tGT\tALLELEID\tCLNDISDB\tCLNREVSTAT\tCLNSIG\tCLNVC\tCLNDN"

    fi

    echo -e "\$f_cols" > "${vcf.getSimpleName()}.${caller}.tsv"
    bcftools query -f "\$f_query" $vcf >> "${vcf.getSimpleName()}.${caller}.tsv"
    """
}


process INTERPRET_VARIANTS {
    input:
    tuple val(sample), file(tsv), val(caller), val(clinical_description)
    val model_config

    output:
    publishDir "$params.outputDir/samples/$sample/interpretation", mode: 'copy', overwrite: true, pattern: '*.json'
    tuple val(sample), file("*.json")

    script:
    """
    interpret.py $tsv $caller "${clinical_description ?: ''}" $model_config
    """
}


process CREATE_REPORTS {
    input:
    tuple val(sample), file(json_files)

    output:
    publishDir "$params.outputDir/samples/$sample", mode: 'copy', overwrite: true, pattern: '*.html'
    file "*.html"

    script:
    """
    reports.py $sample ${json_files.join(' ')}
    """
}

workflow {
    validateParameters()

    sampleSheet = file(params.sampleSheet)
    annots = file(params.clinvarAnnotations)

    // Create the channel for samples from the CSV sample sheet file
    samples = Channel.fromPath(sampleSheet).splitCsv(header: true, sep: ",").map { row -> tuple(row.Sample_Name, file(row.Path), row.Caller, row.Clinical_description ) }

    // Run the annotation process and wait until all annotated VCF files are available
    annots = INDEX_ANNOTS(annots)
    samples = FILTER_PASSING_VARIANTS(samples)
    annotated_vcf = ANNOTATE_VCF(samples, annots)

    // Filter annotated VCF files according to provided pathogenicity level
    filtered_vcf = FILTER_VCF(annotated_vcf)

    // Export selected variant tables
    variant_tables = EXPORT_VARIANT_TABLE(filtered_vcf)

    // LLM
    model_configuration = file(params.model_configuration)
    reports = INTERPRET_VARIANTS(variant_tables, model_configuration)

    // Export reports
    results = channel
    .fromPath("$params.outputDir/samples/*/interpretation/*.json", type: 'file')
    .map { file -> tuple(file.parent.parent.name, file) }
    .groupTuple()

    CREATE_REPORTS(results)
}
