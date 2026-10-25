// dart-nf: DART v2.0.1 end to end, one sample per samplesheet row.
//   profile -> predict (+decoys) -> MMseqs2 search (real + decoy) -> hits2emi (+null) -> damage-annotate -> report
// Every process publishes to ${params.outdir}/<sample>/.

process PROFILE {
    tag { sample }
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(profile_reads)

    output:
    tuple val(sample), path("${sample}.library.tdmg"), emit: tdmg
    path "${sample}.profile.json"
    path "${sample}.profile.log"

    script:
    """
    ${params.dart} profile -i ${profile_reads} -o ${sample}.profile.json \\
        --damage-summary ${sample}.library.tdmg --library-type ${params.library_type} \\
        -t ${task.cpus} > ${sample}.profile.log 2>&1
    """
}

process PREDICT {
    tag { sample }
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(reads), path(tdmg)

    output:
    tuple val(sample), path("${sample}.proteins.faa"), path("${sample}.proteins.decoy.faa"), emit: proteins
    tuple val(sample), path("${sample}.predictions.agd"), emit: agd
    tuple val(sample), path("${sample}.predictions.gff"), emit: gff
    path "${sample}.predict.log"

    script:
    """
    ${params.dart} predict -i ${reads} -o ${sample}.predictions.gff \\
        --fasta-aa ${sample}.proteins.faa --damage-index ${sample}.predictions.agd \\
        --library-profile ${tdmg} --decoy-every ${params.decoy_every} --domain ${params.domain} \\
        -t ${task.cpus} > ${sample}.predict.log 2>&1
    """
}

process MMSEQS_SEARCH {
    tag { sample }
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(proteins), path(decoys)
    path db_files
    val db_name
    path vtml20

    output:
    tuple val(sample), path("${sample}.hits.tsv"), path("${sample}.decoy_hits.tsv"), emit: hits
    path "${sample}.mmseqs.log"

    script:
    """
    mkdir -p tmp
    SEARCH="${params.search_params} --sub-mat ${vtml20} --seed-sub-mat ${vtml20} \\
        --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln \\
        --threads ${task.cpus}"
    ${params.mmseqs} easy-search ${proteins} ${db_name} ${sample}.hits.tsv tmp \$SEARCH > ${sample}.mmseqs.log 2>&1
    ${params.mmseqs} easy-search ${decoys} ${db_name} ${sample}.decoy_hits.tsv tmp \$SEARCH >> ${sample}.mmseqs.log 2>&1
    rm -rf tmp
    """
}

process HITS2EMI {
    tag { sample }
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(hits), path(decoy_hits), path(proteins), path(decoys), path(agd)
    val db_name

    output:
    tuple val(sample), path("${sample}.hits.emi"), path("${sample}.hits.emi.null"), emit: emi
    path "${sample}.hits2emi.log"

    script:
    """
    ${params.dart} hits2emi -i ${hits} -o ${sample}.hits.emi --damage-index ${agd} \\
        --decoy-hits ${decoy_hits} --decoy-fasta ${decoys} --db-label ${db_name} \\
        -t ${task.cpus} > ${sample}.hits2emi.log 2>&1
    """
}

process DAMAGE_ANNOTATE {
    tag { sample }
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(emi), path(null_model), path(agd)

    output:
    tuple val(sample), path("${sample}.annotated.tsv"), emit: annotated
    tuple val(sample), path("${sample}.gene_summary.tsv"), path("${sample}.protein_summary.tsv"), path("${sample}.damage_annotate.log"), emit: summaries

    script:
    """
    mkdir -p em_tape
    ${params.dart} damage-annotate --emi ${emi} --damage-index ${agd} \\
        -o ${sample}.annotated.tsv --gene-summary ${sample}.gene_summary.tsv \\
        --protein-summary ${sample}.protein_summary.tsv \\
        ${null_rule()} ${params.annotate_params} \\
        --em-tape-dir "\$PWD/em_tape" -t ${task.cpus} -v > ${sample}.damage_annotate.log 2>&1
    rm -rf em_tape
    """
}

process REPORT {
    tag { sample }
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(reads), path(gff), path(proteins), path(decoys), path(hits), path(decoy_hits), path(emi), path(null_model), path(annotated), path(gene_summary), path(protein_summary), path(annotate_log)

    output:
    path "${sample}.report.md"
    path "${sample}.versions.txt"

    script:
    """
    { ${params.dart} --version; ${params.mmseqs} version | sed 's/^/mmseqs /'; echo "nextflow ${workflow.nextflow.version}"; echo "dart-nf ${workflow.manifest.version} \$(git -C ${projectDir} rev-parse --short HEAD 2>/dev/null || echo no-git)"; } > ${sample}.versions.txt 2>&1
    ${params.dart} damage-info ${null_model} > null_info.txt 2>&1
    dart_nf_report.py --sample ${sample} --reads ${reads} --gff ${gff} --decoys ${decoys} \\
        --hits ${hits} --decoy-hits ${decoy_hits} --annotated ${annotated} \\
        --gene-summary ${gene_summary} --protein-summary ${protein_summary} \\
        --annotate-log ${annotate_log} --null-info null_info.txt --versions ${sample}.versions.txt \\
        -o ${sample}.report.md
    """
}

process COMPARE_EXPECTED {
    tag { sample }
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(annotated)
    path expected

    output:
    path "${sample}.compare.txt"

    script:
    """
    compare_annotated.py ${annotated} ${expected} | tee ${sample}.compare.txt
    test \${PIPESTATUS[0]} -eq 0
    """
}

// damage-annotate reject rule: the per-bin target-decoy FDR, or the fixed chance rule when --null_p_max is set.
def null_rule() {
    return params.null_p_max != null ? "--null-p-max ${params.null_p_max}" : "--null-fdr ${params.null_fdr}"
}

def resolve_path(base, p) {
    return p.startsWith('/') ? file(p) : base.resolve(p)
}

workflow {
    main:
    if (!params.samplesheet) error "set --samplesheet (CSV: sample,reads[,profile_reads])"
    if (!params.db) error "set --db (protein FASTA or MMseqs2 database prefix)"
    if ((params.decoy_every as int) < 1) error "--decoy_every must be >= 1: the null model needs decoys"

    def sheet = file(params.samplesheet)
    def base = sheet.parent
    def samples = channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row ->
            def reads = resolve_path(base, row.reads)
            def prof = row.profile_reads ? resolve_path(base, row.profile_reads) : reads
            tuple(row.sample, reads, prof)
        }
    def db_files = channel.fromPath("${params.db}*").collect()
    def db_name = file(params.db).name
    def vtml20 = file(params.vtml20)

    PROFILE(samples.map { it -> tuple(it[0], it[2]) })
    PREDICT(samples.map { it -> tuple(it[0], it[1]) }.join(PROFILE.out.tdmg))
    MMSEQS_SEARCH(PREDICT.out.proteins, db_files, db_name, vtml20)
    HITS2EMI(MMSEQS_SEARCH.out.hits.join(PREDICT.out.proteins).join(PREDICT.out.agd), db_name)
    DAMAGE_ANNOTATE(HITS2EMI.out.emi.join(PREDICT.out.agd))

    def report_in = samples.map { it -> tuple(it[0], it[1]) }
        .join(PREDICT.out.gff)
        .join(PREDICT.out.proteins)
        .join(MMSEQS_SEARCH.out.hits)
        .join(HITS2EMI.out.emi)
        .join(DAMAGE_ANNOTATE.out.annotated)
        .join(DAMAGE_ANNOTATE.out.summaries)
    REPORT(report_in)

    if (params.expected) {
        COMPARE_EXPECTED(DAMAGE_ANNOTATE.out.annotated, file(params.expected))
    }
}
