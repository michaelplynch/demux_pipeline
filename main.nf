id_ch = Channel.from(1,2)
id_length_ch = id_ch.count()

seed_ch = Channel.from(2)

//id_ch = read_ch.counts().map { counts-> Channel.from(1..counts)}
id_ch.view()

params.dir="${projectDir}/data/output"

params.sim='./params_ccrcc.csv'

sims_ch = Channel.fromPath(params.sim)
    .splitCsv(header: true)
    .map { row ->
        tuple(row.key, row.pcdoub, tuple(*row.findAll { it.key != 'key' && it.key != 'pcdoub' }.collect { it.value }))
    }
    .view()


log.info """\
    S N P   S I M   P I P E L I N E
    ===================================

    sims: ${params.sim}
    """
    .stripIndent(false)


// Pull sample data from 10X
process PULLDATA {
    publishDir "${projectDir}/data/input"
    container "mplynch28/demux_sim"

    input:

    output:
    path("*/*.tsv"), emit: barcodes
    path("*/*.bam"), emit: bams
    
    shell:
    '''
    mkdir raji
    mkdir jurkat
    cd raji
    curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Raji/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Raji_count_sample_alignments.bam
    curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Raji/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Raji_count_sample_alignments.bam.bai
    curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Raji/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Raji_count_sample_barcodes.csv

    cd ../jurkat
    curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Jurkat/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Jurkat_count_sample_alignments.bam
    curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Jurkat/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Jurkat_count_sample_alignments.bam.bai
    curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Jurkat/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Jurkat_count_sample_barcodes.csv
    
    cd ..
    Rscript !{projectDir}/templates/checks.R '*/*.csv'
    '''
}
// Simulate doublets processes
process PARSE {
    publishDir "${params.dir}"
    container "mplynch28/demux_sim"

    input:
    path read_ch
    val id_ch

    output:
    path "parsed/${read_ch}_parsed.bam"
    
    shell:
    '''
    !{projectDir}/templates/parse_BAM_files.sh "!{read_ch}" "!{read_ch}_parsed.bam" K!{id_ch} parsed
    '''
}

process RENAME_BCS {
    publishDir "${params.dir}"
    container "mplynch28/demux_sim"

    input:
    path barcodes_ch
    val id_ch

    output:
    tuple val(id_ch), path("${barcodes_ch}_rn.tsv")

    shell:
    '''
    !{projectDir}/templates/parse_barcodes.sh "!{barcodes_ch}" !{id_ch}
    '''
}

process SUBSET_BARCODES {
    publishDir "${params.dir}"
    container "mplynch28/demux_sim"

    input:
    tuple val(key), val(doublets), val(ns)
    path barcodes
    val ids

    output:
    tuple val(key), path("*.tsv")


    shell:
    '''
    # Convert Nextflow variables into Bash arrays
    ids_array=(!{ids.join(' ')}) # Convert the `ids` tuple into a space-separated string and then into a Bash array
    ns_array=(!{ns.join(' ')})   # Convert the `ns` tuple into a space-separated string and then into a Bash array
    barcodes_array=(!{barcodes}) # Convert the `barcodes` path into a Bash array
    
    for i in $(seq 0 $((${#ids_array[@]} - 1))); do
        Rscript !{projectDir}/templates/subset_barcodes.R "${barcodes_array[i]}" "!{key}" "${ns_array[i]}" "${ids_array[i]}"
    done
    '''
}

process MERGE {
    cpus 4
    publishDir "${params.dir}"
    container "mplynch28/demux_sim"

    input:
    path parse_ch

    output:
    path "merged/merged.bam"

    shell:
    '''
    !{projectDir}/templates/merge_and_index_BAM.sh "!{parse_ch}" "merged.bam" "merged"
    '''
}

process MERGE_BCS {
    publishDir "${params.dir}"
    container "mplynch28/demux_sim"

    input:
    tuple val(key), path(barcodes)

    output:
    tuple val(key), path("barcodes_merged_${key}.tsv")

    shell:
    '''
    cat !{barcodes} > barcodes_merged_!{key}.tsv
    '''
}

process LOOKUP {
    publishDir "${params.dir}"
    container "mplynch28/demux_sim"

    input:
    tuple val(key), val(pcdoub), path(barcodes)

    output:
    tuple val(key), path("barcodes_merged_ccrcc_${key}_${pcdoub}pc.tsv"), emit: bcs_merged
    tuple val(key), val(pcdoub), path("lookup_table_doublets_ccrcc_${key}_${pcdoub}pc.tsv"), emit: lookup

    shell:
    '''
    Rscript !{projectDir}/templates/generate_awk_lookup_tables_doublets.R "!{projectDir}/data" "!{projectDir}/data" "!{projectDir}/data" "barcodes_merged_!{key}.tsv" !{pcdoub} !{key}
    '''
}

process SIMDOUB {
    publishDir "${params.dir}"
    container "mplynch28/demux_sim"

    input:
    val merge_ch
    tuple val(key), val(doublets), path(lookup_file)
    

    output:
    tuple val(key), val(doublets), path("doublets_ccrcc_${key}_${doublets}pc.bam"), path("doublets_ccrcc_${key}_${doublets}pc.bam.bai")
    
    shell:
    '''
    !{projectDir}/templates/parse_and_index_BAM_doublets.sh "!{projectDir}/data" "!{projectDir}/data" !{doublets} !{key} !{lookup_file} !{merge_ch}
    '''
}

// Run algorithms on sim data
process SOUP {
    cache true
    cpus 10
    publishDir "${params.dir}"
    container "shub://wheaton5/souporcell"
    
    input:
    tuple val(key), val(doublets), path(bam), path(bai), path(lookup)
    val ngroups

    output:
    tuple val(key), path("soup_${key}_${doublets}_1")

    shell:
    '''
    start=`date +%s`
     souporcell_pipeline.py \
        -i !{bam} \
        -b !{lookup} \
        -f !{params.ref} \
        -t 10 \
        -o "soup_!{key}_!{doublets}_1" \
        -k !{ngroups} \
        --common_variants !{params.common_variants} \
        --skip_remap SKIP_REMAP
    end=`date +%s`
    runtime=$((end-start))
    echo "$runtime" > soup_!{key}_!{doublets}_1/runtime.txt
    '''
}

process VIREO {
    cache true
    cpus 5
    publishDir "${params.dir}"
    container "mplynch28/demux_vireo"
    
    input:
    tuple val(key), val(doublets), path(bam), path(bai), path(lookup), val(seed)
    val ngroups

    output:
    path "vireo_${key}_${doublets}_${seed}"

    shell:
    '''
    start=`date +%s`
    cellsnp-lite -s !{bam} \
        -O cellsnp_out \
        -R "!{params.common_variants}" \
        -b !{lookup} \
        -p 5 \
        --gzip
    vireo -c cellsnp_out \
        -N !{ngroups} \
        -o "vireo_!{key}_!{doublets}_!{seed}" \
        -p 5 \
        --randSeed !{seed}
    end=`date +%s`
    runtime=$((end-start))
    echo "$runtime" > vireo_!{key}_!{doublets}_!{seed}/runtime.txt
    '''
}

workflow {
    // create data channels from download or local based on flag
    if (params.download_flag){
        PULLDATA()
        read_ch=PULLDATA.out.bams.flatten()
        barcodes_ch=PULLDATA.out.barcodes.flatten()
    }
    else {
        read_ch = Channel.fromPath(params.bam_path)
        barcodes_ch = Channel.fromPath(params.barcodes_path)
    }

    barcodes_ch.view()
    read_ch.view()

    // parse and merge bams
    parse_ch = PARSE(read_ch, id_ch)
    parse_ch.collect().view()
    merge_ch = MERGE(parse_ch.collect())

    // rename and subset barcodes & create lookup
    rename_bcs_ch = RENAME_BCS(barcodes_ch, id_ch)
    rename_bcs_ch
        .map { it[1] } // Extract only the second value of each tuple
        .collect()
        .set { barcodes_list_ch }

    //sims_ch.combine(id_tuple_ch).view()
    subset_bcs_ch=SUBSET_BARCODES(sims_ch,barcodes_list_ch,id_ch.collect())
    //subset_bcs_ch.view()
    merged_bcs_ch=MERGE_BCS(subset_bcs_ch)
    //merged_bcs_ch.view()
    sims_ch.view()
    sims_ch.map { tuple(it[0],it[1]) }
        .combine(merged_bcs_ch, by: 0)
        .set {doub_bcs_ch}
    //
    LOOKUP(doub_bcs_ch)

    // simulate doublets
    SIMDOUB(merge_ch, LOOKUP.out.lookup)
    test_ch = SIMDOUB.out
        .combine(LOOKUP.out.bcs_merged,by:0)
    runseeds_ch=test_ch.combine(seed_ch)
    //test_ch.view()
    //runseeds_ch.view()
    //id_length_ch.view()
    // run souporcell and Vireo
    SOUP(test_ch,id_length_ch)
    VIREO(runseeds_ch,id_length_ch)
}
