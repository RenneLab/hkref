#!/usr/bin/env nextflow
//Daniel Stribling  |  ORCID: 0000-0002-0649-9506
//Renne Lab, University of Florida
//
//Version: 
//  v01 - 2021-02-15

nextflow.enable.dsl=2

// --------------- Modules Setup ---------------
use_modules = ['helper.nf']
//'helper.nf', 'bowtie2.nf', 'seqkit.nf', 'samtools.nf', 'process_idxstats.py']

params.use_local_modules = true
params.update_modules = false

if( ( !params.use_local_modules || params.update_modules ) 
    && !params.containsKey('main_modules_dir') ) {
    log.error "No local modules defined, or update is requested."
    log.error "Please provide source NF modules directory: 'params.main_modules_dir'"
    exit 1
}
if( params.use_local_modules ) {
    params.modules_dir = "${projectDir}/bin_nf"
} else {
    params.modules_dir = params.main_modules_dir
}  

if( params.use_local_modules && params.update_modules && use_modules ) {   
    println "\n -- Copying Source Modules to: ${params.modules_dir} -- \n"
    exe_dir = file(params.modules_dir)
    exe_dir.mkdir()
    use_modules.each { mod_name -> 
        full_mod_name = "${mod_name}"
        println " -- -- ${mod_name}"
        source_file = file("${params.main_modules_dir}/${full_mod_name}", checkIfExists: true)
        source_file.copyTo(exe_dir)
    }
}
if( use_modules ) {
    println "\n -- Using Modules in: ${params.modules_dir} -- \n"
}


// --------------- Helper Functions Setup ---------------
include { 
    ensure_params_paths;     // Initialize, fill, and check default variables
    get_max_local_cpus;      // Set maximum cpus available per task for local execution
    workflow_details;        // Print initial workflow execution details
    execution_details;       // Print final workflow execution details
    test_params_key_allowed; // Ensure key in params with allowed options
    test_params_key;         // Ensure key in params
    test_params_keys;        // Ensure list of keys in params
    echo_process_details;    // Echo details of process to process stdout
    Standardize_Channel;     // Standardize channel to an N-tuple of items
    check_full_file_path;    // Check a file path for existence with verbose error
    return_as_list           // Ensure variable object is a list
    eexit;                   // Log error message and exit
    glob_to_files;           // Convert "glob" path to file path
} from params.modules_dir + "/helper"

// --------------- Set Default Process Paramter Values ---------------

params.main_prefix       = "hkref"
params.all_publish_dir   = "${projectDir}/${params.main_prefix}_working"
params.raw_publish       = true
params.raw_out_dir       = '01_raw_seqs'
params.clustered_publish = true
params.clustered_out_dir = '02_clustered_seqs'
params.mask_publish      = true
params.mask_out_dir      = '03_masked_seqs'
params.combined_publish  = true
params.combined_out_dir  = '04_combined_seqs'
params.mask_seed         = '11111'
params.final_publish_dir = "${projectDir}"
params.final_publish     = true
params.logs_publish      = true
params.logs_dir          = "${params.all_publish_dir}/logs"


// --------------- Define Processes ---------------


build_db_names = []
params.build_dbs.each { name, details ->
    build_db_names.add(name)
}

Channel.fromList(build_db_names)
      .set { db_names }

// -- Map dbs for Ensembl downloads
db_names.flatMap { db_name ->
    def out_lists = []
    db_details = params.build_dbs[db_name]
    if( db_details.containsKey('ensembl_sets') ) { 
        db_details['ensembl_sets'].each{ set_name ->
            def out_list = [db_name, set_name, db_details.db_settings, params.settings_output]
            out_lists.add(out_list)
        }
    } else {
        def out_list = [db_name, null, db_details.db_settings, params.settings_output]
        out_lists.add(out_list)
    }
    out_lists
    }
.set { ensembl_sets }

// -- Map dbs for Genbank downloads
if( !(params.containsKey('user_email') && params['user_email']) ) {
    eexit('Please Provide Parameter --user_email (params.user_email) for Entrez Downloads')
}
db_names.map { db_name ->
               db_details = params.build_dbs[db_name]
               db_genbank = params.build_dbs[db_name]['genbank_trnas']
               [db_name, params.user_email, db_details.db_settings, params.settings_output, db_genbank]
             }
       .filter { it -> it[4] }
       .set { genbank_sets }

// -- Map dbs for miRBase downloads
db_names.map { db_name ->
    db_details = params.build_dbs[db_name]
    db_file_lines = file(db_details.db_settings).readLines()
    db_file_lines.each{ line ->
        if( line.contains('mirbase_db_url') ) {
            mirbase_db_url = line.stripIndent() - ~/mirbase_db_url:/
            mirbase_db_url = mirbase_db_url.replaceAll('"', '').replaceAll("'", '')
            mirbase_db_url = mirbase_db_url.stripIndent()
        } else if( line.contains('mirbase_db_ver') ) {
            mirbase_db_ver = line.stripIndent() - ~/mirbase_db_ver:/
            mirbase_db_ver = mirbase_db_ver.replaceAll('"', '').replaceAll("'", '')
            mirbase_db_ver = mirbase_db_ver.stripIndent()
        }
    }
    [db_name,
     mirbase_db_ver,
     mirbase_db_url,
     db_details.db_settings, 
     params.settings_output
     ]
    }
.set { mirbase_sets }

// -- Map dbs for External DB downloads
db_names.flatMap { db_name ->
    db_details = params.build_dbs[db_name]
    def out_items = []
    if( db_details.containsKey('ext_dbs') && db_details.ext_dbs ) {
        db_details.ext_dbs.each { ext_db_name ->
            out_items.add([db_name, ext_db_name, db_details.db_settings, params.settings_output])
        }
    }
    out_items
}
.set { ext_db_sets }

// -- Create Settings Channels
output_settings = Channel.fromPath( params.settings_output )
                        .collect()

process DL_Ensembl {
    conda params.python_conda
    tag          { if( db_detail ) { "${db_name}-${db_detail}" } else { db_name } } 
    //label        'multi_cpu'
    label        'cluster'
    cache        'lenient'
    memory       2.GB
    cpus         1
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), val(db_detail), path(db_settings), path(out_settings)

    output:
    tuple val(db_name), path(all_outputs)
    path '.command.log'

    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${out_dir}/*",
               enabled: params.raw_publish
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log       = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    out_dir       = params.raw_out_dir
    all_outputs   = "${out_dir}/*"
    if( db_detail ) {
        use_query_flag = "--use_query ${db_detail}"
    } else {
        use_query_flag = ""
    }

    shell:
    '''
    set -v -H -o history
    echo "Database Name: !{db_name}"
    echo "Database Detail: !{db_detail}"
    echo "Database Settings File: !{db_settings}"
    echo "Output Settings File: !{out_settings}"

    mkdir !{out_dir}
    cd !{out_dir}

    python !{projectDir}/bin/dl_ensembl.py \\
                --db_settings_file ../!{db_settings} \\
                --out_settings_file ../!{out_settings} \\
                --no_cache \\
                !{use_query_flag}

    COMMAND="$( echo !! )"
    set +v +H +o history

    '''
}

process DL_GenBank {
    conda params.python_conda
    tag          { db_name }
    //label        'multi_cpu'
    //label        'cluster'
    cache        'lenient'
    memory       2.GB
    cpus         1
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), val(user_email), path(db_settings), path(out_settings), val(dl_genbank)

    output:
    tuple val(db_name), path(all_outputs)
    path '.command.log'

    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${out_dir}/*",
               enabled: params.raw_publish
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log       = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    out_dir       = params.raw_out_dir
    all_outputs   = "${out_dir}/*"

    shell:
    '''
    set -v -H -o history
    echo "Database Name: !{db_name}"
    echo "Database Settings File: !{db_settings}"
    echo "Output Settings File: !{out_settings}"

    mkdir !{out_dir}
    cd !{out_dir}

    python !{projectDir}/bin/dl_genbank.py \\
                --db_settings_file ../!{db_settings} \\
                --out_settings_file ../!{out_settings} \\
                --user_email !{user_email}

    COMMAND="$(echo !!)"
    set +v +H +o history

    '''
}

process DL_miRBase {
    conda        params.seqkit_conda
    tag          { db_name }
    //label        'multi_cpu'
    //label        'cluster'
    cache        'lenient'
    memory       2.GB
    cpus         1
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), val(mirbase_ver), val(mirbase_url), path(db_settings), path(out_settings)

    output:
    tuple val(db_name), path(all_outputs)
    path '.command.log'

    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${out_dir}/*",
               enabled: params.raw_publish
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log       = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    out_dir       = params.raw_out_dir
    out_file_name = "hkref_${db_name}_miRBase-${mirbase_ver}-full.fa"
    out_file      = "${out_dir}/${out_file_name}"
    all_outputs   = "${out_dir}/*"

    shell:
    '''
    echo "Database Name: !{db_name}"
    echo "Database Fasta: !{mirbase_url}"
    echo "Database Settings File: !{db_settings}"
    echo "Output Settings File: !{out_settings}"

    mkdir !{out_dir}
    ls

    set -v -H -o history
    wget -vc --no-check-certificate !{mirbase_url}
    gunzip -f mature.fa.gz
    mv mature.fa miRBase_mature.fa
    seqkit grep -r -p '^hsa' miRBase_mature.fa | \\
      seqkit seq --rna2dna | \\
      seqkit replace -p ' Homo sapiens ' -r '_miRBase_hsa-' | \\
      seqkit replace -p '.+MIMAT' -r 'MIMAT' | \\
      seqkit replace -p '$' -r '_microRNA' | \\
      seqkit sort -N -o !{out_file}

    SEQKIT_COMMAND="!!"

    echo """
      -- !{out_file_name} -- 
        File: '!{out_file_name}' was created on $(date).
        Mature miRNA sequences (version !{mirbase_ver}) were downloaded from miRBase:
        !{mirbase_url}
        were converted to DNA alphabet, 
        and identifiers were renamed to match Hyb convention.
        Command: ${SEQKIT_COMMAND}
    """ > !{out_file}.notes.txt
   
 
    set +v +H +o history
    '''
}

process DL_Ext_DB {
    conda        params.seqkit_conda
    tag          { "${db_name}-${ext_db_name}" }
    //label        'multi_cpu'
    //label        'cluster'
    cpus         1
    cache        'lenient'
    memory       2.GB
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), val(ext_db_name), path(db_settings), path(out_settings)

    output:
    tuple val(db_name), path(all_outputs)
    path '.command.log'

    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${out_dir}/*",
               enabled: params.raw_publish
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log         = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    out_dir         = params.raw_out_dir
    //out_file_name = "hkref_${db_name}_miRBase-${mirbase_ver}-full.fa"
    //out_file      = "${out_dir}/${out_file_name}"
    out_file_prefix = "${params.main_prefix}_${db_name}"
    all_outputs     = "${out_dir}/*"

    shell:
    '''
    echo "Database Name: !{db_name}"
    echo "External DB Name: !{ext_db_name}"

    echo "Get YAML Details:"
    DB_VER=$(python !{projectDir}/bin/detail_from_yaml.py --yaml_file !{db_settings} \\
        --detail ext_dbs !{ext_db_name} version) 
    echo -e "Ext DB Version: ${DB_VER}"
    DB_SOURCE=$(python !{projectDir}/bin/detail_from_yaml.py --yaml_file !{db_settings} \\
        --detail ext_dbs !{ext_db_name} source) 
    echo -e "Ext DB Source: ${DB_SOURCE}"
    DB_DETAIL=$(python !{projectDir}/bin/detail_from_yaml.py --yaml_file !{db_settings} \\
        --detail ext_dbs !{ext_db_name} detail) 
    echo -e "Ext DB Source: ${DB_DETAIL}"
    DB_URL=$(python !{projectDir}/bin/detail_from_yaml.py --yaml_file !{db_settings} \\
        --detail ext_dbs !{ext_db_name} url) 
    echo -e "Ext DB Details: ${DB_DETAIL}"
    DB_PREP=$(python !{projectDir}/bin/detail_from_yaml.py --yaml_file !{db_settings} \\
        --detail ext_dbs !{ext_db_name} prep_commands) 
    echo -e "Ext DB Prep - Preprocessing Commands: ${DB_PREP}"
    DB_PREPPED_NAME=$(python !{projectDir}/bin/detail_from_yaml.py --yaml_file !{db_settings} \\
        --detail ext_dbs !{ext_db_name} prepped_name) 
    echo -e "Ext DB Prep - Preprocessing Commands: ${DB_PREPPED_NAME}"
    DB_PROC=$(python !{projectDir}/bin/detail_from_yaml.py --yaml_file !{db_settings} \\
        --detail ext_dbs !{ext_db_name} proc_commands) 
    echo -e "Ext DB Prep - Processing Commands: ${DB_PROC}"
    DB_OUT_NAME=$(python !{projectDir}/bin/detail_from_yaml.py --yaml_file !{db_settings} \\
        --detail ext_dbs !{ext_db_name} out_name) 
    echo -e "Ext DB Prep - Output Name: ${DB_OUT_NAME}"

    mkdir !{out_dir}

    echo -e "\\nDownloading Resource\\n"
    set -v -H -o history
    wget -vc --no-check-certificate ${DB_URL}

    ${DB_PREP}

    echo -e "\\nProcessing Resource\\n"

    ls

    PROC_COMMAND="cat ${DB_PREPPED_NAME} | ${DB_PROC} > !{out_dir}/!{out_file_prefix}_${DB_OUT_NAME}"
    echo "${PROC_COMMAND}"
    eval ${PROC_COMMAND}

    echo """
      -- !{out_file_prefix}_${DB_OUT_NAME} -- 
        File: '!{out_file_prefix}_${DB_OUT_NAME}' was created on $(date).
        Sequences were acqured from:
          ${DB_URL}
        Details: ${DB_DETAIL}
        Version: ${DB_VER}
        Source: ${DB_SOURCE}

        Prepared by:
          ${DB_PREP}

          ${PROC_COMMAND}

    """ > !{out_dir}/!{out_file_prefix}_${DB_OUT_NAME}.notes.txt
    set +v +H +o history
    '''
}

process Proc_Transcripts {
    conda        params.seqkit_conda
    tag          { db_name }
    //label        'multi_cpu'
    //label        'cluster'
    cache        'lenient'
    memory       2.GB
    cpus         1
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), path(in_files)

    output:
    tuple val(db_name), path(all_outputs)
    path '.command.log'

    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${out_dir}/*",
               enabled: params.raw_publish
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log       = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    out_dir       = params.raw_out_dir
    all_outputs   = "${out_dir}/*"

    in_fastas      = in_files.findAll {fn -> ( "${fn}".endsWith(".fa") 
                                               || "${fn}".endsWith(".fasta") ) } 
    in_fastas_notes = []
    in_fastas.each { fn -> in_fastas_notes.add("${fn}.notes.txt") }
    in_details = in_files.findAll {fn -> ( "${fn}".endsWith(".csv") ) }
    in_details_notes = []
    in_details.each { fn -> in_details_notes.add("${fn}.notes.txt") }
    out_fasta_name   = "hkref_${db_name}_orig-transcripts.fa"
    out_fasta        = "${out_dir}/${out_fasta_name}"
    out_details_name = "${out_fasta_name}".replace('.fa', '.csv')
    out_details      = "${out_dir}/${out_details_name}"

    shell:
    '''
    set -v -H -o history
    echo "All input files: !{in_files}"
    echo "Input Fastas: !{in_fastas}"
    echo "Input Fastas Notes: !{in_fastas_notes}"
    echo "Input Details: !{in_details}"
    echo "Input Details Notes: !{in_details_notes}"

    mkdir !{out_dir}
    
    # -- Join and sort output sequences
    set -v -H -o history
    cat !{in_fastas.join(' ')} | \\
        seqkit sort --by-name --natural-order | \\
        seqkit sliding -g -s 655000 -W 655000 --suffix "" -o !{out_fasta}
    SEQKIT_COMMAND="!!"
    cat !{in_fastas_notes.join(' ')} > !{out_fasta}.notes.txt
    set +v +H +o history
    
    echo """
    
      -- !{out_fasta_name} --
        File: '!{out_fasta_name}' was created on $(date).
        by concatenating files: !{in_fastas}
        and sorting using seqkit 
        Command:
          ${SEQKIT_COMMAND}
    """ >> !{out_fasta}.notes.txt

    # -- Combine output details
    set -v -H -o history
    cat !{in_details.join(' ')} > !{out_details}
    DETAIL_COMMAND="!!"
    cat !{in_details_notes.join(' ')} > !{out_details}.notes.txt
    set +v +H +o history

    echo """
    
      -- !{out_details_name} --
        File: '!{out_details_name}' was created on $(date).
        by concatenating files: !{in_details}
        and sorting using seqkit 
        Command:
          ${DETAIL_COMMAND}
    """ >> !{out_details}.notes.txt
    echo -e "\\nDone\\n"
    '''
}

process Dedup_Transcripts {
    conda        params.cluster_conda
    tag          { db_name }
    label        'multi_cpu'
    label        'cluster'
    cache        'lenient'
    memory       7.GB
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), path(in_files)

    output:
    tuple val(db_name), path(all_outputs)
    path '.command.log'

    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${out_dir}/*",
               enabled: params.clustered_publish
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log       = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    out_dir       = params.clustered_out_dir
    out_file_name = "hkref_${db_name}_clust.fa"
    out_file      = "${out_dir}/${out_file_name}"
    all_outputs   = "${out_dir}/*"
    mem_flag      = ""

    transcript_fasta    = files_to_fasta(in_files)
    transcript_note     = transcript_fasta + '.notes.txt'
    transcript_details  = "${transcript_fasta}".replaceAll('.fa', '.csv')
    out_file_details    = "${out_file}".replaceAll('.fa', '.csv')

    if( task.memory ) {
        reduced_mem = ( task.memory * 0.9 )
        mem_flag = "-M ${reduced_mem.toUnit("MB")}"
    }

    shell:
    '''
    set -v -H -o history
    echo "All input files: !{in_files}"
    echo "Input Fasta: !{transcript_fasta}"
    echo "Input notes file: !{transcript_note}"

    mkdir !{out_dir}

    echo -e "\nRemoving duplicates and subsequences from file: !{transcript_fasta}"
    set -v -H -o history

    cd-hit-est -i !{transcript_fasta} \\
               -o !{out_file} \\
               -T !{task.cpus} \\
               !{mem_flag} \\
               -c 1.0 -t 1 -uS 0.0 \\
               -sc 1 \\
               -r 0 \\
               -g 0 \\
               -d 100 \\


    CLUSTER_COMMAND="!!"

      # -sc : sort .clustr file by cluster size
      # -c  : required identity proportion
      # -t  : Allowed redundance (default 2)
      # -uS : Maximum percent not aligned of smaller read
      # -r	: 1 or 0, default 1, by default do both +/+ & +/- alignments
      #       if set to 0, only +/+ strand alignment
      # -g  : 1 or 0, if 1 use slow clustering algorithm for best cluster result
      # -d  : length of description in clstr file

    #/apps/cdhit/4.6.8/psi-cd-hit/psi-cd-hit.pl \\
    #              -i !{transcript_fasta} \\
    #              -o !{out_file} \\
    #              -c 1.0 \\
    #              -aS 1.0 \\
    #              -g 1 \\
    #              -prog blastn \\

                  #-sc 1 \\
                  #-T !{task.cpus} \\
                  #!{mem_flag} \\


    # Ex: ./psi-cd-hit.pl -i db.fna -o db90.fna -c 0.9 -G 1 -g 1 -prog blastn -circle 1 -exec local -core 32
    # -aS	alignment coverage for the shorter sequence, default 0
 	# if set to 0.9, the alignment must covers 90% of the sequence
    #-g  (1/0), default 1
    #    by cd-hit's default algorithm, a sequence is clustered to the first 
    #    cluster that meet the threshold (fast cluster). If set to 1, the program
    #    will cluster it into the most similar cluster that meet the threshold
    #    (accurate but slow mode)
    #    but either 1 or 0 won't change the representatives of final clusters

    echo -e """

      -- !{out_file_name} --
        File: '!{out_file_name}' was created on $(date).
        by clustering file: !{transcript_fasta}
        using cd-hit command:
          ${CLUSTER_COMMAND}
    """ >> !{out_file}.notes.txt

    echo -e "Passing through transcript detail files:"
    cp -vPR !{transcript_details} !{out_file_details}
    cp -vPR !{transcript_details}.notes.txt !{out_file_details}.notes.txt

    echo -e "\\nDone\\n"
    '''
}

process Merge_Transcripts {
    conda        params.python_conda
    tag          { db_name }
    //label        'multi_cpu'
    //label        'cluster'
    cache        'lenient'
    memory       7.GB
    cpus         1
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), path(in_files)
    path(out_settings_file)

    output:
    tuple val(db_name), path(all_outputs)
    path '.command.log'

    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${out_dir}/*",
               enabled: params.clustered_publish
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log       = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    out_dir       = params.clustered_out_dir
    out_file_name = "hkref_${db_name}_clust_merged.fa"
    out_file      = "${out_dir}/${out_file_name}"
    all_outputs   = "${out_dir}/*"
    mem_flag      = ""

    transcript_fasta    = files_to_fasta(in_files)
    transcript_note     = transcript_fasta + '.notes.txt'
    transcript_det      = "${transcript_fasta}".replaceAll('.fa', '.csv')
    transcript_det_note = transcript_det + '.notes.txt'
    out_file_details    = "${out_file}".replaceAll('.fa', '.csv')

    shell:
    '''
    echo "All input files:    !{in_files}"
    echo "Input Fasta:        !{transcript_fasta}"
    echo "Input Fasta Notes:  !{transcript_note}"
    echo "Input Details:      !{transcript_det}"
    echo "Input Detail Notes: !{transcript_det_note}"
    echo "Output Settings:    !{out_settings_file}"

    mkdir !{out_dir}

    echo -e "\nMerging split sequences from file: !{transcript_fasta}"
    set -v -H -o history

    python !{projectDir}/bin/merge_split_seqs.py \\
           --in_file !{transcript_fasta} \\
           --in_notes !{transcript_note} \\
           --out_file !{out_file} \\
           --out_notes !{out_file}.notes.txt \\
           --out_settings_file !{out_settings_file} \\
                                        

    COMMAND="!!"

    set +v +H +o history

    echo -e "Passing through transcript detail files:"
    cp -vPR !{transcript_det} !{out_file_details}
    cp -vPR !{transcript_det_note} !{out_file_details}.notes.txt

    echo -e "\\nDone\\n"
    '''
}

process Cluster_miRNAs {
    conda        params.cluster_conda
    tag          { db_name }
    //label        'multi_cpu'
    //label        'cluster'
    cache        'lenient'
    memory       7.GB
    cpus         1
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), path(in_files)

    output:
    tuple val(db_name), path(all_outputs)
    path '.command.log'

    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${out_dir}/*",
               enabled: params.clustered_publish
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log       = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    out_dir       = params.clustered_out_dir
    all_outputs   = "${out_dir}/*"
    mem_flag      = ""

    in_fasta = files_to_fasta(in_files)
    in_note = files_to_note(in_files)
    out_file_name = "${file(in_fasta).getName()}".replaceAll('-full', '-clust')
    out_file      = "${out_dir}/${out_file_name}"

    if( task.memory ) {
        reduced_mem = ( task.memory * 0.9 )
        mem_flag = "-M ${reduced_mem.toUnit("MB")}"
    }

    shell:
    '''
    set -v -H -o history
    echo "All input files: !{in_files}"
    echo "Input Fastas: !{in_fasta}"
    echo "Input notes files: !{in_note}"

    mkdir !{out_dir}

    echo -e "\\nRemoving duplicates and subsequences from file: !{in_fasta}"
    set -v -H -o history
    cd-hit-est -i !{in_fasta} \\
               -o !{out_file} \\
               -T !{task.cpus} \\
               !{mem_flag} \\
               -c 1.0 -t 1 -uS 0.0 \\
               -sc 1 \\
               -r 0 \\
               -d 100 \\
               -p 1 \\
               -g 1 \\

    # -sc : sort .clustr file by cluster size
    # -c  : required identity proportion
    # -t  : Allowed redundance (default 2)
    # -uS : Maximum percent not aligned of smaller read
    # -r	: 1 or 0, default 1, by default do both +/+ & +/- alignments
    #       if set to 0, only +/+ strand alignment
    # -d  : .clustr file description length
    # -p  : add alignment overlap to .clustr file
    # -g  : enable slower algorithm

    #cd-hit -i !{in_fasta} \\
    #       -o !{out_file} \\
    #       !{mem_flag} \\
    #       -T !{task.cpus} \\
    #       -t 1 -d 100 -p 1 -sc 1 \\
    #       -g 1 -c 1 

    # -t  : Allowed redundance (default 2)
    # -d  : .clustr file description length
    # -p  : add alignment overlap to .clustr file
    # -sc : sort .clustr file by cluster size
    # -g  : enable slower algorithm
    # -c  : sequence identity proportion of shorter sequence

    CLUSTER_COMMAND="!!"

    cp -v !{in_note} !{out_file}.notes.txt
    echo """
      -- !{out_file_name} --
        File: '!{out_file_name}' was created on $(date).
        by clustering file: !{in_fasta}
        using cd-hit command:
          ${CLUSTER_COMMAND}
    """ >> !{out_file}.notes.txt


    REPORT_PROPS="0.95 0.90 0.85 0.80" # 0.75"
    echo -e "\\nGenerating clustering reports for proportions: ${REPORT_PROPS}"
    for prop in ${REPORT_PROPS}; do
      prop_name="${prop/\\./}"
      cd-hit-est -i !{in_fasta} \\
                 -o !{out_file_name}.cl${prop_name} \\
                 -T !{task.cpus} \\
                 !{mem_flag} \\
                 -c ${prop} -t 1 \\
                 -sc 1 \\
                 -r 0 \\
                 -d 100 \\
                 -p 1 \\
                 -g 1 \\

      #cd-hit -i !{in_fasta} \\
      #       -o !{out_file_name}.cl${prop_name} \\
      #       !{mem_flag} \\
      #       -T !{task.cpus} \\
      #       -t 1 -d 100 -p 1 -sc 1 \\
      #       -g 1 -c ${prop} 
    done

    cp -v *.clstr !{out_dir}
    
    echo -e "\\nDone\\n"
    '''
}

process Mask_miRNAs {
    conda        params.mask_conda
    tag          { db_name }
    label        'multi_cpu'
    label        'cluster'
    cache        'lenient'
    memory       7.GB
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), path(mirna_files), path(transcript_files)

    output:
    tuple val(db_name), path(mask_outputs)
    tuple val(db_name), path(combined_outputs)
    path '.command.log'

    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${mask_out_dir}/*",
               enabled: params.mask_publish
    publishDir "${params.all_publish_dir}", mode: 'copy',
               pattern: "${combined_out_dir}/*",
               enabled: params.combined_publish
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log          = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    mask_out_dir     = params.mask_out_dir
    combined_out_dir = params.combined_out_dir
    mask_outputs     = "${mask_out_dir}/*"
    combined_outputs = "${combined_out_dir}/*"
    seed             = params.mask_seed

    mirna_fasta         = files_to_fasta(mirna_files)
    mirna_note          = files_to_note(mirna_files)
    transcript_fasta    = files_to_fasta(transcript_files)
    transcript_note     = transcript_fasta + '.notes.txt'
    transcript_details  = "${transcript_fasta}".replaceAll('.fa', '.csv')
    mask_fasta          = "${mask_out_dir}/${file(transcript_fasta).getBaseName()}.fa"
    mask_fasta          = "${mask_fasta}".replaceAll('orig-transcripts', 'masked')
    mask_fasta_name     = file(mask_fasta).getName()
    mask_details        = "${mask_fasta}".replaceAll('.fa', '.csv')
    combined_fasta      = "${combined_out_dir}/${file(transcript_fasta).getBaseName()}.fa"
    combined_fasta      = "${combined_fasta}".replaceAll('orig-transcripts', 'masked-combined')
    combined_fasta_name = file(combined_fasta).getName()
    combined_details    = "${combined_fasta}".replaceAll('.fa', '.csv')

    bt2_db         = "${transcript_fasta}" - ~/.fa$/
    mirna_alns_sam = ("${mirna_fasta}" - ~/.fa$/) + ".sam"
    mirna_alns_bed = mirna_alns_sam.replace('.sam', '.bed')

    if( task.memory ) {
        reduced_mem = ( task.memory * 0.9 )
        split_mem = ( reduced_mem.div(task.cpus) )
        mem_flag = "-m ${split_mem.toUnit("MB")}"
    }

    shell:
    '''
    set -v -H -o history
    echo "All input mirna files: !{mirna_files}"
    echo "All input transcript files: !{transcript_files}"

    mkdir !{mask_out_dir}
    mkdir !{combined_out_dir}

    unset BOWTIE2_INDEXES 

    echo -e "Aligning miRNAs: !{mirna_fasta}"
    echo -e "  to Transcripts: !{transcript_fasta}"
    echo -e "  to find overlap."

    set -v -H -o history
    bowtie2-build --threads !{task.cpus} !{transcript_fasta} !{bt2_db}
    ALL_COMMANDS="!!"

    bowtie2 -x !{bt2_db} -f !{mirna_fasta} -S !{mirna_alns_sam} \\
             --threads !{task.cpus} \\
             --no-unal --xeq \\
             --norc -D 20 -R 3 -N 0 -L 16 -a --local -i S,1,0.50 \\
             --score-min L,17,0 --ma 1 --np 0 --mp 2,2 --rdg 5,1 --rfg 5,1

    COMMAND="!!"
    ALL_COMMANDS="${ALL_COMMANDS}\\n      ${COMMAND}"
    set +v +H +o history

    # ---- Create bed file from alignments ----
    echo -e "\\nProcessing alignments into BED file"

    set -v -H -o history
    sam2bed < !{mirna_alns_sam} | cut -f1-6 > !{mirna_alns_bed}
    COMMAND="!!"
    ALL_COMMANDS="${ALL_COMMANDS}\\n    ${COMMAND}"
    set +v +H +o history

    #echo "reducing size of miRNA-regions by 6 nt"
    #set -v -H -o history
    #bedtools slop -i ${ALN_BED} -g ${SEQ_GENOME} -b "-6" > ${ALN_BED_TRIM} #/.bed/_unsort.bed}
    #COMMAND="!!"
    #ALL_COMMANDS="${ALL_COMMANDS}\\n${COMMAND}"
    #set +v +H +o history

    # ---- Create masked fasta ----
    set -v -H -o history
    bedtools maskfasta -fi !{transcript_fasta} -bed !{mirna_alns_bed} -fo !{mask_fasta}
    COMMAND="!!"
    ALL_COMMANDS="${ALL_COMMANDS}\\    ${COMMAND}"
    set +v +H +o history

    # ---- Create notes file for !{mask_fasta} ----
    cp -v !{transcript_fasta}.notes.txt !{mask_fasta}.notes.txt 
    
    echo -e "Passing through transcript detail files:"
    cp -vPR !{transcript_details} !{mask_details}
    cp -vPR !{transcript_details}.notes.txt !{mask_details}.notes.txt

    echo -e """
    
      -- !{mask_fasta_name} --
        File: '!{mask_fasta_name}' was created on $(date).
        by commands: 
          ${ALL_COMMANDS} 
    """ >> !{mask_fasta}.notes.txt

    # ---- Create combined fasta file ----
    cat !{mask_fasta} !{mirna_fasta} > !{combined_fasta}
    
    # ---- Create notes file for combined fasta file ----
    cat !{mask_fasta}.notes.txt !{mirna_note} > !{combined_fasta}.notes.txt
    echo -e """
    
      -- !{combined_fasta_name} --
        File: '!{combined_fasta_name}' was created on $(date).
        by concatenating files: !{mask_fasta} !{mirna_fasta} 
    """ >> !{combined_fasta}.notes.txt

    echo -e "Passing through transcript detail files:"
    cp -vPR !{transcript_details} !{combined_details}
    cp -vPR !{transcript_details}.notes.txt !{combined_details}.notes.txt

    #combined_nocoord_fasta_name

    echo -e "\\nDone\\n"
    '''
}

process Zip_Outputs {
    //conda        params.cluster_conda
    tag          { db_name }
    //label        'multi_cpu'
    //label        'cluster'
    cache        'lenient'
    memory       7.GB
    //cpus         { task.label.contains('multi_cpu') ? params.max_local_cpus : 1 }
    beforeScript { echo_process_details(task) }
    //echo         true
    //stageInMode  'copy'
    
    input:
    tuple val(db_name), path(in_files), path(db_settings), path(out_settings)

    output:
    tuple val(db_name), path(all_outputs)
    path '.command.log'

    publishDir "${params.final_publish_dir}", mode: 'copy',
               pattern: "${all_outputs}",
               enabled: params.final_publish
               
    publishDir "${params.logs_dir}", mode: 'copy',
               pattern: ".command.log", saveAs: { out_log },
               enabled: params.logs_publish

    script:
    out_log       = "${task.tag}.${task.process}.nflog.txt".replaceAll(":","_")
    out_dir       = params.clustered_out_dir
    all_outputs   = "*.{zip,gz,notes.txt}"

    in_fasta       = files_to_fasta(in_files)
    in_fasta_note  = "${in_fasta}.notes.txt"
    in_detail      = "${in_fasta}".replaceAll('.fa', '.csv')
    in_detail_note = "${in_detail}.notes.txt"
    final_fasta       = "${params.main_prefix}_${db_name}.fa"
    final_fasta_note  = "${params.main_prefix}_${db_name}.fa.notes.txt"
    final_detail      = "${params.main_prefix}_${db_name}.csv"
    final_detail_note = "${params.main_prefix}_${db_name}.csv.notes.txt"

    shell:
    '''
    set -v -H -o history
    echo "All input files: !{in_fasta}"

    set -v -H -o history
    gzip --force --best --to-stdout !{in_fasta} > !{final_fasta}.gz 
    gzip --force --best --to-stdout !{in_detail} > !{final_detail}.gz 

    cp -v !{in_fasta_note} !{final_fasta_note}
    cp -v !{in_detail_note} !{final_detail_note}

    ls
    
    echo -e "\\nDone\\n"
    '''
}

workflow {
    // -- Transcript Processing
    DL_Ensembl(ensembl_sets)
    //DL_Ensembl(Channel.empty())
    DL_GenBank(genbank_sets)
    DL_Ext_DB(ext_db_sets)
    DL_Ensembl.out[0]
        .mix(DL_GenBank.out[0], DL_Ext_DB.out[0])
        .groupTuple()
        .map { it ->
            def all_files = []
            it[1].each { files -> 
                all_files.addAll(files) 
                }
            [it[0], all_files]
        }
        .set { joined_transcripts }
    Proc_Transcripts(joined_transcripts)
    Dedup_Transcripts(Proc_Transcripts.out[0])
    Merge_Transcripts(Dedup_Transcripts.out[0], output_settings)

    DL_miRBase(mirbase_sets)
    Cluster_miRNAs(DL_miRBase.out[0])
                .set { clustered_mirnas } 

    DL_miRBase.out[0].mix ( clustered_mirnas )
             .cross ( Merge_Transcripts.out[0] )
             .map { items ->
                 db_name = items[0][0]
                 mirna_files = items[0][1]
                 transcript_files = items[1][1]
                 [db_name, mirna_files, transcript_files] 
             } \
    | Mask_miRNAs
    Mask_miRNAs.out[1].map { db_name, out_files ->
        db_details = params.build_dbs[db_name]
        [db_name, out_files, db_details.db_settings, params.settings_output]
    } \
      | Zip_Outputs
        

}


def files_to_file( in_files, suffixes, descriptor ) {
    def specific_files = []
    suffixes.each { suffix ->
        def find_files = in_files.findAll {fn -> ( "${fn}".endsWith(suffix) ) }
        specific_files.addAll(find_files)
    }
    if( specific_files.size() < 1 ) { 
        eexit("<1 ${descriptor}: ${in_files}")
    } else if( specific_files.size() > 1 ) { 
        eexit(">1 ${descriptor}: ${in_files}")
    }
    specific_files[0]
}

def files_to_fasta( in_fastas ) {
    files_to_file( in_fastas, ['.fasta', '.fa'], 'Fasta Files' )
}

def files_to_note( in_notes ) {
    files_to_file( in_notes, ['.notes.txt', '.notes'], 'Notes Files' )
}




