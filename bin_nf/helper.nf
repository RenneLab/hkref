#!/usr/bin/env nextflow
//Daniel Stribling  |  ORCID: 0000-0002-0649-9506
//Renne Lab, University of Florida
//
//Version: 
//  v20210201

nextflow.enable.dsl=2

// --------------- Setup Default Variables: ---------------
params.main_dir = "/blue/renne/${workflow.userName}"
params.proj_dir = "${params.main_dir}"
params.task_dir = "${workflow.launchDir}"
params.verbose  = false
params.debug    = false
params.err_on_dup = true
params.ch_card    = 5
params.front_pad  = 4
params.desc_pad   = 17
params.mem_reduce = 0.9
params.seed_chars = 8
params.local_module_dir = "bin_nf"

// --------------- Nextflow Variable Manipulation Functions: ---------------

// Return object as a list, wrapping in list as needed
def return_as_list(item) {
    if( item instanceof List ) {
        return item
    } else {
        return [item]
    }
}

// Reduce Memory Value
def reduce_memory(mem_str) {
    def split_string = mem_str.split(' ')
    def orig_val = split_string[0] as int
    def final_val = (orig_val * params.mem_reduce) as int
    return "${final_val} ${split_string[1]}"
}

// Divide Memory Value by Given Number
def divide_memory(mem_str, divisor_str) {
    def divisor = divisor_str as int
    def split_string = mem_str.split(' ')
    def orig_val = split_string[0] as int
    //def final_val = (orig_val  params.mem_reduce) as int
    def final_val = orig_val.intdiv(divisor) as int
    return "${final_val} ${split_string[1]}"
}

// Report error message and exit
def eexit(message) {
    println ""
    log.error message
    println ""
    sleep(1000)
    exit 1
}

// Convert "glob" path to file path
def glob_to_files(file_glob) {
    def files_list = file(file_glob)
    if( files_list.isEmpty() ) {
        log.error "Files not found for glob:"
        log.error "-   ${file_glob}"
        println ""
        exit 1
    }
    return files_list.sort()
}

// Check for the existence of a file by path, and verbosely report if absent
def check_full_file_path(file_name) {
    if( !file(file_name).exists()) {
        log.error ""
        log.error "File Not Found:  ${file_name}"
        log.error ""
        build_subpath = ""
        // Regex to remove starting and trailing slashes of path
        use_string = file_name - ~/^\// - ~/\/$/
        use_string.split('/').each { seg ->
            build_subpath += "/" + seg
            if( file(build_subpath, checkIfExists: false).exists() ) {
                log.error "    Exists:         " + build_subpath //+ "\n"
            } else {
                log.error "    Does not Exist: " + build_subpath //+ "\n"
            }
        }
        println ""
        exit(1)
    }
}

// Print a padded item linewise with description
def pad_print_item_description(data, description) {
    // params.front_pad is initial padding
    // params.desc_pad is description padding
    if( data instanceof List ) {
        def message = ""
        message += "-".multiply(params.front_pad)
        message += "${description}:".padRight(params.desc_pad)
        log.info message + data[0]
        if( data.size() > 1 ) {
            data.subList(1, data.size()).each { item ->
                log.info "-".multiply(params.front_pad) + " ".multiply(params.desc_pad) + item
            }
        }
    } else {
        log.info "${description}:".padRight(params.desc_pad) + data
    }
}

// Convert input to a string, hash string, and return the first N bytes as a numeric seed
def generate_seed (in_var) {
    def ret_str = ''
    "${in_var}".digest("SHA-1").decodeHex().each { it
        ret_str += "${it}" - ~/-/
    }
    Long.valueOf(ret_str[0..(params.seed_chars - 1)])
}

// --------------- Nextflow Script Startup Functions: ---------------

// Create mapping with expected parameters with form: tuple('key', def_value, 'def_value_description')
def return_default_params() {
    def default_params = [
        tuple('main_dir', workflow.homeDir, 'environment default'),
        tuple('proj_dir', workflow.launchDir, 'script launch directory'),
        tuple('task_dir', workflow.launchDir, 'script launch directory'),
        tuple('in_dir', params.task_dir, 'match task directory'),
        //tuple('pub_mode_mid', 'symlink', 'Nextflow default'),
        //tuple('pub_mode_final', 'symlink', 'Nextflow default'),
    ]
    return default_params
}

// Fill missing default paramaters in params.
def fill_default_params(params, default_params) {
    default_params.each { key, def_value, description ->
        if( !params.containsKey(key) ) {
            if( params.verbose ) { log.info "--Setting params.${key} to ${description}: ${def_value}" }
            params[key] = def_value
        } else if( params.verbose ) {
            log.info "--Param   params.${key} provided as: ${params[key]}"
        }
    }
    if( params.verbose ) { println "" }
}

//Ensure existence of standard paths in params.
def ensure_standard_paths(params) {
    ['main_dir', 'proj_dir', 'task_dir', 'in_dir'].each { dir ->
        if( params.containsKey(dir) ) { 
            check_full_file_path(params[dir])
        }
    }
    if( params.containsKey('task_dir') && "${params.task_dir}" != "${launchDir}") {
        eexit("Missmatched task_dir and launchDir variables: ${params.task_dir} ${launchDir}")
    }
    if( params.containsKey('in_files') ) { 
        params.in_files.sort().each{ in_file -> 
            check_full_file_path(in_file) 
        }
    }
    if( params.containsKey('aux_files') ) { 
        params.aux_files.sort().each{ aux_file -> 
            check_full_file_path(aux_file) 
        }
    }
}

// Perform all startup variable checks.
def ensure_params_paths(params) {
    fill_default_params(params, return_default_params())
    ensure_standard_paths(params)
}

// Get maximum available local cpus for processes, based on custom divisor
//def get_max_local_cpus(params, divisor) {
//    if( params.containsKey('max_local_cpus') ) {
//        return params.max_local_cpus
//    } else {
def get_max_local_cpus(divisor) {
    max_cpus = Runtime.getRuntime().availableProcessors()
    if( divisor > max_cpus ) {
        //System.err << "WARNING: Cannot Divide Available CPUS:${max_cpus}"
        //System.err << " into requested number of divisions: ${factor}\n"
        //System.err << "    Defaulting to max_local_cpus = 1\n"
        return 1
    } else {
        return max_cpus.intdiv(divisor)
    }
}

def publish_code(workflow, wf_params) {
    if( wf_params.containsKey('publish_dir') 
        && wf_params.containsKey('code_publish') 
        && wf_params['code_publish']
        && workflow.success
      ) {
        local_module_dir = file("${workflow.projectDir}/${params.local_module_dir}")
        publish_dir = file("${wf_params.publish_dir}/nf_code")
        publish_module_dir = file("${wf_params.publish_dir}/nf_code/${params.local_module_dir}")
        publish_files = []
        publish_module_files = []
        publish_files.add(file(workflow.scriptFile)) 
        workflow.configFiles.each { fn ->
            publish_files.add(file(fn))
        }
        file(workflow.launchDir).eachFileMatch(~/.*\.sbatch/) { fn ->
            publish_files.add(fn)
        }
        if( local_module_dir.exists() && local_module_dir.isDirectory() ) {
            //publish_files.add(local_module_dir) 
            local_module_dir.eachFile { fn ->
                publish_module_files.add(fn)
            }
        }
        println "\n-- Publishing main code and task config to publish directory:"
        println "     ${publish_dir}"
        publish_dir.mkdirs()
        publish_files.each { fn ->
            if( params.verbose ) {
                println "-- --${fn}"
            }
            fn.copyTo(publish_dir)    
        }
        if( publish_module_files ) {
            publish_module_dir.mkdirs()
            publish_module_files.each { fn ->
                println "-- --${fn}"
                fn.copyTo(publish_module_dir)    
            }
        }
        println ""
    }
}

// --------------- Dependency Checking / Fetching Functions: ---------------

// Get resources from params, if existing
def get_resource(name, params, type) {
    def check_str = "${name}_${type}"
    if( params.containsKey(check_str) ) {
        return params[check_str]
    } else {
        return null 
    }
}

// Get container/conda/module from params, if existing.
def get_container(name, params) { get_resource(name, params, 'container') }
def get_conda(name, params)     { get_resource(name, params, 'conda') }
def get_module(name, params)    { get_resource(name, params, 'module') }


// --------------- "params" Checking Functions: ---------------

// Check that a key is present in params, and is among allowed options.
def test_params_key_allowed(params, key, allowed_opts) {
    // Check key exists
    if( !params.containsKey(key) ) {
        log.error "Required parameter key not provided:"
        log.error "    ${key}"
        println ""
        exit 1
    }
    // if allowed is 'nonblank', check key is nonblank
    if( allowed_opts == "nonblank" ) { 
        if( params[(key)] == null || params[(key)] == "" ) {
            log.error "Value of key cannot be blank:"
            log.error "    ${key}"
            println ""
            exit 1
        }
    // If allowed_opts is provided, check that value is among allowed options.
    } else if( allowed_opts ) { 
        return_as_list(params[key]).each {
            if( !(allowed_opts.contains(opt)) ) {
                log.error "Parameter: '${key}' does not contain an allowed option:"
                log.error "  Provided: ${opt}"
                log.error "  Allowed:  ${allowed_opts}"
                println ""
                exit 1
            }
        }
    }
}             

// Check that multiple keys are present in params with allowed options
def test_params_keys_allowed(params, test_keys) {
    test_keys.each{keyopts -> test_params_key_allowed(params, *keyopts) }
}

// Check that a key is present in params.
def test_params_key(params, key) {
    test_params_key_allowed(params, key, null)
}

// Check that multiple keys are present in params.
def test_params_keys(params, test_keys) {
    test_keys.each{key -> test_params_key(params, key) }
}

// Construction of functions uncertain, currently deprecated.
/*
def test_params_file(params, in_test_file) {
    def test_file = in_test_file[0]
    if( !params.containsKey(test_file) ) {
        log.error "Required file parameter key not provided:"
        log.error "    ${test_file}"
        log.error ""
        exit 1
    }
    def this_file = file(params[test_file], checkIfExists: false)
    if( this_file instanceof List ) { this_file = this_file[0] }
    if( !this_file.exists() ) {
        log.error "Required file parameter '${test_file}' does not exist:"
        log.error check_full_file_path("${this_file}")
        exit 1
    }
} 

def test_params_files(params, test_files) {
    test_files.each{file_key -> test_params_file(params, file_key) }
}
*/


// --------------- Workflow Detail Printing Functions: ---------------

// Print command used to invoke nextflow
def print_command ( command ) {
    println ""
    log.info "Nextflow Command:"
    log.info "    ${command}"    

    // If command is extensive, print individual parameter details..
    command_list = "${command}".split()
    if( command_list.size() > 5 ) {
        log.info ""
        log.info "Nextflow Command Details:"
        message = "    " + command_list[0]
        last_command = null
        [command_list].flatten().subList(1, command_list.size()).each {
            if( it == "run" ) {
                if( last_command != "--mode" ) { log.info message ; message = "   "}
                message += " run"
            } else if( it.startsWith('-') || it.startsWith('--') ) {
                log.info message
                message = "    $it"
            } else {
                message += " $it"
            }
            last_command = it 
        }
        log.info message
    }
    println ""
}

// Print manifest details
def print_manifest_details(workflow, params) {
    def manifest_keys = [
        'name': 'Name',
        'version': 'Version',
        'description':'Description',
    ]
    if( params.verbose && manifest_keys.any{key, value -> workflow.manifest[key]} ) {
        log.info "-- Project Description:"
        manifest_keys.each{key, description -> 
            if( workflow.mainfest.containsKey(key) ) {
                log.info "${description}".padRight(params.desc_pad) + ": ${workflow.manifest[key]}"
            }
        }
    }
}

// If in_files and/or out_files provided in params, print values.
def print_use_files(params) {
    def proj_file_keys = ['in_files': 'In File(s)', 'aux_files': 'Aux. File(s)']
    
    if( proj_file_keys.any{key, value -> params.containsKey(key)} ) {
        println  ""
        log.info "Project Files:"
    }
    proj_file_keys.each { key, description ->
        if( params.containsKey(key) ) {
            pad_print_item_description(params[key], description)
        }
    }
}

// Print config file(s)
def print_config_files(workflow) {
    def config_file = ""
    config_file += '-'.multiply(params.front_pad) 
    config_file += 'NF Config Files'.padRight(params.desc_pad) 
    config_file += "${workflow.configFiles[0]}"
    log.info config_file
    workflow.configFiles.subList(1, workflow.configFiles.size()).each {
        log.info '-'.multiply(params.front_pad) + ' '.multiply(params.desc_pad) + it
    }
}

// Print workflow properties
def print_properties(workflow, wf_params) {
    [ 
        'NF Config Prof.': workflow.profile,
        'NF Script'      : workflow.scriptFile,
        'NF Launch Dir'  : workflow.launchDir,
        'NF Script'      : workflow.scriptFile,
        'NF Config Prof.': workflow.profile,
        'NF Work Dir'    : workflow.workDir,
        'NF Launch Dir'  : workflow.launchDir,
        'User'           : workflow.userName,
        'User Home Dir'  : workflow.homeDir,
        'User Main Dir'  : wf_params.main_dir,
        'Project Dir'    : wf_params.proj_dir,
        'Task Dir'       : wf_params.task_dir,
        'In Dir'         : wf_params.in_dir,
    ].each{key, value -> 
        log.info '-'.multiply(params.front_pad) + key.padRight(params.desc_pad) + value
    }
}

// Print all details about workflow
def workflow_details(workflow, params) {
    if( params.verbose ) { print_command(workflow.commandLine) }
    print_manifest_details(workflow, params)
    print_config_files(workflow)
    print_properties(workflow, params)
    print_use_files(params)
    if( params.containsKey('test_inputs') && params.test_inputs ) { eexit('Test Complete.') }
} 

// Print execution summary following pipeline
def execution_details(workflow) {
    def message = """\
        -   Nextflow Pipeline Execution Summary
        -----------------------------------
        Completed At: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        """.stripIndent().replace('\n', '\n    ')
    if( "${workflow.errorMessage}" == "SIGINT" ) {
        message += """\
            Exit Status : CANCELLED
            
            -----------------------------------
            """.stripIndent().replace('\n', '\n    ')
    } else if( workflow.exitStatus ) {
        message += """\
            Exit Status : ERROR : ${workflow.exitStatus}
            Error Text  : ${workflow.errorMessage}
            -----------------------------------
            """.stripIndent().replace('\n', '\n    ')
    } else if( !workflow.success ) {
        message += """\
            Exit Status : ERROR : ${workflow.exitStatus}
            Error Text  : ${workflow.errorMessage}
            -----------------------------------
            """.stripIndent().replace('\n', '\n    ')
    } else {
        message += """\
            Exit Status : GOOD : ${workflow.exitStatus}
                      
            -----------------------------------
            """.stripIndent().replace('\n', '\n    ')
    }
    sleep(100)
    println "\n" 
    log.info message
    println "\n\n"
    sleep(100)
}

// --------------- Workflow Process Functions: ---------------

// Begin task stdout with ouptut of relevant task settings.
def echo_process_details(task) {
    def resource_string = ""
    if( task.container ) { 
        resource_string = "Container: '${task.container}'"
    } else if( task.module ) { 
        resource_string = "Module(s): '${task.module.join(':')}'"
    } else if( task.conda ) {
        resource_string = "Conda-env: '${task.conda.toString()}'"
    } else {
        resource_string = 'Resources: None'
    }
    def label_str = ( "${task.label}" == "null" ? "" : "${return_as_list(task.label).join(', ')}" )
    def time_str  = ( "${task.time}" == "null" ? "" : "${task.time}" )
    def mem_str   = ( "${task.memory}" == "null" ? "" : "${task.memory}" )
    def queue_str = ( "${task.queue}" == "null" ? "" : "${task.queue}" )
    def scratch_str = ( "${task.scratch}" == "null" ? "" : "${task.scratch}" )
    def ret_string = """
    echo -e "\n${task.name}.${task.tag}"
    echo    "  -  Label(s):  ${label_str}"
    echo    "  -  Executor:  ${task.executor}"
    echo    "  -  CPUs:      ${task.cpus}"
    echo    "  -  Time:      ${time_str}"
    echo    "  -  Mem:       ${mem_str}"
    echo    "  -  ${resource_string}"
    echo    "  -  Queue:     ${queue_str}"
    echo -e "  -  Scratch:   ${scratch_str}\\n"
    echo    "  -  WorkDir:   \${PWD}"
    """.stripIndent()
    ret_string
}

// --------------- Channel Manipulation Functions: ---------------
//   Covnert channel input into N-item tuple.
def Cardinalize_Channel(in_channel) {
    def out_channel = in_channel.map { items ->
        if( !(items instanceof List) ) { 
            [items] + [null].multiply(params.ch_card - 1)
        } else {
            items + [null].multiply(params.ch_card - items.size())
        }
    }
    out_channel
}

//   Check Channel for Duplicate Names
def Check_Channel_Names(in_channel) {
    in_channel.map{ it -> [it[0], "yes"] }
             .groupTuple()
             .toSortedList()
             .subscribe onNext: { names ->
                 def duplicates = names.findAll { name, counts -> counts.size() > 1 }
                 if( duplicates ) { 
                     if( params.err_on_dup ) {
                         log.error "Duplicate Name Fields Detected:"
                         duplicates.each{ name, counts ->
                             log.error " -- ${name} (${counts.size()})"
                         }
                         exit 1 
                     } else {
                         log.warn "Duplicate Name Fields Detected:"
                         duplicates.each{ name, num ->
                             log.warn " -- ${name} (${counts.size()})"
                         }
                         log.warn ""
                     }
                 }
             }
}

//   Covnert channel input into N-item tuple.
def Debug_Transform_Channel(in_channel) {
    def out_channel = null
    if( params.debug ) {
        out_channel = in_channel
                       .take(2)
                       .toSortedList()
                       .flatMap{ in_list -> in_list }
                       .view(it -> "DEBUG - Limited Input: ${it[0]}")
    } else {
        out_channel = in_channel
    }
    out_channel
}

//   Perform debug-limiting, standardization, and checks on input channel before use
def Prepare_Channel ( in_channel ) {
    def temp_channel = Standardize_Channel(in_channel)
    Debug_Transform_Channel(temp_channel)
}

//   Perform standardization and checks on input channel before use
def Standardize_Channel ( in_channel ) {
    def cardinal_channel = Cardinalize_Channel(in_channel)
    Check_Channel_Names(cardinal_channel)
    cardinal_channel
}

//   Subsample fastq if debugging
def Prepare_Fastq(in_channel, Boolean pe) {
    def out_channel = null
    if( params.debug ) {
        def temp_channel_1 = Prepare_Channel(in_channel)
        def temp_channel_2 = temp_channel_1
                              .map{ in_data -> 
                                    [file(in_data[1][0]).getSimpleName()] + in_data
                              }
        
        out_channel = temp_channel_1 
                       .map{ in_data -> in_data[1] }
                       .splitFastq(file: true, by:100, compress:true, limit:100, pe:pe)
                       .map { in_fpair -> [file(in_fpair[0]).getSimpleName() - ~/\.1/, in_fpair] }
                       .join(temp_channel_2)
                       .map { in_data ->
                           [in_data[2], in_data[1]] + in_data[4..-1]
                       }
                       .toSortedList()
                       .flatMap { in_list -> in_list }
                       .view(it -> "DEBUG - Subsampled:    ${it[0]}")
    } else { 
        out_channel = Prepare_Channel(in_channel)
    }    
    out_channel
}

def Prepare_FastqPE(in_channel) {
    Prepare_Fastq(in_channel, true)
}

def Prepare_FastqSE(in_channel) {
    Prepare_Fastq(in_channel, false)
}

def Add_Parent_To_Label(in_channel) {
    def out_channel = in_channel
    .toSortedList()
    .flatMap { in_items ->
        def all_parents = ['a', 'a', 'b']
        def out_items = []
        def add_num = 0
        while( all_parents.unique(false).size() != all_parents.size() ) {
            all_parents = []
            out_items = []
            add_num += 1
            in_items.each { in_data ->
                def name = in_data[0]
                def file_names = in_data[1]
                def remainder = []
                if( in_data.size() > 2 ) {
                    remainder = in_data[2..-1]
                }
                def parent_str = ""
                def current_loop = 1
                def use_parent_path = file(file_names[0]).getParent()
                def parent_paths
                //while( "${use_parent_path}" != "/" && current_loop <= max_loops ) {
                while( "${use_parent_path}" != "/" ) {
                    add_parent = file(use_parent_path).getBaseName()
                    parent_str = "${add_parent}--${parent_str}"
                    use_parent_path = file(use_parent_path).getParent()
                    current_loop += 1
                }
                def new_name = "${parent_str}${name}"
                all_parents.add(parent_str)
                out_items.add([new_name, file_names] + remainder)
            }
        }
        out_items
    }
    out_channel
}

def Prepare_Nonunique(in_channel) {
    def temp_channel = Add_Parent_To_Label(in_channel)
    Prepare_Channel(temp_channel)
}

//   Add a sequential numeric index to first value of each item in channel
def Index_Channel_Names(in_channel) {
    in_channel.toSortedList()
             .flatMap { in_list ->
                 def out_list = []
                 def counter = 1
                 in_list.each { item ->
                     def count_str = "${counter}"
                     if( counter < 10 ) {
                         count_str = "0${count_str}"
                     }
                     use_name = "${counter}_${item[0]}"
                     out_list.add([use_name] + item[1..-1])
                     counter += 1
                 }
                 out_list
             }
}



