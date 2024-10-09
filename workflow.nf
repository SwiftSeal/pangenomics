EDTA_DIRECTORY = "/mnt/shared/scratch/msmith/apps/EDTA"

process Helixer {
    container 'docker://gglyptodon/helixer-docker:helixer_v0.3.0_cuda_11.2.0-cudnn8'
    publishDir 'output'
    cpus 4
    memory '32 GB'
    queue 'gpu'
    clusterOptions '--gpus=1'
    containerOptions '--nv'
    input:
    tuple val(genome), path(fasta)
    path model
    output:
    tuple val(genome), path("${genome}.gff")
    script:
    """
    Helixer.py \
      --species ${genome} \
      --fasta-path ${fasta} \
      --gff-output-path ${genome}.gff \
      --model-filepath ${model} \
    """
}

process GetPeptide {
    container 'https://depot.galaxyproject.org/singularity/agat:1.3.3--pl5321hdfd78af_0'
    publishDir 'output'
    cpus 1
    time '1h'
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    script:
    input:
    tuple val(genome), path(gff), path(fasta)
    output:
    tuple val(genome), path("${genome}.pep")
    """
    agat_sp_extract_sequences.pl \
      --cfs \
      -g ${gff} \
      -f ${fasta} \
      -p -o ${genome}.pep
    """
}

process Resistify {
    container 'docker://quay.io/biocontainers/resistify:0.3.0--pyhdfd78af_0'
    publishDir 'output'
    cpus 4
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    input:
    tuple val(genome), path(peptide)
    output:
    path "${genome}"
    script:
    """
    resistify --threads $task.cpus ${peptide} ${genome}
    """
}

process Orthofinder {
    container 'docker://davidemms/orthofinder:2.5.2'
    publishDir 'output'
    time '6d'
    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    input:
    path peptides
    output:
    path "orthofinder"
    script:
    """
    mkdir -p peptides
    cp *.pep peptides/
    orthofinder -f peptides -t 16 -o orthofinder
    """
}

process Edta {
    container 'docker://quay.io/biocontainers/edta:2.0.0--hdfd78af_0'
    publishDir 'output', mode: 'copy'
    cpus 16
    time '6d'
    memory { 128.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    input:
    tuple val(genome), path(fasta)
    output:
    tuple val(genome), path("${fasta}"), path("${fasta}.mod"), path("${fasta}.mod.EDTA.anno"), path("${fasta}.mod.EDTA.combine"), path("${fasta}.mod.EDTA.final"), path("${fasta}.mod.EDTA.intact.gff3"), path("${fasta}.mod.EDTA.raw"), path("${fasta}.mod.EDTA.TEanno.gff3"), path("${fasta}.mod.EDTA.TEanno.sum"), path("${fasta}.mod.EDTA.TElib.fa"), path("${fasta}.mod.MAKER.masked")
    script:
    """
    unset which # really???
    EDTA.pl \
      --genome ${fasta} \
      --threads $task.cpus \
      --anno 1
    """
}

workflow {
    model = file("/mnt/shared/home/msmith/.local/share/Helixer/models/land_plant/land_plant_v0.3_a_0080.h5")
    genomes = Channel.of(
      "annuum",
      "bulbocastanum",
      "candolleanum",
      "chilense",
      "chmielewskii",
      "corneliomulleri",
      "etuberosum",
      "galapagense",
      "habrochaites",
      "lycopersicoides",
      "lycopersicum",
      "neorickii",
      "peruvianum",
      "pimpinellifolium",
      "tuberosum_phureja_E4_63",
      "tuberosum_phureja_E86_69",
      "tuberosum_stenotomum_A6_26",
      "tuberosum_stenotomum_PG6359",
      "tuberosum_tuberosum_RH10_15",
      "tuberosum_tuberosum_RH",
      "verrucosum",
    ) map { genome -> [genome, file("genomes/${genome}.fa")] }

    helixer = Helixer(genomes, model)

    peptide = GetPeptide(helixer.combine(genomes, by: 0))

    Resistify(peptide)

    orthofinder = peptide
      .map { genome, pep -> pep }
      .collect()
      | Orthofinder

    edta = Edta(genomes)
}
