EDTA_DIRECTORY = "/mnt/shared/scratch/msmith/apps/EDTA"

process Helixer {
    container 'docker://gglyptodon/helixer-docker:helixer_v0.3.0_cuda_11.2.0-cudnn8'
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
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    queue 'short'
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

process GetBed {
    container 'https://depot.galaxyproject.org/singularity/agat:1.3.3--pl5321hdfd78af_0'
    publishDir 'output'
    cpus 1
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    queue 'short'
    script:
    input:
    tuple val(genome), path(gff)
    output:
    tuple val(genome), path("${genome}.bed")
    """
    agat_convert_sp_gff2bed.pl \
      --gff ${gff} \
      -o temp.bed

    awk -F'\t' '{print \$1, \$2, \$3, \$4}' temp.bed > ${genome}.bed
    """
}

process Resistify {
    conda 'resistify'
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
    resistify --threads 4 ${peptide} ${genome}
    """
}

process Orthofinder {
    container 'docker://davidemms/orthofinder:2.5.2'
    publishDir 'output'
    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    queue 'long'
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
    conda 'edta.yml'
    publishDir 'output', mode: 'copy'
    cpus 16
    memory { 64.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    queue 'long'
    input:
    tuple val(genome), path(fasta)
    output:
    tuple val(genome), path("${fasta}"), path("${fasta}.mod"), path("${fasta}.mod_divergence_plot.pdf"), path("${fasta}.mod.EDTA.anno"), path("${fasta}.mod.EDTA.combine"), path("${fasta}.mod.EDTA.final"), path("${fasta}.mod.EDTA.intact.fa"), path("${fasta}.mod.EDTA.intact.gff3"), path("${fasta}.mod.EDTA.raw"), path("${fasta}.mod.EDTA.TEanno.gff3"), path("${fasta}.mod.EDTA.TEanno.sum"), path("${fasta}.mod.EDTA.TElib.fa"), path("${fasta}.mod.MAKER.masked"), path("${fasta}.mod.RM2.raw.fa")
    script:
    """
    ${EDTA_DIRECTORY}/EDTA.pl \
      --genome ${fasta} \
      --threads 16 \
      --anno 1
    """
}

//process BedtoolsIntersect {
//    container 'quay.io/biocontainers/bedtools:2.30.0--hdb8b3b0_0'
//    cpus 1
//    memory '2 GB'
//    queue 'short'
//    input:
//    path bed1
//    path bed2
//    output:
//    path "${bed1.baseName}.intersect.${bed2.baseName}.bed"
//    script:
//    """
//    bedtools intersect -a ${bed1} -b ${bed2} > ${bed1.baseName}.intersect.${bed2.baseName}.bed
//    """
//}

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
    helixerBed = GetBed(helixer)

    Resistify(peptide)

    orthofinder = peptide
      .map { genome, pep -> pep }
      .collect()
      | Orthofinder

    edta = Edta(genomes)

    //bedtoolsIntersect = BedtoolsIntersect(helixerBed, edtaBed)
}