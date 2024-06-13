process Helixer {
    container 'docker://gglyptodon/helixer-docker:helixer_v0.3.0_cuda_11.2.0-cudnn8'
    cpus 4
    memory '16 GB'
    queue 'gpu'
    clusterOptions '--gpus=1'
    containerOptions '--nv'
    input:
    path genome
    path model
    output:
    path "${genome.baseName}.gff"
    script:
    """
    Helixer.py \
      --species ${genome.baseName} \
      --fasta-path ${genome} \
      --gff-output-path ${genome.baseName}.gff \
      --model-filepath ${model} \
    """
}

process GetPeptide {
    container 'https://depot.galaxyproject.org/singularity/agat:1.3.3--pl5321hdfd78af_0'
    publishDir 'output'
    cpus 1
    memory '2 GB'
    queue 'short'
    script:
    input:
    path gff
    output:
    path "${gff.baseName}.pep"
    path "${gff.baseName}.bed"
    """
    agat_sp_extract_sequences.pl \
      --cfs \
      -g ${gff} \
      -f ${genome} \
      -p -o ${genome.baseName}.pep

    agat_convert_sp_gff2bed.pl \
      --gff ${gff} \
      -o ${gff.baseName}.bed

    awk -i inplace '{{OFS="\t"; print $1, $2, $3, $4}}' ${gff.baseName}.bed
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
    path peptide
    output:
    path "${peptide.baseName}"
    script:
    """
    resistify ${peptide} ${peptide.baseName}
    """
}

process Orthofinder {
    container 'docker://davidemms/orthofinder:2.5.2'
    publishDir 'output'
    cpus 16
    memory '16 GB'
    queue 'medium'
    input:
    path genomes
    output:
    path "orthofinder"
    script:
    """
    orthofinder -f ${genomes} -t 16 -o orhtofinder
    """
}

process Edta {
    scratch true
    cpus 8
    memory '32 GB'
    queue 'long'
    input:
    path genome
    output:
    path "${genome.baseName}.TEanno.gff3"
    path "${genome.baseName}.TElib.fa"
    script:
    """
    source activate edta
    EDTA.pl \
      --genome ${genome} \
      --species ${genome.baseName} \
      --threads 8 \
      --anno 1
    """
}

process EdtaBed {
    container 'https://depot.galaxyproject.org/singularity/agat:1.3.3--pl5321hdfd78af_0'
    cpus 1
    memory '2 GB'
    queue 'short'
    input:
    path gff
    output:
    path "${gff.baseName}.bed"
    script:
    """
    agat_convert_sp_gff2bed.pl \
      --gff ${gff} \
      -o ${gff.baseName}.bed
    """
}

process BedtoolsIntersect {
    container 'quay.io/biocontainers/bedtools:2.30.0--hdb8b3b0_0'
    cpus 1
    memory '2 GB'
    queue 'short'
    input:
    path bed1
    path bed2
    output:
    path "${bed1.baseName}.intersect.${bed2.baseName}.bed"
    script:
    """
    bedtools intersect -a ${bed1} -b ${bed2} > ${bed1.baseName}.intersect.${bed2.baseName}.bed
    """
}

workflow {
    model = file("/mnt/shared/home/msmith/.local/share/Helixer/models/land_plant/land_plant_v0.3_a_0080.h5")
    genomes = Channel.fromPath("genomes/*.fna")

    gff = Helixer(genomes, model)

    peptide = GetPeptide(gff)

    Resistify(peptide)

    orthofinder = Orthofinder(peptide.collect())

    edta = Edta(genomes)

    edtaBed = EdtaBed(edta)

    bedtoolsIntersect = BedtoolsIntersect(edtaBed, orthofinder.out)
}