version 1.0

workflow findMendelianVariations {
    input {
        Array[File] gvcfTrioArray
        Array[File] gvcfTrioArrayIndex
	File pedigree
        File referenceFasta
        File refIndex
        File refDict
        String region
        String dockerImage = "broadinstitute/gatk:4.1.3.0"
        String outputGVCFName
        String outputVCFName
	String outputMetrics
    }

    call combineGVCFs {
        input:
            gvcfTrioArray=gvcfTrioArray,
            gvcfTrioArrayIndex=gvcfTrioArrayIndex,
            referenceFasta=referenceFasta,
            refIndex=refIndex,
            refDict=refDict,
            region=region,
            outputGVCFName=outputGVCFName,
            dockerImage=dockerImage
    }

    call genotypeGVCFs {
        input:
            inputGVCF=combineGVCFs.outputCombinedGVCF,
            inutGVCFIndex=combineGVCFs.outputCombinedGVCFIndex,
            referenceFasta=referenceFasta,
            refIndex=refIndex,
            refDict=refDict,
            region=region,
            outputVCFName=outputVCFName,
            dockerImage=dockerImage
    }

    call findMendelianViolations {
        input:
            inputVCF=genotypeGVCFs.outputCombinedVCF,
            inutVCFIndex=genotypeGVCFs.outputCombinedVCFIndex,
            pedigree=pedigree,
            outputMetrics=outputMetrics,
            dockerImage=dockerImage
    }

    output {
        File combinedGVCF    = combineGVCFs.outputCombinedGVCF
        File combinedGVCFIdx = combineGVCFs.outputCombinedGVCFIndex

        File combinedVCF     = genotypeGVCFs.outputCombinedVCF
        File combinedVCFIdx  = genotypeGVCFs.outputCombinedVCFIndex

        File Metrics   = findMendelianViolations.outputMetricsFile
    }
}


task combineGVCFs {

    input {
        Array[File] gvcfTrioArray
        Array[File] gvcfTrioArrayIndex
        File referenceFasta
        File refIndex
        File refDict
        String region
        String outputGVCFName
        String outputGVCFIndexName = "${outputGVCFName}.idx"

        Int memSizeGB = 4
        Int diskSizeGB = 128
        String dockerImage
    }


    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        gatk CombineGVCFs \
            -R ~{referenceFasta} \
            -V ~{sep=' -V'  gvcfTrioArray} \
            -L ~{region} \
            --allow-old-rms-mapping-quality-annotation-data \
            -O ~{outputGVCFName}
    >>>

    output {
        File outputCombinedGVCF = outputGVCFName
        File outputCombinedGVCFIndex = outputGVCFIndexName
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
	preemptible: 1
    }
}

task genotypeGVCFs {

    input {
        File inputGVCF
        File inutGVCFIndex
        File referenceFasta
        File refIndex
        File refDict
        String region
        String outputVCFName
        String outputVCFIndexName = "${outputVCFName}.idx"

        Int memSizeGB = 4
        Int diskSizeGB = 128
        String dockerImage
    }


    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace


        gatk GenotypeGVCFs \
            -R ~{referenceFasta} \
            -V ~{inputGVCF} \
            -L ~{region} \
            --allow-old-rms-mapping-quality-annotation-data \
            -O ~{outputVCFName}
    >>>

    output {
        File outputCombinedVCF = outputVCFName
        File outputCombinedVCFIndex = outputVCFIndexName
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
	preemptible: 1
    }
}

task findMendelianViolations {

    input {
        File inputVCF
        File inutVCFIndex
        File pedigree
        String outputMetrics

        Int memSizeGB = 4
        Int diskSizeGB = 128
        String dockerImage
    }


    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace


        gatk FindMendelianViolations \
            -I ~{inputVCF} \
            -PED ~{pedigree} \
            -MIN_DP 10 \
            -O ~{outputMetrics}
    >>>

    output {
        File outputMetricsFile = outputMetrics
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

