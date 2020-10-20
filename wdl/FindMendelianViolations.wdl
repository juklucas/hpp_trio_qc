version 1.0

workflow findMendelienVariations {
    input {
        Array[File] gvcfTrioArray
        Array[File] gvcfTrioArrayIndex
        File referenceFasta
        File refIndex
        File refDict
        String region
        String dockerImage = "broadinstitute/gatk:4.1.3.0"
        String ouputGVCFName
        String ouputVCFName
    }

    call combineGVCFs {
        input:
            gvcfTrioArray=gvcfTrioArray,
            gvcfTrioArrayIndex=gvcfTrioArrayIndex,
            referenceFasta=referenceFasta,
            refIndex=refIndex,
            refDict=refDict,
            region=region,
            ouputGVCFName=ouputGVCFName,
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
            ouputVCFName=ouputVCFName,
            dockerImage=dockerImage
    }

    output {
        File combinedGVCF    = combineGVCFs.outputCombinedGVCF
        File combinedGVCFIdx = combineGVCFs.outputCombinedGVCFIndex

        File combinedVCF     = genotypeGVCFs.outputCombinedVCF
        File combinedVCFIdx  = genotypeGVCFs.outputCombinedVCFIndex
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
        String ouputGVCFName
        String ouputGVCFIndexName = "${ouputGVCFName}.idx"

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
            -O ~{ouputGVCFName}
    >>>

    output {
        File outputCombinedGVCF = ouputGVCFName
        File outputCombinedGVCFIndex = ouputGVCFIndexName
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
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
        String ouputVCFName
        String ouputVCFIndexName = "${ouputVCFName}.idx"

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
            -O ~{ouputVCFName}
    >>>

    output {
        File outputCombinedVCF = ouputVCFName
        File outputCombinedVCFIndex = ouputVCFIndexName
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}






