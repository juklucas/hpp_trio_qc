# hpp_trio_qc

This repository holds a workflow to check for Mendelian violations in trio GVCF data from the [1000G project's 30X dataset](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38). The trio is tested in this workflow using the GATK/Picard tool [FindMendelianViolations](https://gatk.broadinstitute.org/hc/en-us/articles/360040098772-FindMendelianViolations-Picard-)

------------------

### Steps taken by the WDL

**Requires GVCF inputs for trio (along with valid pedigree file)**
* CombineGVCFs
* GenotypeGVCFs
* FindMendelianViolations
