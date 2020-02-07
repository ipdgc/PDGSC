library(seqminer)

london_casecontrol_assoc <- Sys.getenv("london_casecontrol_assoc")
london_caseonly_assoc <- Sys.getenv("london_caseonly_assoc")
# london_casecontrol_tbi <- Sys.getenv("london_casecontrol_tbi")
# london_caseonly_tbi <- Sys.getenv("london_caseonly_tbi")

# merck_casecontrol_tbi <- Sys.getenv("merck_casecontrol_tbi")
# merck_caseonly_tbi <- Sys.getenv("merck_caseonly_tbi")
merck_casecontrol_assoc <- Sys.getenv("merck_casecontrol_assoc")
merck_caseonly_assoc <- Sys.getenv("merck_caseonly_assoc")

# download.annotation.resource("~/seqminer.annotation")

param <- makeAnnotationParameter(list(reference = "~/seqminer.annotation/hs37d5.fa", geneFile = "~/seqminer.annotation/refFlat_hg19.txt.gz", codonFile = "~/seqminer.annotation/codon.txt", priorityFile = "~/seqminer.annotation/priority.txt" ))

annotatePlain(london_casecontrol_assoc, paste("anno", london_casecontrol_assoc, sep="_"), param)
annotatePlain(london_caseonly_assoc, paste("anno", london_caseonly_assoc, sep="_"), param)
annotatePlain(merck_casecontrol_assoc, paste("anno", merck_casecontrol_assoc, sep="_"), param)
annotatePlain(merck_caseonly_assoc, paste("anno", merck_caseonly_assoc, sep="_"), param)
