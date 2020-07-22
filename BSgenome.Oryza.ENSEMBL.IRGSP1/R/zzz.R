###
###

.pkgname <- "BSgenome.Oryza.ENSEMBL.IRGSP1"

.seqnames <- NULL

.circ_seqs <- NULL

.mseqnames <- NULL

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Oryza sataiva japonica",
        common_name="Oryza",
        provider="ENSEMBL",
        provider_version="IRGSP1.0",
        release_date="November, 2018",
        release_name="Oryza sataiva japonica IRGSP1.0",
        source_url="ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/oryza_sativa/dna/",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "Oryza"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

