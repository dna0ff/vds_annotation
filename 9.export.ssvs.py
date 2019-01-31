#!python

# Exports VDS to tsv summary stats

from hail import *
hc = HailContext(log = 'log/09.export.ssvs.log', tmp_dir = 'tmp/hail')

vds = hc.read('../cohort.annotated.nHetHom.no_samples.vds')

if not vds.was_split():
    vds = vds.split_multi()

vds = vds.annotate_variants_expr('va.allele = va.alleles[va.aIndex-1]')

kt=vds.variants_table()

kt = kt.annotate('VARIANT=v, \
CHROMOSOME=v.contig, \
START=v.start, \
REF=v.ref, \
ALT=v.alt, \
RSID=orElse(va.allele.annotation.rsid, "\\\N"), \
AC=va.allele.metrics.allele_counts.total.toInt(), \
AF=if (isMissing(va.allele.metrics.allele_frequencies.total)) "\\\N" else str(va.allele.metrics.allele_frequencies.total),\
nHomRef=va.extra.nHomRef.toInt(), \
nHet=va.extra.nHet.total.toInt(), \
nHomVar=va.extra.nHomVar.toInt(), \
TYPE=v.altAllele.category(), \
CATO=if (isMissing(va.allele.annotation.predictions.CATO.score)) "\\\N" else str(va.allele.annotation.predictions.CATO.score), \
eigen=if (isMissing(va.allele.annotation.predictions.Eigen.EigenRaw)) "\\\N" else str(va.allele.annotation.predictions.Eigen.EigenRaw), \
polyPhen=if (isMissing(va.allele.annotation.vep.predictions.PolyPhen)) "\\\N" else str(va.allele.annotation.vep.predictions.PolyPhen), \
sift=if (isMissing(va.allele.annotation.vep.predictions.SIFT)) "\\\N" else str(va.allele.annotation.vep.predictions.SIFT), \
tgpAF=if (isMissing(va.allele.annotation.freqs.tgp.AF)) "\\\N" else str(va.allele.annotation.freqs.tgp.AF), \
hrcAF=if (isMissing(va.allele.annotation.freqs.hrc.AF)) "\\\N" else str(va.allele.annotation.freqs.hrc.AF), \
gnomadAF=if (isMissing(va.allele.annotation.freqs.gnomad.AF)) "\\\N" else str(va.allele.annotation.freqs.gnomad.AF), \
gnomadAF_AFR=if (isMissing(va.allele.annotation.freqs.gnomad.AF_AFR)) "\\\N" else str(va.allele.annotation.freqs.gnomad.AF_AFR), \
gnomadAF_AMR=if (isMissing(va.allele.annotation.freqs.gnomad.AF_AMR)) "\\\N" else str(va.allele.annotation.freqs.gnomad.AF_AMR), \
gnomadAF_ASJ=if (isMissing(va.allele.annotation.freqs.gnomad.AF_ASJ)) "\\\N" else str(va.allele.annotation.freqs.gnomad.AF_ASJ), \
gnomadAF_EAS=if (isMissing(va.allele.annotation.freqs.gnomad.AF_EAS)) "\\\N" else str(va.allele.annotation.freqs.gnomad.AF_EAS), \
gnomadAF_FIN=if (isMissing(va.allele.annotation.freqs.gnomad.AF_FIN)) "\\\N" else str(va.allele.annotation.freqs.gnomad.AF_FIN), \
gnomadAF_NFE=if (isMissing(va.allele.annotation.freqs.gnomad.AF_NFE)) "\\\N" else str(va.allele.annotation.freqs.gnomad.AF_NFE), \
gnomadAF_OTH=if (isMissing(va.allele.annotation.freqs.gnomad.AF_OTH)) "\\\N" else str(va.allele.annotation.freqs.gnomad.AF_OTH), \
ensemblId=orElse(va.allele.annotation.vep.cross_references.Ensembl_gene, "\\\N"), \
consequences=orElse(va.allele.annotation.vep.consequences.consequences.mkString(","), "\\\N"), \
geneSymbol=orElse(va.allele.annotation.vep.cross_references.symbol, "\\\N"), \
clinvar=orElse(va.allele.annotation.clinvar.ClinicalSignificance.mkString(","), "\\\N")')

kt = kt.select([u'VARIANT', u'CHROMOSOME', u'START', u'REF', u'ALT', u'RSID', u'AC', u'AF', u'nHomRef', u'nHet', u'nHomVar', u'TYPE', 'CATO', 'eigen', 'sift', 'polyPhen', 'tgpAF', 'hrcAF', 'gnomadAF', 'gnomadAF_AFR', 'gnomadAF_AMR', 'gnomadAF_ASJ', 'gnomadAF_EAS', 'gnomadAF_FIN', 'gnomadAF_NFE','gnomadAF_OTH','ensemblId', 'consequences', 'geneSymbol','clinvar'])

kt.export("../export.ssvs.tsv")

