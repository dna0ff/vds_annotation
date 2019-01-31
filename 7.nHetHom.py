#!python
# to fix org.codehaus.janino.JaninoRuntimeException:
# Code of method "(Lorg/apache/spark/sql/catalyst/expressions/GeneratedClass;[Ljava/lang/Object;)V" of class "org.apache.spark.sql.catalyst.expressions.GeneratedClass$SpecificUnsafeProjection" grows beyond 64 KB --> upgrade to Spark 2.1.0+

# Annotates the existing VDS with nHom and nHet

from hail import *
hc = HailContext(log = 'log/07.nHetHom.log', tmp_dir = 'tmp/hail')

vds = hc.read('../cohort.annotated.vds')
vds = vds.split_multi()

vds = vds.annotate_variants_expr([
        'va.allele = va.alleles[va.aIndex-1]',
        'va.extra.wasSplit = va.wasSplit',
        'va.extra.nHomVarFemale = gs.filter(s => sa.pheno.isFemale).filter(g => g.isHomVar()).count()',
        'va.extra.nHomVarMale = gs.filter(s => !sa.pheno.isFemale).filter(g => g.isHomVar()).count()',
        'va.extra.nHomRefFemale = gs.filter(s => sa.pheno.isFemale).filter(g => g.isHomRef()).count()',
        'va.extra.nHomRefMale = gs.filter(s => !sa.pheno.isFemale).filter(g => g.isHomRef()).count()'
    ]).annotate_variants_expr([
        'va.extra.nHomVar = va.extra.nHomVarMale + va.extra.nHomVarFemale',
        'va.extra.nHomRef = va.extra.nHomRefMale + va.extra.nHomRefFemale'
    ])

vds = vds.annotate_variants_expr('''
    va.extra.nHet = {
        total: 
            if (v.inXNonPar())
                (2 * (va.allele.metrics.allele_counts.male - va.extra.nHomVarMale)) + (va.allele.metrics.allele_counts.female - (2 * va.extra.nHomVarFemale))
            else if (v.inYNonPar())
                0
            else
                va.allele.metrics.allele_counts.total - (2 * va.extra.nHomVar)   
            }
    ''')

# vds = vds.filter_variants_expr('va.allele.metrics.allele_counts.total.toInt() <= 0', keep=False)

vds.write('../cohort.annotated.nHetHom.vds')
