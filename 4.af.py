#!python
import hail

hc = hail.HailContext(log = 'log/04af.log', tmp_dir = 'tmp/hail')

# Remove related and 45-and-up cancer samples.  Also drop variants 
# with no sample support, and calculate sex-aware variant frequencies.
#
# We do not drop alleles without sample support, for rather involved
# reasons around the immutability of the Hail TGenotype type not 
# interacting well with filter_alleles.  Essentially, following a 
# filter_alleles stage, the v.oneHotGenotypes and v.oneHotAlleles
# methods do not work properly for variants that have had alleles
# dropped; this means that downstream gtCounts and alleleCounts
# calculations fail.  Performing allelic filtering *after* gtCounts
# and alleleCounts calculation doesn't help: there is no easy way
# to filter gtCounts to match the new alleles, and in rare cases
# to do with variants on the Y chromosome only in females, and the
# genotype-changing nature of the subset=True mode of filter_alleles,
# would result in the alleleCounts potentially not agreeing with the
# genotypes.
#
# There may be a way to resolve this through recreating the TGenotypes
# during allele subsetting, or performing the subsetting 'manually' in
# an annotate_genotypes_expr stage, but I'm not sure how to achieve this.
#
# Note about va.info.gtCounts and va.info.alleleCounts:
#   * va.info.gtCounts is raw counts split by male and female; no attempt
#     is made to correct for the different interpretation of the gonosomes
#     between sexes.  Downstream logic should take this into account if
#     needed.
#   * va.info.alleleCounts *does* correct for gonosomes.  Specifically,
#     allele counts for Y loci are set to zero for females, and allele counts
#     for each male sample are set to sum to one on both X and Y loci.
#     As it was judged too complicated to identify for male heterozygous
#     gonosome genotypes the most likely allele present, the weight is split
#     across the two component alleles.  This means that total allele counts
#     are on average correct for males on the X and Y chromosomes, but might
#     take on non-integer values in some cases.


vds_mgrb_filtered = hc.read('../cohort.minrep.vds')

vds_mgrb_filtered = (vds_mgrb_filtered
    .annotate_variants_expr('''
        va.gatk = va.info, 
        va.info = {NS:
            if (v.inYNonPar())
                gs.filter(s => !sa.pheno.isFemale).filter(g => g.isCalled()).count()
            else
                gs.filter(g => g.isCalled()).count()
        }''')
    .filter_variants_expr('va.info.NS > 0')
    .annotate_variants_expr('''
        va.info.gtCounts = {
            female: gs.filter(s => sa.pheno.isFemale).map(g => if (g.isCalled()) g.oneHotGenotype(v) else range(v.nGenotypes()).map(x => 0)).sum(),
            male: gs.filter(s => !sa.pheno.isFemale).map(g => if (g.isCalled()) g.oneHotGenotype(v) else range(v.nGenotypes()).map(x => 0)).sum()},
        va.info.alleleCounts = {
            female:
                if (v.inYNonPar())
                    gs.filter(s => sa.pheno.isFemale).map(g => range(v.nAlleles()).map(x => 0.0)).sum()
                else
                    gs.filter(s => sa.pheno.isFemale).map(g => if (g.isCalled()) g.oneHotAlleles(v)*1.0 else range(v.nAlleles()).map(x => 0.0)).sum(),

            male:
                if (v.inXNonPar() || v.inYNonPar())
                    gs.filter(s => !sa.pheno.isFemale).map(g => if (g.isCalled()) g.oneHotAlleles(v)*0.5 else range(v.nAlleles()).map(x => 0.0)).sum()
                else
                    gs.filter(s => !sa.pheno.isFemale).map(g => if (g.isCalled()) g.oneHotAlleles(v)*1.0 else range(v.nAlleles()).map(x => 0.0)).sum()}''')
    .annotate_variants_expr('''
        va.info.gtCounts.total = va.info.gtCounts.female + va.info.gtCounts.male,
        va.info.alleleCounts.total = va.info.alleleCounts.female + va.info.alleleCounts.male''')
    .annotate_variants_expr('''
        va.info.AC = va.info.alleleCounts.total,
        va.info.AF = va.info.alleleCounts.total / va.info.alleleCounts.total.sum()''')
)

print("After variant removal:")
print(vds_mgrb_filtered.count())

vds_mgrb_filtered.write('tmp/04.cohort.minrep.af.vds')
