#! python
import hail
from hail.expr import TString, TBoolean, TFloat, TInt

hc = hail.HailContext(log = 'log/1.load.log', tmp_dir = 'tmp/hail')

vds = hc.import_vcf('../cohort.vcf')

sample_table = (hc
    .import_table('../metadata/cohort_metadata.csv', delimiter=',', types={
        'sampleID': TString(),
        'cohort': TString(),
        'YOB':TInt(),
        'isFemale': TBoolean()
    })
    .key_by('sampleID')
)

vds = vds.annotate_samples_table(sample_table, root='sa.pheno')

# Note the use of min_rep only here: vt norm is not performed (as it was for phase 1).
# In phase 1 we observed that vt norm almost never changed GATK's representation (rate
# < 1e-5), so this extra step was considered not worthwhile.  As it is, in the presence
# of complex and multi-allelic variants, simple variant normalization is not especially
# useful.
vds.min_rep().repartition(500).write('../cohort.minrep.vds')
