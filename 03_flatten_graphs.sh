cd /lustre/scratch/jmanthey/18_contopus_sv/00_data


awk '/^S/{print ">"$2;print $3}' C_sordidulus__KU39667.asm.bp.hap1.p_ctg.gfa > ../01_haplotypes/C_sordidulus__KU39667__hap1.fa

awk '/^S/{print ">"$2;print $3}' C_sordidulus__KU39667.asm.bp.hap2.p_ctg.gfa > ../01_haplotypes/C_sordidulus__KU39667__hap2.fa


awk '/^S/{print ">"$2;print $3}' C_virens__KU39615.asm.bp.hap1.p_ctg.gfa > ../01_haplotypes/C_virens__KU39615__hap1.fa

awk '/^S/{print ">"$2;print $3}' C_virens__KU39615.asm.bp.hap2.p_ctg.gfa > ../01_haplotypes/C_virens__KU39615__hap2.fa


