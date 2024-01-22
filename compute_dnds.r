source('KaKs_DnDs/script/kaks.R')

# compute ka/ks for a given alignment
p1 <- kaks("core_gene_alignment_test.fasta", genetic_code_path = 'KaKs_DnDs/data/genetic_code.txt', sn_table_path = 'KaKs_DnDs/data/SN_codonTable.txt')

# compare to 1
dnds_all <- p1$KaKs
wilcox.test(dnds_all, mu=1)