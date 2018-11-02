#ifndef BZ_SITE__H
#define BZ_SITE__H

blast_table_t *siteblast_table_new(SEQ *seq, int max_seq_length, bz_flags_t *bz_flags);

msp_table_t *siteblast_search(SEQ *seq1, SEQ *seq2, blast_table_t *bt, enum blast_table_type btt, 
			      msp_table_t *msp,
			      ss_t ss, int X, int K, int P, int T);
msp_table_t *process_siteblast_hits(SEQ *seq1, SEQ *seq2, blast_table_t *bt, 
				    enum blast_table_type btt,
				    msp_table_t *msp_tab,
				    ss_t ss, int X, int K, int P, int T);
void correct_bt_for_2nd_search(int length , blast_table_t *bt,enum blast_table_type btt);

void siteblast_table_i_free(blast_table_t *bt);
void siteblast_table_I_free(blast_table_t *bt);
void siteblast_table_p_free(blast_table_t *bt);

#endif
