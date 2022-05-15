#!/bin/bash

./feature_among_read_variation_scoring_func_DMRs_sets.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'cmt3_BSseq_Rep1_and_met1_BSseq_Rep1_hypoCHG' 'bodies'
./feature_among_read_variation_scoring_func_DMRs_sets.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'cmt3_BSseq_Rep1_and_met1_BSseq_Rep1_hypoCHG' 'bodies'

./feature_among_read_variation_scoring_func_DMRs_sets.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'cmt3_BSseq_Rep1_not_met1_BSseq_Rep1_hypoCHG' 'bodies'
./feature_among_read_variation_scoring_func_DMRs_sets.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'cmt3_BSseq_Rep1_not_met1_BSseq_Rep1_hypoCHG' 'bodies'

./feature_among_read_variation_scoring_func_DMRs_sets.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'met1_BSseq_Rep1_and_cmt3_BSseq_Rep1_hypoCG' 'bodies'
./feature_among_read_variation_scoring_func_DMRs_sets.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'met1_BSseq_Rep1_and_cmt3_BSseq_Rep1_hypoCG' 'bodies'

./feature_among_read_variation_scoring_func_DMRs_sets.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CpG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'met1_BSseq_Rep1_not_cmt3_BSseq_Rep1_hypoCG' 'bodies'
./feature_among_read_variation_scoring_func_DMRs_sets.R Col_0_deepsignalDNAmeth_30kb_90pc t2t-col.20210610 CHG 0.50 1.00 'Chr1,Chr2,Chr3,Chr4,Chr5' 'met1_BSseq_Rep1_not_cmt3_BSseq_Rep1_hypoCG' 'bodies'

