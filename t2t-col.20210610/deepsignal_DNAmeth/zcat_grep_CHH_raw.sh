#!/bin/bash

zcat Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw.tsv.gz | grep -P "^Chr1\t" - > Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw_Chr1.tsv

zcat Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw.tsv.gz | grep -P "^Chr2\t" - > Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw_Chr2.tsv

zcat Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw.tsv.gz | grep -P "^Chr3\t" - > Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw_Chr3.tsv

zcat Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw.tsv.gz | grep -P "^Chr4\t" - > Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw_Chr4.tsv

zcat Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw.tsv.gz | grep -P "^Chr5\t" - > Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CHH_raw_Chr5.tsv
