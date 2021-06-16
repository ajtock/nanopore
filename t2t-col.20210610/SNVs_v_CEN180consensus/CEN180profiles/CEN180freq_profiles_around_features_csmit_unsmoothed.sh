#!/bin/bash

csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr1 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;
csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr2 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;
csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr3 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;
csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr4 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;
csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr5 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;

csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr1 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;
csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr2 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;
csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr3 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;
csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr4 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;
csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr5 180 2000 2kb 10 10bp genomewide CEN180" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide CENsoloLTR" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr2 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr2 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr2 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr3 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr3 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr3 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr4 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr4 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr4 2000 2000 2kb 10 10bp genomewide CENsoloLTR" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr4 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr5 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr5 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr5 2000 2000 2kb 10 10bp genomewide CENsoloLTR" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./CEN180freq_profiles_around_features.R Chr5 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide CENsoloLTR" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr1 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr2 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr2 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr2 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr3 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr3 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr3 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr4 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr4 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr4 2000 2000 2kb 10 10bp genomewide CENsoloLTR" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr4 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;

#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr5 2000 2000 2kb 10 10bp genomewide CENgap" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr5 2000 2000 2kb 10 10bp genomewide CENAthila" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr5 2000 2000 2kb 10 10bp genomewide CENsoloLTR" & sleep 10;
#csmit -m 130G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./SNVs_v_CEN180consensus_profiles_around_features.R Chr5 2000 2000 2kb 10 10bp genomewide nonCENAthila" & sleep 10;
