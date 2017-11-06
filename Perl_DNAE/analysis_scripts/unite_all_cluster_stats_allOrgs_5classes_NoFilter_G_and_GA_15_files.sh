find */RC -maxdepth 2 -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_RC.txt & 
find */DNA -maxdepth 2 -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_DNA.txt & 
find */LINE -maxdepth 2 -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_LINE.txt &
find */LTR -maxdepth 2 -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_LTR.txt &
find */SINE -maxdepth 2 -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_SINE.txt &

find */SINE/results/Clusters_80_percent_length_G -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_SINE_G.txt &
find */LINE/results/Clusters_80_percent_length_G -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_LINE_G.txt &
find */LTR/results/Clusters_80_percent_length_G -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_LTR_G.txt &
find */DNA/results/Clusters_80_percent_length_G -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_DNA_G.txt &
find */RC/results/Clusters_80_percent_length_G -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_RC_G.txt &

find */SINE/results/Clusters_80_percent_length_GA -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_SINE_GA.txt &
find */LINE/results/Clusters_80_percent_length_GA -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_LINE_GA.txt &
find */LTR/results/Clusters_80_percent_length_GA -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_LTR_GA.txt &
find */DNA/results/Clusters_80_percent_length_GA -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_DNA_GA.txt &
find */RC/results/Clusters_80_percent_length_GA -name cluster_stats* -exec grep -H "" {} \; | sed 's/\//\t/g' | sed 's/\tresults.*:/\t/g' > ../DNA_editing_results/cluster_stats_all_orgs_RC_GA.txt &

 