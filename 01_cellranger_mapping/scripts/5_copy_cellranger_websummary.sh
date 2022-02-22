#!/bin/bash

# Copies all the web_summary.html report files 
# from all cellranger runs into a `results/web_summary` folder
# and rename each file including the gem_id.

project_path=$(cd ../ && pwd)
cd "$project_path"
mkdir "results"
mkdir "results/web_summary"
subprojects_path="${project_path}/subprojects"
cd "$subprojects_path"
for subproject in ./*; do
	subproject_path="${subprojects_path}/${subproject}"
	jobs_path="${subproject_path}/jobs"
	cd "$jobs_path"
	for gem_id in *; do
		cp "${gem_id}/${gem_id}/outs/per_sample_outs/${gem_id}/web_summary.html" "${project_path}/results/web_summary/${subproject}__${gem_id}__web_summary.html"
	done
	cd "$subprojects_path"	
done
