#!/bin/bash

# Submit one cellranger job per library,
# where cellranger maps the fastq files to a reference
# and calls the feature-barcode matrices.

project_path=$(cd ../ && pwd)
subprojects_path="${project_path}/subprojects"
cd "$subprojects_path"
for subproject in ./*; do
	subproject_path="${subprojects_path}/${subproject}"
	jobs_path="${subproject_path}/jobs"
	cd "$jobs_path"
	for gem_id in ./*; do
		cd "$gem_id"
		if [[ ! -d "${gem_id}" ]]
		then
			sbatch "${gem_id}.cmd"
		fi
		cd "$jobs_path"
	done
	cd "$subprojects_path"
done
