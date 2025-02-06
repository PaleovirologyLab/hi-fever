shopt -s extglob

rm -rf .nextflow.log* report.html logs output/

cd work

rm -rf !(apptainer)
