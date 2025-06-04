shopt -s extglob

rm -rf report.html logs output/

if [ -d "work" ]; then
	cd work
	rm -rf  !(apptainer|conda)
else
  echo "Work directory does not exist."
fi
