# bash script to execute FEA python files   

# input file name (without file ending)

inputfile_name=plane_bridge_structure_LC_2


# run python 1D FE solver

python run_linear_static.py $inputfile_name.txt | tee ./results/$inputfile_name.log

echo
echo 'Done bash script.'   