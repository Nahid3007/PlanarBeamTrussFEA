# bash script to execute FEA python files   

# ENTER file path and input filename (without file ending)

input_path='./input_files/rod_truss_test/'

filename=(
# "rod_fixed_displacement"
# "rod_fixed_point_load"
"plane_truss_triangle_F"
"plane_truss_triangle_U"
# "plane_bridge_structure_LC_1"
# "plane_bridge_structure_LC_2"
# "truss_structure_w_six_members"
# "ws_11_12"
# "beam_w_square_cross_section_LC_1"
# "beam_w_square_cross_section_LC_2"
# "beam_w_square_cross_section_LC_3"
# "beam_w_square_cross_section_LC_4"
)

# run python 1D FE solver

for file in "${filename[@]}"; do
    python main_sol_linear_static.py $input_path$file.txt | tee ./results/$file.log
done

echo 'Done bash run.'   