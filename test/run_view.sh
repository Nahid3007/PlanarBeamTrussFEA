# bash script to execute FEA python files   

# ENTER input file name (without file ending)

input_file_name=(
# "rod_fixed_displacement"
# "rod_fixed_point_load"
# "plane_truss_triangle_F"
# "plane_truss_triangle_U"
# "plane_bridge_structure_LC_1"
# "plane_bridge_structure_LC_2"
# "truss_structure_w_six_members"
# "ws_11_12"
"beam_w_square_cross_section_LC_1"
"beam_w_square_cross_section_LC_2"
"beam_w_square_cross_section_LC_3"
"beam_w_square_cross_section_LC_4"
)

# run preview input file

for infile in "${input_file_name[@]}"; do
    python view_input_file.py $infile.txt
done

echo 'Done bash run.'   