# bash script to execute FEA python files   

# ENTER file path and input filename (without file ending)

input_path_1='./input_files/truss_test/'
input_path_2='./input_files/beam_test/'
input_path_3='./input_files/truss_beam_test/'

filename=(
# Truss Models
"$input_path_1 rod_fixed_displacement"
"$input_path_1 rod_fixed_point_load"
"$input_path_1 plane_truss_triangle_F"
"$input_path_1 plane_truss_triangle_U"
"$input_path_1 plane_bridge_structure_LC_1"
"$input_path_1 plane_bridge_structure_LC_2"
"$input_path_1 truss_structure_w_six_members"
# Beam Models
"$input_path_2 beam_w_square_cross_section_LC_1"
"$input_path_2 beam_w_square_cross_section_LC_2"
"$input_path_2 beam_w_square_cross_section_LC_3"
"$input_path_2 beam_w_square_cross_section_LC_4"
"$input_path_2 beam_w_distributed_load"
"$input_path_2 portal_frame_w_distr_load"
"$input_path_2 plane_bridge_structure"
# Truss/Beam Models
"$input_path_3 ws_11_12"
"$input_path_3 ss_09"
"$input_path_3 ss_10"
"$input_path_3 ws_09_10"
)

for line in "${filename[@]}"; do
  read -r path file <<< $line
  # echo $path $file
  python3 test_txt2vtk.py --infile $path$file/$file.txt --outpath $path$file/

  python3 test_solve_LinearStatic.py --infile $path$file/$file.txt --outpath $path$file/

  python3 test_vtkpost.py --vtuFile $path$file/ELEMENTS.vtu --resultsFile $path$file/$file.h5 --outputFile $path$file/ELEMENTS_WITH_RESULTS.vtu
done

echo 'Done bash run.'