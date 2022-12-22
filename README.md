# Overview

Simple python scripts to calculate one dimensional planar truss and beam structures based on the finite element approach [[1]](#1).

Five main python scripts have been developped at the moment:
- `fem_classes.py`
- `fem_functions`
- `parse_input_file_linear_static.py`
- `visualize_input_file.py`
- `main_sol_linear_static.py`

<br />

# Input File Strucutre

All input files have to be created manually and must be stored in the 
```
./input_files
```

directory.

### **Node Definition**

The nodal inputs are defined after using the `*node` keyword. 

```
*node
node_id_1, x_1, y_1
node_id_2, x_2, y_2
...
```

### **Element Definition**

The element inputs are defined after using the `*element` keyword. 

- Rod element types are defined as followed:

```
*element, type=rod
element_id_1, node_id_1, node_id_2
...
```

- Beam element types are defined as followed:

```
*element, type=beam
element_id_2, node_id_2, node_id_3
...
```
### **Property Definition**

The element property inputs are defined after using the `*property` keyword. 

- Rod properties are defined as followed:

```
*property rod
first_element_id, last_element_id, E, A
```

- Beam properties are defined as followed:

```
*property beam
first_element_id, last_element_id, E, A, I, h_max
```

So for example if 
- element id 1, 2, 4 

are defined as rod elements and 
- element id 3, 5, 6 

are defined as beam elements and E, A, I and h_max remains constant through out all rod/beam elements the property block looks as follwed:

```
*property rod
1, 2, 70000., 40.
4, 4, 70000., 40.
*property beam
3, 3, 210000., 100., 3333.33, 10.
5, 6, 210000., 100., 3333.33, 10.
```

### **Loads Definition**

The load inputs are defined after using the `*load` keyword.

```
*load
node_id, dof, value 
```

So a force applied onto node id 2 in Y direction with a value of -10 (N) and a moment applied onto node id 4 around the plane normal with a value of 50 (Nmm) is defined as followed:

```
*load
2, 2, -10.
4, 3, 50. 
```

### **Constraints Definition**

The constraints inputs are defined after using the `*boundary` keyword.

```
*boundary
node_id, first_dof, last_dof, value 
```
- Fixed constraints

A fixed dof has a value of 0. So for example an encastred boundary condition for a beam node (e.g. id 2) is defined as followed:

```
*boundary
2, 1, 3, 0. 
```

- Enforced displacement

Enforced displacements are also defined within the constraints definition. If a node (e.g. id 3) is remoted with a displacement value of -1 (mm) in X direction the input is defined as followed:

```
*boundary
3, 1, 1, -1. 
```

<br />

# Visualize the Input File

To visualize the finite element model defined by the input file the `run_view.sh` bash script can be used.

The bash script - as an example - looks as followed:

```
# bash script to execute FEA python files   

# ENTER file path and input filename (without file ending)

input_path='./input_files/rod_truss_test/'

filename=(
"rod_fixed_displacement"
# "rod_fixed_point_load"
)

# run preview input file

for file in "${filename[@]}"; do
    python visualize_input_file.py $input_path$file.txt
done

echo 'Done bash run.' 
```

As mentioned all input files have to be stores in the `./input_files` directory. When input files are stored in subfolder as for example in `rod_truss_test` in the example bash file, this folder has to be explictly specified in the `input_path` variable.

The `filename` array contains all input filenames without the file ending. Multiple filename can be specified by entering it row by row (*do not forget the quotations*). If the a row starts with a hash (`#`) sign this line will be igored.

To execute the bash script, run in the terminal the following command

```
./run_view.sh
```

and it will automatically run the `view_input_file.py` script. The output `*_FE_Model.png` file will be stored in the specified input file directory.

<br />

Example output result of the input file visualization:

![Alt text](/test/input_files/rod_truss_test/truss_structure_w_six_members_FE_Model.png "truss_structure_w_six_members")

<br />

# Solve the FE Model

To solve the FE model input file the `run_solve.sh` bash script can be used.

The bash script - as an example - looks as followed:

```
# bash script to execute FEA python files   

# ENTER file path and input filename (without file ending)

input_path='./input_files/rod_truss_test/'

filename=(
"rod_fixed_displacement"
# "rod_fixed_point_load"
)

# run python 1D FE solver

for file in "${filename[@]}"; do
    python main_sol_linear_static.py $input_path$file.txt | tee ./results/$file.log
done

echo 'Done bash run.'    
```

The bash script inputs are similar to `run_view.sh` and is described above.

To execute the bash script, run in the terminal the following command

```
./run_solve.sh
```

and it will automatically run the `run_linear_static.py` script. 

Outputs of the bash script run are:
- a log file (`*.log`)
- a result file (`res_input_file_name.txt`)

and they will be stored in
```
./results
```

The terminal outputs will be stored in a `*.log` file, and can be used for example for debugging purposes.

<br />


Example output log file result:

```
[INF] Run Linear Static Solver
[INF] Parsing input file: truss_structure_w_six_members.txt
      Number of nodes: 5
      Number of elements: 6
      Number of SPCs: 4
      Number of loads: 1
[INF] Element properties (eid, n1, n2, l, alpha, E, A, I)
      1, 1, 2, 2.00, 0.00, 30.0, 5.0
      2, 3, 4, 2.00, 0.00, 30.0, 5.0
      3, 1, 4, 2.83, 315.00, 30.0, 5.0
      4, 2, 4, 2.00, 270.00, 30.0, 5.0
      5, 4, 5, 2.00, 0.00, 30.0, 5.0
      6, 2, 5, 2.83, 315.00, 30.0, 5.0
[INF] SPCs (nid, first DOF, last DOF, value)
      1, 1, 2, 0.0
      3, 1, 2, 0.0
      5, 1, 1, 0.2
      5, 2, 2, -0.4
[INF] Loads (nid, DOF, value)
      4, 2, -30.0
[INF] Map nodal local DOF to global DOF
      Nodal DOF mapping: {1: [1, 2], 2: [3, 4], 3: [5, 6], 4: [7, 8], 5: [9, 10]}
      Element DOF mapping: {1: [1, 2, 3, 4], 2: [5, 6, 7, 8], 3: [1, 2, 7, 8], 4: [3, 4, 7, 8], 5: [7, 8, 9, 10], 6: [3, 4, 9, 10]}
[INF] Assemble global stiffness matrix K
      Global stiffness matrix size: (10, 10)
[INF] Solve linear system of equations K * u = f
[INF] Postprocess results
      Results written to ./results/res_truss_structure_w_six_members.txt

Done.
```

In the results file the solved nodal displacements and reaction forces and elemental strains and stresss results will be printed out.

<br />

Example output of the result file:

```
 
GUID : 6e174ef1-4746-495d-9429-b530eedde05d
 
                                      R e s u l t s                                       
 
              Linear Static Analysis for Planar Truss and Beam Structures
 
 
 Input file: truss_structure_w_six_members.txt
 
 
                        S u m m a r y   F E   M o d e l   S i z e
 
                              Number of Nodes       5
                              Number of Elements    6
                                     incl.  Rods    6
                                     incl. Beams    0
                              Total Number of DOFs  10
                              Number of SPCs        4
                              Number of Loads       1
                              Location COG          (1.6, -1.2)
 

 
                          N o d a l    D i s p l a c e m e n t s

nid     u_x             u_y             u_R            
1       0.0000e+00      0.0000e+00
2      -7.7404e-02     -8.9633e-01
3       0.0000e+00      0.0000e+00
4      -6.1298e-02     -9.7374e-01
5       2.0000e-01     -4.0000e-01

 
                        N o d a l   R e a c t i o n   F o r c e s

nid     RF_x            RF_y            RM             
1      -1.8389e+01      2.4195e+01
2       0.0000e+00      1.2434e-14
3       4.5974e+00      0.0000e+00
4       0.0000e+00     -3.0000e+01
5       1.3792e+01      5.8053e+00

 
                            E l e m e n t a l   S t r a i n s

eid     Pt. 1           Pt. 2          
1      -3.8702e-02     -3.8702e-02
2      -3.0649e-02     -3.0649e-02
3       2.2811e-01      2.2811e-01
4       3.8702e-02      3.8702e-02
5       1.3065e-01      1.3065e-01
6      -5.4733e-02     -5.4733e-02

 
                           E l e m e n t a l   S t r e s s e s

eid     Pt. 1           Pt. 2          
1      -1.1611e+00     -1.1611e+00
2      -9.1947e-01     -9.1947e-01
3       6.8433e+00      6.8433e+00
4       1.1611e+00      1.1611e+00
5       3.9195e+00      3.9195e+00
6      -1.6420e+00     -1.6420e+00
 

                               E n d   o f   R e s u l t s 
```

<br />

# Summary and Outlook

The FE results provided by the python code have been tested and verified for now for planar truss strucutres (see files in `./test` directory). The FE models for the truss and beam strucutres and their results for verification are taken from [[2]](#2).

Beam and Truss/Beam structures will be tested.

At the moment, only the input file model is visualized by using the `view_input_file.py` script. For future crosschecks of the strucutre behaviour, the deformed shape of the structure will also be available to plot. 

<br />

# References
[[1]](#1) Finite Element Method Lecture Notes, Frankfurt University of Applied Sciences

[[2]](#2) Öchsner, Andreas & Öchsner, Marco. (2018). A First Introduction to the Finite Element Analysis Program MSC Marc/Mentat. 10.1007/978-3-319-71915-3. 
