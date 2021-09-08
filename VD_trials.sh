F_DEBUG="F_DEBUG=1"
if [ $# == 2 ]; then
	F_DEBUG="F_DEBUG="$2;
fi

########################################
# SCOREC Interface tests:
########################################
# From SCOREC tests. Generate a 2-dim mesh, and add the 3rd 
# dimension. Add a number of triangles and a tetrahedron bounded by these. 

if [ $1 == 10 ]; then
	echo "SCOREC 1: Add entities"
	rm ./newdim
	make newdim $F_DEBUG
	./newdim
 exit;
fi 

# Using the SCOREC interface, load a gmsh mesh into the program. Verify and exit.

if [ $1 == 11 ]; then
	echo "SCOREC 2: Load mesh"
	make sim
  ./vdsim mshfiles/boxnbox.tess mshfiles/boxnbox.msh
 exit;
fi

# Using the SCOREC interface, load a gmsh mesh into the program. Refine the mesh
# Using meshadapt. Verify and exit.
if [ $1 == 12 ]; then
	echo "SCOREC 2: Refine the mesh 2x"
	make refine2x
	make split
  #mpirun -n 1 ./refine2x mshfiles/beforecol.dmg mshfiles/beforecol.smb 2 output/afterrefine
  mpirun -n 1 ./refine2x mshfiles/boxnbox.tess 2 output/afterrefine
  mpirun -n 2 ./sp_tess mshfiles/boxnbox.tess output/afterrefine 2
 exit;
fi

########################################
# VDlib
########################################

########################################
# To generate the article data:
########################################

##################
# Method
##################

if [ $1 == 21 ]; then
	echo "PRM 2021 seven grain simulation (spurious insertions)"
	make disptess $F_DEBUG
	mkdir -p ./sims
	mkdir -p ./sims/method1
	mkdir -p ./sims/method1/mshfiles
	cp vddisp ./sims/method1/vddisp
	cp Makefile_sim ./sims/method1/Makefile
	cp ./mshfiles/method/boxnbox.tess ./sims/method1/mshfiles/boxnbox.tess
	cp ./mshfiles/eqn_mason_fixed_rk2.txt ./sims/method1/mshfiles/eqn_mason_fixed_rk2.txt
	cp ./mshfiles/method/time_opts_boxnbox.txt ./sims/method1/mshfiles/time_opts_boxnbox.txt
	cp ./mshfiles/ext_opts_none.txt ./sims/method1/mshfiles/ext_opts_none.txt
	cp ./mshfiles/benchmark/adapt_opts_1cell_3.txt ./sims/method1/mshfiles/adapt_opts_1cell_3.txt
	cd ./sims/method1
	make folder
  #./vddisp tess_file / extract_options / equations_of_motion / time_options / mesh_adaptation_options / insertion_flag / collapse_flag / insert2refinement_length_ratio / vtk_sub_save_flag 
	./vddisp ./mshfiles/boxnbox.tess ./mshfiles/ext_opts_none.txt ./mshfiles/eqn_mason_fixed_rk2.txt ./mshfiles/time_opts_boxnbox.txt ./mshfiles/adapt_opts_1cell.txt 1 1 5. 0 > verb.txt
#./vddisp ./mshfiles/boxnbox.tess ./mshfiles/ext_opts_none.txt ./mshfiles/eqn_mason_fixed_rk2.txt ./mshfiles/time_opts_boxnbox.txt ./mshfiles/adapt_opts_boxnbox_4.txt 1 1 5. 0 > verb.txt
#./vddisp ./mshfiles/boxnbox.tess ./mshfiles/ext_opts_none.txt ./mshfiles/eqn_mason_shell_rk2.txt ./mshfiles/time_opts_boxnbox.txt ./mshfiles/adapt_opts_boxnbox_4.txt 1 1 5. 0 > verb.txt
	cd ../..
 exit;
fi 

if [ $1 == 22 ]; then
	echo "PRM 2021 insertion results generation (Fig. 16). It also generates the results for two different insertion energy calculations."
	make edisc4tessdiss $F_DEBUG
	mkdir -p ./sims
	mkdir -p ./sims/method2
	mkdir -p ./sims/method2/mshfiles
	cp vdedisc ./sims/method2/vdedisc
	cp Makefile_sim ./sims/method2/Makefile
	cp ./m/diss_comp.m ./sims/method2/diss_comp.m
	cp ./mshfiles/method/Tetrahedral2c.tess ./sims/method2/mshfiles/Tetrahedral2c.tess
	cd ./sims/method2
	make folder
	./vdedisc mshfiles/Tetrahedral2c.tess 7 "2.1" 0 0 1 "1.4" > verb.txt
	cd ../..
 exit;
fi 

if [ $1 == 23 ]; then
	echo "PRM 2021 100 grain simulation (Fig. 19)"
	make disptess $F_DEBUG
	mkdir -p ./sims
	mkdir -p ./sims/method3
	mkdir -p ./sims/method3/mshfiles
	cp vddisp ./sims/method3/vddisp
	cp Makefile_sim ./sims/method3/Makefile
	cp ./mshfiles/method/gigi100.tess ./sims/method3/mshfiles/gigi100.tess
	cp ./mshfiles/method/gigi100_voronoi.tess ./sims/method3/mshfiles/gigi100_voronoi.tess
	cp ./mshfiles/ext_opts.txt ./sims/method3/mshfiles/ext_opts.txt
	cp ./mshfiles/eqn_mason_shell_rk2.txt ./sims/method3/mshfiles/eqn_mason_shell_rk2.txt
	cp ./mshfiles/method/time_opts_gigi100.txt ./sims/method3/mshfiles/time_opts.txt
	cp ./mshfiles/ext_opts.txt ./sims/method3/mshfiles/ext_opts.txt
	cp ./mshfiles/adapt_opts_1cell.txt ./sims/method3/mshfiles/adapt_opts_1cell.txt
	cd ./sims/method3
	make folder
	# For regularized grains:
  #./vddisp tess_file / extract_options / equations_of_motion / time_options / mesh_adaptation_options / insertion_flag / collapse_flag / insert2refinement_length_ratio / vtk_sub_save_flag 
	./vddisp ./mshfiles/gigi100.tess ./mshfiles/ext_opts.txt ./mshfiles/eqn_mason_shell_rk2.txt ./mshfiles/time_opts.txt ./mshfiles/adapt_opts_1cell.txt 1 1 5. 0 > verb.txt
	cd ../..
	# For the initial condition used in the paper:
	#./vddisp ./mshfiles/gigi100_voronoi.tess ./mshfiles/ext_opts.txt ./mshfiles/eqn_mason_shell_rk2.txt ./mshfiles/time_opts.txt ./mshfiles/adapt_opts_1cell.txt 1 1 5. > verb.txt
 exit;
fi 
#./vddisp ./mshfiles/sph_fibo_sph0.dmg ./mshfiles/sph_fibo_sph0.smb ./mshfiles/ext_opts_sph_count.txt ./mshfiles/eqn_mason_NBC_shell_rk2.txt ./mshfiles/time_opts_0.txt ./mshfiles/adapt_opts_0.txt 0 0 50.
#	make disptess
#./vddisp ./mshfiles/gigi100.tess ./mshfiles/ext_opts.txt ./mshfiles/eqn_mason_shell_rk2.txt ./mshfiles/time_opts.txt ./mshfiles/adapt_opts_1cell.txt 1 1 5.

##################
#	CoM
##################

if [ $1 == 24 ]; then
	echo "PRB 2021 Relax the Kelvin cell configuration."
	make benchrcbound $F_DEBUG
	mkdir -p ./sims
	mkdir -p ./sims/CoM
	mkdir -p ./sims/CoM/mshfiles
	cp vdbenche ./sims/CoM/vdbenche
	cp Makefile_sim ./sims/CoM/Makefile
	cp ./mshfiles/CoM/troctahedra_unit.tess ./sims/CoM/mshfiles/troctahedra_unit.tess
	cd ./sims/CoM
	make folder
  #./vdbenche f_name / max_ref / mult / dt / g_tol
	./vdbenche ./tempmesh/gauss 0.6 2. 500. 0.001
	cd ../..
 exit;
fi

##################
# Benchmarks
##################
#	Benchmarks
# Optional
# Generates the initial conditions for the benchmark cases in https://arxiv.org/abs/2203.03167.
# Or download the mesh files from https://github.com/erdemeren/VDtrials_meshes.

# To simulate the three cases for different levels of refinement:
# 4XY: X for the case (0: shrinking sphere, 1: triple line, 2: quadruple point)
#      Y for the refinement: 0-5 for increasing levels of refinement. Please 
# check ./mshfiles/benchmark/adapt_opts_X_Y.txt for the refinement edge length. 
if [[ $1 =~ 40([1-5]+) ]]; then
  ref_lev=${BASH_REMATCH[1]}
  sim_dir="./sims/benchmark/trial_sph_"$ref_lev;
	echo "CMS simulate the shrinking sphere case for refinement level "$ref_lev
	make dispsmb $F_DEBUG
	mkdir -p ./sims
	mkdir -p ./sims/benchmark
	mkdir -p $sim_dir
	mkdir -p $sim_dir/mshfiles
	cp vddisp $sim_dir/vddisp
	cp Makefile_sim $sim_dir/Makefile
	cp ./mshfiles/benchmark/sph_fibo_sph0.dmg $sim_dir/mshfiles/sph_fibo_sph0.dmg
	cp "./mshfiles/benchmark/sph_fibo_sph"$ref_lev"0.smb" $sim_dir"/mshfiles/sph_fibo_sph"$ref_lev"0.smb"
	cp ./mshfiles/eqn_mason_NBC_shell_rk2.txt $sim_dir/mshfiles/eqn_mason_NBC_shell_rk2.txt
	cp ./mshfiles/benchmark/time_opts_sph.txt $sim_dir/mshfiles/time_opts_sph.txt
	cp ./mshfiles/benchmark/ext_opts_sph.txt $sim_dir/mshfiles/ext_opts_sph.txt
	cp ./mshfiles/benchmark/move_list_sph.txt $sim_dir/mshfiles/move_list_sph.txt
	cp ./mshfiles/benchmark/ref_list_sph.txt $sim_dir/mshfiles/ref_list_sph.txt
	cp ./mshfiles/benchmark/adapt_opts_sph_$ref_lev.txt $sim_dir/mshfiles/adapt_opts_sph_$ref_lev.txt
	cd $sim_dir
	make folder

  #./vddisp tess_file / extract_options / equations_of_motion / time_options / mesh_adaptation_options / insertion_flag / collapse_flag / insert2refinement_length_ratio / vtk_sub_save_flag 
	./vddisp ./mshfiles/sph_fibo_sph0.dmg ./mshfiles/sph_fibo_sph$ref_lev.smb ./mshfiles/ext_opts_sph.txt ./mshfiles/eqn_mason_NBC_shell_rk2.txt ./mshfiles/time_opts_sph.txt ./mshfiles/adapt_opts_sph_$ref_lev.txt 0 0 25. 0 > verb.txt
 exit;
fi

if [[ $1 =~ 41([1-5]+) ]]; then
  ref_lev=${BASH_REMATCH[1]}
  sim_dir="./sims/benchmark/trial_trip_"$ref_lev;
	echo "CMS simulate the triple line case for refinement level "$ref_lev
	make dispsmb $F_DEBUG
	mkdir -p ./sims
	mkdir -p ./sims/benchmark
	mkdir -p $sim_dir
	mkdir -p $sim_dir/mshfiles
	cp vddisp $sim_dir/vddisp
	cp Makefile_sim $sim_dir/Makefile
	cp ./mshfiles/benchmark/trip_fibo_tri0.dmg $sim_dir/mshfiles/trip_fibo_tri0.dmg
	cp "./mshfiles/benchmark/trip_fibo_tri"$ref_lev"0.smb" $sim_dir"/mshfiles/trip_fibo_tri"$ref_lev"0.smb"
	cp ./mshfiles/eqn_mason_NBC_shell_rk2.txt $sim_dir/mshfiles/eqn_mason_NBC_shell_rk2.txt
	cp ./mshfiles/benchmark/time_opts_trip.txt $sim_dir/mshfiles/time_opts_trip.txt
	cp ./mshfiles/benchmark/ext_opts_trip.txt $sim_dir/mshfiles/ext_opts_trip.txt
	cp ./mshfiles/benchmark/move_list_trip.txt $sim_dir/mshfiles/move_list_trip.txt
	cp ./mshfiles/benchmark/ref_list_trip.txt $sim_dir/mshfiles/ref_list_trip.txt
	cp ./mshfiles/benchmark/adapt_opts_trip_$ref_lev.txt $sim_dir/mshfiles/adapt_opts_trip_$ref_lev.txt
	cd $sim_dir
	make folder

  #./vddisp tess_file / extract_options / equations_of_motion / time_options / mesh_adaptation_options / insertion_flag / collapse_flag / insert2refinement_length_ratio / vtk_sub_save_flag 
	./vddisp ./mshfiles/trip_fibo_tri0.dmg ./mshfiles/trip_fibo_tri$ref_lev.smb ./mshfiles/ext_opts_trip.txt ./mshfiles/eqn_mason_NBC_shell_rk2.txt ./mshfiles/time_opts_trip.txt ./mshfiles/adapt_opts_trip_$ref_lev.txt 0 0 25. 0 > verb.txt
 exit;
fi

if [[ $1 =~ 42([1-5]+) ]]; then
  ref_lev=${BASH_REMATCH[1]}
  sim_dir="./sims/benchmark/trial_hex_"$ref_lev;
	echo "CMS simulate the hexagonal case for refinement level "$ref_lev
	make dispsmb $F_DEBUG
	mkdir -p ./sims
	mkdir -p ./sims/benchmark
	mkdir -p $sim_dir
	mkdir -p $sim_dir/mshfiles
	cp vddisp $sim_dir/vddisp
	cp Makefile_sim $sim_dir/Makefile
	cp ./mshfiles/benchmark/hex_fibo_tri0.dmg $sim_dir/mshfiles/hex_fibo_tri0.dmg
	cp "./mshfiles/benchmark/hex_fibo_tri"$ref_lev"0.smb" $sim_dir"/mshfiles/hex_fibo_tri"$ref_lev"0.smb"
	cp ./mshfiles/eqn_mason_NBC_shell_rk2.txt $sim_dir/mshfiles/eqn_mason_NBC_shell_rk2.txt
	cp ./mshfiles/benchmark/time_opts_hexa.txt $sim_dir/mshfiles/time_opts_hexa.txt
	cp ./mshfiles/benchmark/ext_opts_hexa.txt $sim_dir/mshfiles/ext_opts_hexa.txt
	cp ./mshfiles/benchmark/move_list_hex.txt $sim_dir/mshfiles/move_list_hex.txt
	cp ./mshfiles/benchmark/ref_list_hex.txt $sim_dir/mshfiles/ref_list_hex.txt
	cp ./mshfiles/benchmark/adapt_opts_hex_$ref_lev.txt $sim_dir/mshfiles/adapt_opts_hex_$ref_lev.txt
	cd $sim_dir
	make folder

  #./vddisp tess_file / extract_options / equations_of_motion / time_options / mesh_adaptation_options / insertion_flag / collapse_flag / insert2refinement_length_ratio / vtk_sub_save_flag 
	./vddisp ./mshfiles/hex_fibo_tri0.dmg ./mshfiles/hex_fibo_tri$ref_lev.smb ./mshfiles/ext_opts_hexa.txt ./mshfiles/eqn_mason_NBC_shell_rk2.txt ./mshfiles/time_opts_hexa.txt ./mshfiles/adapt_opts_hex_$ref_lev.txt 0 0 25. 0 > verb.txt
 exit;
fi

#if [ $1 == 31 ]; then
#fi 

########################################
# Individual tests
########################################

# Using the VDlib simulation interface, load a gmsh mesh into the program. 
# Deform it under grain growth until a cell collapse is necessary. 
# Collapse the cell.
# Verify and exit.

