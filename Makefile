F_DEBUG ?= 1
ifeq ($(F_DEBUG), 1)
  CXXFLAGS= -gdwarf-2 -g3 -g -pg -std=gnu++11 -pedantic-errors -Wuninitialized -no-pie
else
  CXXFLAGS=-O2 -lmpi -std=gnu++11 -Wreturn-type -no-pie
endif

CC=mpicxx
#CXXFLAGS= -O2 -pg -std=gnu++11 -pedantic-errors -Wuninitialized -no-pie
SCORECDIR = <SCOREC DIRECTORY>
ZOL_PARMETISDIR = <ZOLTAN/PARMETIS DIRECTORY>
GSL_DIR = <GSL DIRECTORY>

VDLIBS = -I$(VDDIR)/src  -L$(VDDIR)/lib -I$(GSL_DIR)/include -L$(GSL_DIR)/lib -lvdtopo -lgsl -lgslcblas
#SCORECLIBS = -I$(SCORECDIR)/include -I$(ZOL_PARMETISDIR)/include -L$(SCORECDIR)/lib -L$(ZOL_PARMETISDIR)/lib -lma -lparma -lapf_zoltan -lzoltan -lparmetis -lmetis -lmds -lapf -lgmi -lpcu -lmpi

SCORECLIBS = -I$(SCORECDIR)/include -I$(ZOL_PARMETISDIR)/include -L$(SCORECDIR)/lib -L$(ZOL_PARMETISDIR)/lib -lma -lparma -lapf_zoltan -lmds -lapf -lgmi -llion -lcrv -lmth -lph -lsam -lspr -lpumi -lpcu -lzoltan -lparmetis -lmetis

binaries=vdsim vdcirc vdedisc vdedisc newdim vdcol vddisp refine2x vdcirc vdsim2 vdcol

# --- targets

folder:
	mkdir -p ./output
	mkdir -p ./output/sub
	mkdir -p ./outputiter
	mkdir -p ./outputtemp
	mkdir -p ./movie
	mkdir -p ./movie/gif
	mkdir -p ./tempmesh

all: 
	$(CC) ./src/VD_sim.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdsim
	$(CC) ./src/VD_circ.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdcirc
	$(CC) ./src/VD_edisc.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
	$(CC) ./src/newdim.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o newdim
	$(CC) ./src/VD_sim_col.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdcol
	$(CC) ./src/VD_sim_disp.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

VTK: 
	$(CC) ./src/VD_vtk_mesh.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdvtk

VTKtess: 
	$(CC) ./src/VD_vtk_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdvtk

adapt:
	$(CC) ./src/VD_sim_adapt.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdadapt

vector:
	$(CC) ./src/apfvector.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdvect

refine2x:
	$(CC) ./src/refine2x.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o refine2x

sim: 
	$(CC) ./src/VD_sim.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdsim
sim2: 
	$(CC) ./src/VD_sim2.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdsim2
custom: 
	$(CC) ./src/VD_sim_custom.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdsim

cbase: 
	$(CC) ./src/VD_sim_cbase.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdcbase

elist: 
	$(CC) ./src/VD_sim_elist.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdsim

#####################################
#          Collapse tests           #
#####################################
colstress: 
	$(CC) ./src/VD_sim_collapse_stress.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdstress

#####################################
#         Insertion tests           #
#####################################
edisc: 
	$(CC) ./src/VD_edisc.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc2: 
	$(CC) ./src/VD_edisc2.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc3: 
	$(CC) ./src/VD_edisc3.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc4: 
	$(CC) ./src/VD_edisc4.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc4tess: 
	$(CC) ./src/VD_edisc4_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc4tessstr: 
	$(CC) ./src/VD_edisc4_tess_str.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc4tessdiss:
	$(CC) ./src/VD_edisc4_tess_diss.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc

gmireloadtest: 
	$(CC) ./src/APF_model_change_test.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc

edisc5: 
	$(CC) ./src/VD_edisc5.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc6: 
	$(CC) ./src/VD_edisc6.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc6tess:
	$(CC) ./src/VD_edisc6_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc7: 
	$(CC) ./src/VD_edisc7.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc7smb: 
	$(CC) ./src/VD_edisc7_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc8smb: 
	$(CC) ./src/VD_edisc8_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc
edisc9smb: 
	$(CC) ./src/VD_edisc9_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdedisc

newdim: 
	$(CC) ./src/newdim.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o newdim

col:
	rm -f $(binaries)
	$(CC) ./src/VD_sim_col.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdcol

colr: 
	rm -f $(binaries)
	$(CC) ./src/VD_sim_col_rand.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdcol

disp: 
	$(CC) ./src/VD_sim_disp.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

# Compare dissipation rate calculations vs. the actual dissipation rate, for a
# quadruple point for changing edge angles, including tetrahedral angle.
diss: 
	$(CC) ./src/VD_sim_bench_diss.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddiss

#####################################
#             Extract MS            #
#####################################
MSsmb: 
	$(CC) ./src/VD_sim_MS_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdMS
MSevolve: 
	$(CC) ./src/VD_sim_MS_evolve.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdMSevol

extcasessmb:
	$(CC) ./src/VD_sim_ext_cases_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdextfiles

#####################################
#        Extract file list          #
#####################################
extfiles: 
	$(CC) ./src/VD_ext_files.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdextfiles

#####################################
#         Evolution programs        #
#####################################
dispsmb: 
	$(CC) ./src/VD_sim_disp_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

# Variable time step. Last 1/10th of the sim. at second time step.
dispsmbvt: 
	$(CC) ./src/VD_sim_disp_smb_vt.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

dispsmbevol:
	$(CC) ./src/VD_sim_disp_smb_evol.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

disptessevol:
	$(CC) ./src/VD_sim_disp_tess_evol.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

disptess: 
	$(CC) ./src/VD_sim_disp_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

dispo: 
	$(CC) ./src/VD_sim_disp_open.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

dispoc: 
	$(CC) ./src/VD_sim_disp_open_col.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

eomtess:
	$(CC) ./src/VD_sim_eom_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

eomtessiter:
	$(CC) ./src/VD_sim_eom_tess_bench.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddispiter

kuprattess:
	$(CC) ./src/VD_sim_disp_tess_kuprat.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vddisp

#####################################
#         Generate test cases       #
#####################################
# TODO Clear out the outdated tests. Awful naming convention due to multiple
# subsequent updates to the way meshes were generated. Clean it up after the 
# benchmark articles.
dispevolboundtess: 
	$(CC) ./src/VD_sim_evol_bound_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdevol

disprelaxtriptess: 
	$(CC) ./src/VD_sim_relax_trip_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtrip

disprelaxquadtess: 
	$(CC) ./src/VD_sim_relax_quad_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdquad

disprefinetripsmb: 
	$(CC) ./src/VD_sim_refine_trip_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtrip

disprefinequadsmb: 
	$(CC) ./src/VD_sim_refine_quad_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdquad

disprelaxtripspurtess: 
	$(CC) ./src/VD_sim_relax_spur_trip_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtrip

disprelaxquadspurtess: 
	$(CC) ./src/VD_sim_relax_spur_quad_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdquad

disprelaxtripspursphtess:
	$(CC) ./src/VD_sim_relax_spur_trip_sph_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtrip

disprelaxquadspursphtess: 
	$(CC) ./src/VD_sim_relax_spur_quad_sph_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdquad

disprelaxquadspur1sphtess: 
	$(CC) ./src/VD_sim_relax_spur_quad_sing_sph_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdquad

disprelaxtripspur1sphtess: 
	$(CC) ./src/VD_sim_relax_spur_trip_sing_sph_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtrip

disprelaxtripspur1sphntess: 
	$(CC) ./src/VD_sim_relax_spur_trip_sing_sph_new_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtrip

# 1/8th sphere
disprelaxsphtess: 
	$(CC) ./src/VD_sim_relax_spur_prism_sing_sph_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdsph

disprelaxhexsphtess: 
	$(CC) ./src/VD_sim_relax_spur_hex_sing_sph_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdhex

#####################################
#         Evolve test cases         #
#####################################
dispbenchsmb: 
	$(CC) ./src/VD_sim_benchmark_cases_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbench

# Using the vd_sim interfaces.
dispbenchsimsmb: 
	$(CC) ./src/VD_sim_benchmark_cases_sim_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbench

# Using the vd_sim interfaces. Assumes the first mesh is a spherical shell with voids.
dispbenchsimsmbtri: 
	$(CC) ./src/VD_sim_benchmark_cases_sim_smb_tri.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbench

# Using the vd_sim interfaces. Assumes the first mesh is a spherical shell with voids. Assumes meshes are sphere, triple line, hexagonal, triple line flat, hexagonal flat
dispbenchsimsmbtrihex: 
	$(CC) ./src/VD_sim_benchmark_cases_sim_smb_tri_hex.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbench

# Only process the hausdorff distance calculation from existing meshes.
dispbenchsimsmbtrihexhaus: 
	$(CC) ./src/VD_sim_benchmark_cases_sim_smb_tri_hex_hausonly.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbench2

# Using the vd_sim interfaces. Assumes the first mesh is a spherical shell with voids. Assumes meshes are sphere, new triple line, hexagonal, triple line flat, hexagonal flat
dispbenchsimsmbsthsph:
	$(CC) ./src/VD_sim_benchmark_cases_smb_sth_sph.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbench3

# Using the vd_sim interfaces. Assumes the first mesh is a spherical shell with voids. Assumes meshes are sphere, new triple line, hexagonal, 1/8th sphere. 
# Doesn't evolve the sphere.
dispbenchsimsmbsthsph2:
	$(CC) ./src/VD_sim_benchmark_cases_smb_sth_sph2.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbench3

# Extract the off file from a single mesh.
dispextoff:
	$(CC) ./src/VD_sim_ext_off.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o extoff

# Benchmark evolution of a single MS.
dispbenchsmb_sing: 
	$(CC) ./src/VD_sim_bench_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbench

#####################################
#        Constant of motion         #
#####################################
# Given a mesh and extraction options for constant of motion, extract the relevant quantities.
# COM: Complete microstructure
# COM_GRAIN: Per grain
benchsingle: 
	$(CC) ./src/VD_sim_bench_smb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbench

benchmulti: 
	$(CC) ./src/VD_sim_bench_constant.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenchm

# Conjugate gradient based surface minimization
benchconj: 
	$(CC) ./src/VD_sim_bench_GAUSS.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenchg

# Given a maximum adaptation parameter, evolve the surface until the maximum of 
# the norms of the velocities is below a threshold. Lower the time step as the
# maximum of the norm decreases. First, start with adaptation parameter at 2.
# increase by set step until it exceeds the maximum adaptation parameter.
# 
benchevol: 
	$(CC) ./src/VD_sim_bench_GAUSS_evol.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

# Mirror BC. Specifically used to find the equilibrium metastable shape of 
# the truncated octahedron by reflecting the energetic triangles by the exterior 
# surfaces of the boundary vertices.
benchmirror:
	$(CC) ./src/VD_sim_bench_mirror_shell.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

# Mirror BC. Specifically used to find the equilibrium metastable shape of 
# the truncated octahedron by reflecting the energetic triangles by the exterior 
# surfaces of the boundary vertices. Also, the configuration of S0 can be rotated 
# to be matched with the first S0. By using these rotations, motion of the 
# interior S0 are matched and the average motion is calculated to be applied to
# all interior S0.
benchmirrors0:
	$(CC) ./src/VD_sim_bench_mirror_shell_s0.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

# Fixed interior s0.
benchmirrors0fix:
	$(CC) ./src/VD_sim_bench_mirror_shell_s0_fix.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

# Fixed interior s0 without mesh adaptation at each substep. Relies on massaging 
# the interior vertex locations.
benchmirrors0fixrc:
	$(CC) ./src/VD_sim_bench_mirror_shell_s0_fix_rc.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

# Fixed interior s0 without mesh adaptation at each substep. Relies on massaging 
# the interior vertex locations. The surface minimization is conjugate gradient 
# based, as well.
benchrc:
	$(CC) ./src/VD_sim_bench_conj_grad.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

benchrcbound:
	$(CC) ./src/VD_sim_bench_conj_grad_bound.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

benchrcboundlens:
	$(CC) ./src/VD_sim_bench_conj_grad_bound_lens.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

# Fixed interior s1 without mesh adaptation at each substep. Relies on massaging 
# the interior vertex locations. The surface minimization is conjugate gradient 
# based, as well.
benchrctrip:
	$(CC) ./src/VD_sim_bench_conj_grad_trip.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

# Mirror rotation. Specifically used to find the equilibrium metastable shape of 
# the truncated octahedron by only moving the local environment of a single 0-
# stratum
benchmirrot:
	$(CC) ./src/VD_sim_bench_mir_rot_local.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenche

#####################################
#          GRAPH TESTERS            #
#####################################
circ: 
	$(CC) ./src/VD_circ.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdcirc

graph: 
	$(CC) ./src/VD_graph.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdgraph

graphtess: 
	$(CC) ./src/VD_graph_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdgraph

tmeshtess: 
	$(CC) ./src/VD_tmesh_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtmesh

tmeshcb: 
	$(CC) ./src/VD_tmesh_cb.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtmesh

emeshtess:
	$(CC) ./src/VD_sim_ext_mesh.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdemesh

#####################################
#    Generate sphere in a volume    #
#####################################
sphere:
	$(CC) ./src/VD_sim_load_sph.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o sph

# A sphere just consisting of boundary 2cell triangles and tets bounded by them.
spheretri:
	$(CC) ./src/VD_sim_load_sph_tri.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o sph

sphereref:
	$(CC) ./src/VD_sim_refine_sph.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o sph

sphevolve: 
	$(CC) ./src/VD_sim_evolve_sphere.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o sphevol

#####################################
#          Test programs            #
#####################################
topo:
	$(CC) ./src/VD_topo.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtopo  

tess:
	$(CC) ./src/VD_sim_load_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtess

ltess:
	$(CC) ./src/VD_load_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdtess

matest: 
	$(CC) ./src/ma_test.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o matest  

malarge: 
	$(CC) ./src/ma_test_large_angle.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o malarge

macrash: 
	$(CC) ./src/ma_crash.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o macrash

memleak: 
	rm -rf memleak
	$(CC) ./src/memleak.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o memleak

cinsinit: 
	rm -rf memleak
	$(CC) ./src/c_ins_init.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o cinsinit

repart: 
	$(CC) ./src/repart_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o repart
reposition: 
	$(CC) ./src/reposition.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o reposition

disprepart: 
	$(CC) ./src/VD_sim_disp_tess_repart.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o repart

split: 
	$(CC) ./src/split_tess.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o sp_tess

split_ref: 
	$(CC) ./src/split_refine.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o sp_ref

split_eqn: 
	$(CC) ./src/VD_eqn_mason_dist.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o sp_tess

split_eqn_kuprat: 
	$(CC) ./src/VD_eqn_kuprat_gsl.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o sp_tess

refine_field: 
	$(CC) ./src/refine_fieldcollect.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o sp_tess

apfvec:
	$(CC) ./src/VD_apfVector.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdapfvec

apftrifield:
	$(CC) ./src/apf_tri_field.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o trifield

apfjacobian:
	$(CC) ./src/apf_integrate_test.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o jacob

cutplane: 
	$(CC) ./src/VD_cut_mesh_plane.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdcut

cutsph: 
	$(CC) ./src/VD_cut_mesh_sph.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdcut

benchCS: 
	$(CC) ./src/VD_bench_CS.cc $(CXXFLAGS) $(VDLIBS) $(SCORECLIBS) -o vdbenchCS

clean:
	rm -rf ./*~ ./src/*~ ./output/* ./outputiter/* ./outputtemp/* ./tempmesh/*
	mkdir ./output/sub
	echo "" > ./verb.txt
cleanbuild:
	rm -f $(binaries)
cleanmovie:
	rm -f ./movie/*.png

