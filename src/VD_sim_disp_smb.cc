// Adapted from ma_test.cc in SCOREC test directory.

#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>

#include <PCU.h>

#include <apfNumbering.h>
#include <apfShape.h>

// VD related modules:
/* Simple operators for random number calculations, etc. Currently generates 
   random normalized vector:  */
#include "topo_rand.h"
/* Topology manipulator:  */
#include "topo_manip.h"
#include "topo_disp.h"
/* Topology extractor:  */
#include "topo_extinfo.h"
/* Topology write:  */
#include "topo_write.h"
/* Topology geometry:  */
#include "topo_geom.h"
/* Topology field:  */
#include "topo_field.h"
/* Topology lens:  */
#include "topo_lens.h"
#include "topo_glens.h"

#include "topo_vd.h"

#include "topo_topo.h"
#include "topo_graph.h"

int main(int argc, char** argv)
{

  assert(argc==10);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();

  const char* extFile = argv[3];
  const char* eqnFile = argv[4];
  const char* timeFile = argv[5];

  bool ins_flag = (atoi(argv[6]) != 0);
  bool col_flag = (atoi(argv[7]) != 0);

  std::string extStr(extFile);
  std::string eqnStr(eqnFile);
  std::string timeStr(timeFile);

  double ratio_ins_in = atof(argv[8]);
  double ratio_col_in = 2.6*ratio_ins_in;
  double ratio_col_sh_in = 1.3*ratio_ins_in;

  double ad_param = atof(argv[9]);

  // NEW PART
  vd_sim sim_trial;

  sim_trial.set_ins_flag(ins_flag);
  sim_trial.set_col_flag(col_flag);

  sim_trial.set_adapt(1);
  sim_trial.set_adapt_param(ad_param);
  //sim_trial.set_adapt_type(ADAPT_TYPE::ADAPT_STEP);
  sim_trial.set_adapt_type(ADAPT_TYPE::ADAPT_STEP_1CELL);
  //sim_trial.set_adapt_type(ADAPT_TYPE::ADAPT_BOUND);

  sim_trial.set_mesh_smb(modelFile, meshFile);
  //sim_trial.reload_mesh();

  sim_trial.set_field_calc(VEL_TYPE::MASON);
  sim_trial.set_vec_sp_calc(PROJ_TYPE::FIXED);
  sim_trial.set_integ_type(INTEG_TYPE::RK2);

  sim_trial.set_ins_rat(ratio_ins_in);
  sim_trial.set_col_rat(ratio_col_sh_in, ratio_col_in);

  printf("Set mesh, reloading mesh...\n");
  sim_trial.set_mov_flag(false);

  std::string ext_opts(extFile);
  sim_trial.set_extract(true, ext_opts);

  std::vector<std::vector<std::string> > temp_str_vec(0, std::vector<std::string> (0,""));
  ReadNames(eqnStr, ";", temp_str_vec);
  sim_trial.set_equations(temp_str_vec.at(0));

  ReadNames(timeStr, ";", temp_str_vec);
  sim_trial.set_time(temp_str_vec.at(0));

  sim_trial.start_sim();

  printf("started sim, cleaning up...\n");
  sim_trial.clean_up();

  PCU_Comm_Free();
  MPI_Finalize();
}

