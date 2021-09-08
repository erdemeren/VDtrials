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

  assert(argc==8);
  const char* modelFile = argv[1];
  //const char* meshFile = argv[2];

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();

  int cell_0 = atoi(argv[2]);

  // NEW PART
  // vd_sim sim_trial;
  //sim_trial.set_mesh(modelFile, meshFile);
  //apf::Mesh2* m = sim_trial.get_mesh();
  vd_sim sim_trial;
  sim_trial.set_adapt_param(atof(argv[3]));
  sim_trial.set_adapt_type(ADAPT_TYPE::ADAPT_STEP);

  sim_trial.set_mesh(modelFile);
  printf("Set mesh, reloading mesh...\n");
  //sim_trial.set_mesh_smb(modelFile, meshFile);
  sim_trial.set_end(0.0);
  sim_trial.reload_mesh();

  sim_trial.set_field_calc(VEL_TYPE::MASON);
  sim_trial.set_vec_sp_calc(PROJ_TYPE::EXT_SHELL);

  //sim_trial.start_sim();
  sim_trial.adapt();
  sim_trial.save_vtk_name("./output/before_all");

  apf::Mesh2* m = sim_trial.get_mesh();

  double x = atof(argv[4]);
  double y = atof(argv[5]);
  double z = atof(argv[6]);

  double shr1 = atof(argv[7]);

  m = sim_trial.get_mesh();
  apf::Vector3 ax;
  double v[3] = {x,y,z};
  ax.fromArray(v);

  vd_entlist* e_list = sim_trial.get_elist();
  cell_base* c_base = sim_trial.get_c_base();
  field_calc* f_calc = sim_trial.get_f_calc();


  sim_trial.save_vtk_name("./output/tetra_base");
  e_list->refresh();
  vd_edisc* e_d = new vd_edisc(m, c_base, e_list);
  e_d->set_proj(PROJ_TYPE::EXT_SHELL);
  e_d->set_vdpar(f_calc->vdparam);
  e_d->set_field_calc(*f_calc);

  e_d->set_len_sh(0.01);
  e_d->set_rho_rat(4, 1);
  e_d->set_0cell(cell_0, false);
  e_d->wg_tag = WG_TYPE::TRI;
  e_d->calc_max_diss_wg();
  e_d->write_diss_csv("./output/w_t_1.csv");
  e_d->wg_tag = WG_TYPE::EDGE;
  e_d->calc_max_diss_wg();
  e_d->write_diss_csv("./output/w_e_1.csv");
  e_d->calc_max_diss_trial_wg();
  e_d->write_diss_csv("./output/w_s_1.csv");
  e_d->try_1cell();
  e_d->try_2cell();
  e_d->write_diss_csv("./output/w_v_1.csv");
  e_d->write_diss_exp_csv("./output/we_v_1.csv");
  e_d->write_ei_csv("./output/e_v_1.csv");

  delete e_d;


  std::vector<double> shr_rat({0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 2});
  double shr_last = 1.;
  for(int i = 0; i < shr_rat.size(); i++) {
    double shr_curr = 1./shr_last*shr_rat.at(i)*shr1;
    shr_last = shr_rat.at(i);
    vd_shr_axis(m, ax, shr_curr);
    std::string temp("");
    temp = "./output/tetra_" + std::to_string(i);
    sim_trial.save_vtk_name(temp.c_str());
    e_list->refresh();
    vd_edisc* e_d = new vd_edisc(m, c_base, e_list);
    e_d->set_proj(PROJ_TYPE::EXT_SHELL);
    e_d->set_vdpar(f_calc->vdparam);
    e_d->set_field_calc(*f_calc);

    e_d->set_len_sh(0.01);
    e_d->set_rho_rat(4, 1);
    e_d->set_0cell(cell_0, false);
    e_d->wg_tag = WG_TYPE::TRI;
    e_d->calc_max_diss_wg();
    temp = "./output/w_t_" + std::to_string(i) + ".csv";
    e_d->write_diss_csv(temp.c_str());
    e_d->wg_tag = WG_TYPE::EDGE;
    e_d->calc_max_diss_wg();
    temp = "./output/w_e_" + std::to_string(i) + ".csv";
    e_d->write_diss_csv(temp.c_str());
    e_d->calc_max_diss_trial_wg();
    temp = "./output/w_s_" + std::to_string(i) + ".csv";
    e_d->write_diss_csv(temp.c_str());
    e_d->try_1cell();
    e_d->try_2cell();
    temp = "./output/w_v_" + std::to_string(i) + ".csv";
    e_d->write_diss_csv(temp.c_str());
    temp = "./output/we_v_" + std::to_string(i) + ".csv";
    e_d->write_diss_exp_csv(temp.c_str());
    temp = "./output/e_v_" + std::to_string(i) + ".csv";
    e_d->write_ei_csv(temp.c_str());

    delete e_d;
  }

  printf("started sim, cleaning up...\n");
  sim_trial.clean_up();
  PCU_Comm_Free();
  MPI_Finalize();
}

