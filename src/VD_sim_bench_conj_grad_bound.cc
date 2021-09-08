// Adapted from ma_test.cc in SCOREC test directory.

#include <cmath>

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
#include "topo_vd.h"

#include "topo_topo.h"
#include "topo_graph.h"

#include "topo_cg.h"

#include "topo_ma.h"

#define PI_L std::acos(-1.)  /* pi */

int lookup_tri_e_local [3][2] = {{0,2},{0,1},{1,2}};
int lookup_tri_vd_local [3][2] = {{1,2},{0,2},{1,0}};
int lookup_v_tri_local [4] = {2,3,1,0};
apf::Matrix3x3 const M_EYE(1,0,0,0,1,0,0,0,1);

// Given an string vector, some entities with the name "[file].csv" replace the 
// entity with "[file]_suffix.csv".
void add_suffix_csv_sub(std::vector<std::string> &opts, std::string &suffix) {
  std::smatch sm;
  std::string r_exp(".csv");
  for(int i = 0; i < opts.size(); i++) {
    if(match_regex(opts.at(i), sm, r_exp)) {
      opts.at(i) = insert_between(opts.at(i), r_exp, suffix);
    }
  }
}

void add_suffix_csv(std::vector<std::vector<std::string> > &opts, 
                                                  std::string suffix) {
  for(int i = 0; i < opts.size(); i++) {
    add_suffix_csv_sub(opts.at(i), suffix);
  }
}

// Given a csv file, append the opts to the csv.
void append_opts_csv(std::string &filename, std::vector<std::string> &opts) {
  csvfile csv_curr(filename, std::string(";"), std::string(""));

  for(int i = 0; i < opts.size(); i++)
    csv_curr << opts.at(i).c_str();
  csv_curr << endrow;
}

// Create a csv file from a vector of vector of strings. 
void create_opts_csv(std::string &filename, 
                              std::vector<std::vector<std::string> > &opts) {
  std::ofstream ofs;
  ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
  //csvfile csv_curr(filename);
  //csv_curr.flush();
  
  for(int i = 0; i < opts.size(); i++)
    append_opts_csv(filename, opts.at(i));
}

// Reading from csv filename with some entities with the name "[file].csv", 
// replace those with "[file]_suffix.csv" and create a new csv with name 
// filename_out.
void copy_opts_csv(std::string &filename, std::string &filename_out, 
                                                    std::string &suffix) {
  std::vector<std::vector<std::string> > opts(0, std::vector<std::string> (0,""));
  ReadNames(filename, ";", opts);
  add_suffix_csv(opts, suffix);
  create_opts_csv(filename_out, opts);
}

// Conjugate gradient based surface minimization. 
// Sequentilly minimize surfaces and improve element quality by moving interior 
// vertices until the surface energy is arbitrarily close to a minimum.
class vd_evolve_conj_grad{
  private:
  public:
    apf::Mesh2* mesh;
    cell_base* c_base;
    field_calc* f_calc;
    bool ext_cell;
    bool upd_shell;

    // Boundary and interior vertices:
    std::vector<apf::MeshEntity*> v_b;
    std::vector<apf::MeshEntity*> v_i;

    // All tetrahedra:
    std::vector<apf::MeshEntity*> et_3c;
    // Tetrahedra and triangles adjacent to boundary vertices:
    std::vector<apf::MeshEntity*> et_b;
    std::vector<apf::MeshEntity*> es_b;

    std::vector<apf::MeshEntity*> ee;
    std::vector<apf::MeshEntity*> es;
    std::vector<apf::MeshEntity*> et;

    std::vector<std::vector<apf::MeshEntity*> > tet_v;
    std::vector<std::vector<apf::MeshEntity*> > tri_v;

    std::map<apf::MeshEntity*, apf::Vector3> pos_ori;
    std::map<apf::MeshEntity*, apf::Vector3> v_pos;

    apf::Field* f_field;

    // Most negative volume, used to scale w, st. w*vol_t/vol_scale = 100;
    double vol_neg;
    double v_th_tol;
    double g_tol;
    double g_max;
    double epsilon;

    // Conjugate gradient related:
    std::vector<apf::Vector3> x0;
    std::vector<apf::Vector3> x1;

    std::vector<apf::Vector3> g0;
    std::vector<apf::Vector3> g1;
    std::vector<apf::Vector3> dir;
    std::vector<apf::Vector3> ndir;

    bool inv_flag;

    int iter_limit;
    int iter;
    double g_norm;
    double xa, xb, xc, xd, xe, xu;
    double fa, fb, fc, fd, fe, fu;
    double p, q, r, s, t, m, tol, tol2;
    bool inv_quad_step;
    double Phi0, Phi1;

    // Destructor:
    vd_evolve_conj_grad(apf::Mesh2* m_in, cell_base* cb_in, field_calc* f_calc_in,
                             bool ext_cell_in, 
                             bool upd_shell_in, 
                             double g_tol_in = 0.00001, 
                             double v_th_tol_in = 0.00001);

    void collect_ent();
    bool update_vol_inv_grad();

    // Calculate the potential.
    double calc_energy();
    void calc_Cs();
    void upd_grad(std::vector<apf::Vector3> &g_curr);

    void shift_v_sp_pos(std::vector<apf::Vector3> &x_curr,
                        double x_in);
    // Revert the split vertices to their original positions.
    void revert_ori();

    // Relax the interior vertices to prevent element inversions.
    bool relax_int();
    // Try to find a hull for non-inverting tets by conjugate gradient method.
    double relax();

    ~vd_evolve_conj_grad();
};

void vd_evolve_conj_grad::collect_ent() {
  vd_entlist e_list(mesh, c_base);

  int sz = e_list.e.at(3).size();
  int count = 0;
  for(int i = 0; i < sz; i++) {
    count = count +  e_list.e.at(3).at(i).at(0).size();
  }
  v_i.clear();
  v_i.reserve(count);

  for(int i = 0; i < sz; i++) {
    v_i.insert(v_i.end(), e_list.e.at(3).at(i).at(0).begin(), 
                                             e_list.e.at(3).at(i).at(0).end());
  }

  count = 0;
  for(int i = 0; i < sz; i++) {
    count = count +  e_list.e.at(3).at(i).at(3).size();
  }
  et_3c.clear();
  et_3c.reserve(count);

  for(int i = 0; i < sz; i++) {
    et_3c.insert(et_3c.end(), e_list.e.at(3).at(i).at(3).begin(), 
                                             e_list.e.at(3).at(i).at(3).end());
  }

  count = 0;
  for(int dim = 1; dim < 3; dim++) {
    sz = e_list.e.at(dim).size();
    for(int i = 0; i < sz; i++) {
      count = count +  e_list.e.at(dim).at(i).at(0).size();
    }
  }
  if(upd_shell) {
    ext_shell* e_sh = f_calc->get_e_sh();
    sz = e_list.e.at(0).size();
    for(int i = 0; i < sz; i++) {
      if(c_base->get_cell_ext(0, i)) {
        count = count +  e_list.e.at(0).at(i).at(0).size();
      }
    }
  }

  v_b.clear();
  v_b.reserve(count);

  for(int dim = 1; dim < 3; dim++) {
    sz = e_list.e.at(dim).size();
    for(int i = 0; i < sz; i++) {
      v_b.insert(v_b.end(), e_list.e.at(dim).at(i).at(0).begin(), 
                                           e_list.e.at(dim).at(i).at(0).end());
    }
  }

  if(upd_shell) {
    ext_shell* e_sh = f_calc->get_e_sh();
    sz = e_list.e.at(0).size();
    for(int i = 0; i < sz; i++) {
      if(c_base->get_cell_ext(0, i)) {
        v_b.insert(v_b.end(), e_list.e.at(0).at(i).at(0).begin(), 
                                             e_list.e.at(0).at(i).at(0).end());
      }
    }
  }

  x0.resize(v_b.size());
  x1.resize(v_b.size());
  g0.resize(v_b.size());
  g1.resize(v_b.size());
  dir.resize(v_b.size());
  ndir.resize(v_b.size());

  vd_set_up(mesh, &v_b, &ee);
  vd_set_up(mesh, &ee, &es_b);
  vd_set_up(mesh, &es_b, &et_b);

  tet_v.resize(et_b.size());
  for(int i = 0; i < et_b.size(); i++) {
    tet_v.at(i).resize(4);
  }

  tri_v.resize(es_b.size());
  for(int i = 0; i < es_b.size(); i++) {
    tri_v.at(i).resize(3);
  }
}

vd_evolve_conj_grad::vd_evolve_conj_grad(apf::Mesh2* m_in, cell_base* cb_in, 
                             field_calc* f_calc_in, 
                             bool ext_cell_in, 
                             bool upd_shell_in,
                             double g_tol_in, 
                             double v_th_tol_in) : 
                  // Initialization of mesh related information 
                  mesh(m_in), c_base(cb_in),
                  f_calc(f_calc_in), ext_cell(ext_cell_in),
                  upd_shell(upd_shell_in),
                  // Initialization of entity vectors and maps 
                  ee(0), es(0), et(0),
                  v_b(0), v_i(0),
                  et_3c(0), et_b(0), es_b(0),

                  pos_ori{}, v_pos{},
                  tet_v (0, std::vector<apf::MeshEntity*>(4)),
                  tri_v (0, std::vector<apf::MeshEntity*>(3)),

                  f_field(NULL),
                  // Conjugate gradient initialization
                  x0(0, apf::Vector3(0,0,0)),
                  x1(0, apf::Vector3(0,0,0)),
                  g0(0, apf::Vector3(0,0,0)),
                  g1(0, apf::Vector3(0,0,0)),
                  dir(0, apf::Vector3(0,0,0)),
                  ndir(0, apf::Vector3(0,0,0)),
                  vol_neg(0),
                  epsilon(std::numeric_limits<double>::min()*10e10),
                  v_th_tol(v_th_tol_in), g_tol(g_tol_in),
                  g_max(0), 

                  inv_flag(false),
                  iter_limit(400),
                  iter(0), g_norm(0), 
                  xa(0), xb(0), xc(0), xd(0), xe(0), xu(0),
                  fa(0), fb(0), fc(0), fd(0), fe(0), fu(0),
                  p(0), q(0), r(0), s(0), t(0), m(0), tol(0), tol2(0),
                  inv_quad_step(false),
                  Phi0(0), Phi1(0) {

  f_field = vd_att_vv_field(mesh, "f_restore");
  /////////////////////////////////////////////
  // Collect adjacency related information:  //
  /////////////////////////////////////////////
  collect_ent();

  apf::Downward d_v;
  apf::Downward d_e;
  apf::Vector3 temp(0,0,0);
  std::vector<apf::Vector3> pts (4, apf::Vector3(0,0,0));
  std::vector<apf::MeshEntity*>::iterator it;

  // Reserve the vertex-tetrahedra lists:
  for(int i = 0; i < v_b.size(); i++) {
    // In case exterior, there could be one more vertex:
    mesh->getPoint(v_b.at(i), 0, temp);
    pos_ori[v_b.at(i)] = temp;
    x0.at(i) = temp;
  }

  for (int i = 0; i < et_b.size(); i++) {
    mesh->getDownward(et_b.at(i), 0, d_v);
    for(int j = 0; j < 4; j++) {
      tet_v.at(i).at(j) = d_v[j];
      mesh->getPoint(d_v[j], 0, temp);
      v_pos[d_v[j]] = temp;
    }
  }

  for (int i = 0; i < es_b.size(); i++) {
    mesh->getDownward(es_b.at(i), 0, d_v);
    for(int j = 0; j < 3; j++) {
      tri_v.at(i).at(j) = d_v[j];
    }
  }

  apf::writeVtkFiles("./output/Glens_inv", mesh);
}

bool vd_evolve_conj_grad::update_vol_inv_grad() {
  double min_vol = epsilon;
  std::vector<apf::Vector3> pts(4, apf::Vector3(0,0,0));

  inv_flag = false;
  vol_neg = 0.;

  for (int i = 0; i < et_b.size(); i++) {
    for(int j = 0; j < 4; j++) {
      pts.at(j) = v_pos[tet_v.at(i).at(j)];
    }
    double vol_curr = vd_volume_tet(&pts);
    if(vol_curr < std::numeric_limits<double>::min() and
       vol_curr < vol_neg)
      vol_neg = vol_curr;

    if(vol_curr < min_vol) {
      inv_flag = true;
    }
    assert(!std::isnan(vol_curr));
  }
  return inv_flag;
}

// Calculate the potential.
double vd_evolve_conj_grad::calc_energy() {
  std::vector<apf::Vector3> pts(3, apf::Vector3(0,0,0));

  double Phi = 0;
  for (int i = 0; i < es_b.size(); i++) {
    for(int j = 0; j < 3; j++) {
      pts.at(j) = v_pos[tri_v.at(i).at(j)];
    }
    double gam_curr = f_calc->gam2(mesh, es_b.at(i));
    Phi = Phi + vd_area_out_n(pts.at(0), pts.at(1), pts.at(2)).getLength()
                                                                    *gam_curr;
  }

  return Phi;
}

// TODO this works well enough, but it doesn't consider the cross terms
// grad_i phi_l and only considers grad_i phi_i. See SM of the method paper.
void vd_evolve_conj_grad::upd_grad(std::vector<apf::Vector3> &g_curr) {

  for (int i = 0; i < v_b.size(); i++) {
    g_curr.at(i) = f_calc->vd_calc_force(mesh, v_b.at(i))*(-1);
    assert(!std::isinf( g_curr.at(i).getLength() ));

    apf::setVector(f_field, v_b.at(i), 0, g_curr.at(i));
    //std::cout << " f_i: " << g_curr << std::endl;
  }
}

void vd_evolve_conj_grad::shift_v_sp_pos(std::vector<apf::Vector3> &x_curr,
                    double x_in) {
  apf::Vector3 pos(0,0,0);
  //apf::Vector3 temp(0,0,0);
  for(int i = 0; i < v_b.size(); i++) {
    if (f_calc->chk_vert_special(mesh, v_b.at(i)))
      pos = f_calc->get_vec_special(mesh, v_b.at(i), ndir.at(i)*x_in);
    else
      pos = ndir.at(i)*x_in;

    assert(!std::isnan(pos.getLength()));
    pos = x_curr.at(i) + pos;
    assert(!std::isnan(pos.getLength()));
    mesh->setPoint(v_b.at(i), 0, pos);
    v_pos[v_b.at(i)] = pos;
  }
  update_vol_inv_grad();
}

// Revert the split vertices to their original positions.
void vd_evolve_conj_grad::revert_ori() {
  for (int i = 0; i < v_b.size(); i++) {
    mesh->setPoint(v_b.at(i), 0, pos_ori[v_b.at(i)]);
  }
}

bool vd_evolve_conj_grad::relax_int() {

  std::map<apf::MeshEntity*, int> e_map{};
  std::vector<std::vector<apf::MeshEntity*> > ent(4, 
                        std::vector<apf::MeshEntity*>(0));

  double r_edge_min = 0.01;
  apf::Vector3 midpoint(0,0,0);

  vd_relax_cont rc(mesh, c_base, f_calc, 0, 0, 
                         v_i, ent, et_3c, e_map,
                         false, false, midpoint, r_edge_min, 
                                            v_th_tol);

  return rc.relax();
}

double vd_evolve_conj_grad::relax() {
  double STEP_FIRST = 10e-2;

  Phi0 = calc_energy();
  upd_grad(g0);

  for (int i = 0; i < v_b.size(); i++) {
    dir.at(i) = g0.at(i)*(-1);
  }
  g_norm = 0;
  g_max = 0;
  for (int i = 0; i < v_b.size(); i++) {
    double g_curr = dir.at(i) * dir.at(i);
    g_norm = g_norm + g_curr;
    if(g_curr > g_max)
      g_max = g_curr;
  }
  g_norm = std::sqrt(g_norm);

  for (int i = 0; i < v_b.size(); i++) {
    ndir.at(i) = dir.at(i) / g_norm;
  }

  for (int b = 0; b < MAX_CG_ITER; ++b) {

    inv_flag = update_vol_inv_grad();
    ////////////////////////////////////////////////////////////
    // bounds the line search, with ax < bx < cx and fa > fb < fc
    ////////////////////////////////////////////////////////////
    if(inv_flag) {
      if(!relax_int())
        break;

      Phi0 = calc_energy();
      upd_grad(g0);

      for (int i = 0; i < v_b.size(); i++) {
        dir.at(i) = g0.at(i)*(-1);
      }
      g_norm = 0;
      g_max = 0;
      for (int i = 0; i < v_b.size(); i++) {
        double g_curr = dir.at(i) * dir.at(i);
        g_norm = g_norm + g_curr;
        if(g_curr > g_max)
          g_max = g_curr;
      }
      g_norm = std::sqrt(g_norm);

      if(std::fabs(g_norm) < std::numeric_limits<double>::min())
        g_norm = 1;

      for (int i = 0; i < v_b.size(); i++) {
        ndir.at(i) = dir.at(i) / g_norm;
      }
    }

    xa = 0.;
    fa = Phi0;

    xb = STEP_FIRST;
    shift_v_sp_pos(x0, xb);
    fb = calc_energy();

    while (fb > fa && xb > EPS) {
      // decrease step until energy is decreasing along ndir
      xb /= 10.;
      shift_v_sp_pos(x0, xb);
      fb = calc_energy();
    }

    // try inverse quadratic interpolation
    p = 0;
    for (int i = 0; i < v_b.size(); i++) {
      p = p - g0.at(i)*ndir.at(i);
    }
    p = p * xb;
    q = (fb - fa) + p;
    inv_quad_step = false;
    if (q > EPS) {
      // parabola is concave up, find the minimum
      xc = (p * xb) / (2. * q);
      if (xc > (MAX_BND_STEP + 1.) * xb) {
        // maximum step length
        inv_quad_step = true;
        // MAX_BND_STEP could be too large when xb is not reduced.
        //if (xb/STEP_FIRST < 0.01)
        //  xc = (100 + 1.) * xb;
        //else
        //  xc = (MAX_BND_STEP + 1.) * xb;
        xc = (MAX_BND_STEP + 1.) * xb;

        shift_v_sp_pos(x0, xc);
        fc = calc_energy();
      } else if (xc > (BND_MAG + 1.) * xb) {
        // normal step
        inv_quad_step = true;
        shift_v_sp_pos(x0, xc);
        fc = calc_energy();
        if (fc < fb) {
          // try to step past minimum
          SHIFT2(xa, xb, xc);
          xc = xb + BND_MAG * (xb - xa);
          SHIFT2(fa, fb, fc);
          shift_v_sp_pos(x0, xc);
          fc = calc_energy();
        }
      } else if (xc > xa + SQRT_EPS && xc < xb - SQRT_EPS) {
        // minimum falls in (ax, bx)
        shift_v_sp_pos(x0, xc);
        fc = calc_energy();
        if (fc < fb) {
          // found bracket, all done
          inv_quad_step = true; 
          std::swap(xb, xc);
          std::swap(fb, fc);
        }
      }
    }
    if (!inv_quad_step) {
      // quadratic interpolation failed, conservative step
      xc = (BND_MAG + 1.) * xb;
      shift_v_sp_pos(x0, xc);
      fc = calc_energy();
    }
    while (fc < fb) {
      // try inverse quadratic interpolation
      p = xc - xb;
      q = xa - xb;
      r = (fa - fb) * p;
      s = (fc - fb) * q;
      t = r - s;
      inv_quad_step = false;
      if (t > EPS) {
        // parabola is concave up, find the minimum
        xd = xb + (r * p - s * q) / (2. * t);
        if (xd > xc + MAX_BND_STEP * p) {
          // maximum step length
          inv_quad_step = true;
          xd = xc + MAX_BND_STEP * p;

          shift_v_sp_pos(x0, xd);
          fd = calc_energy();
        } else if (xd > xc + BND_MAG * p) {
          // normal step
          inv_quad_step = true;
          shift_v_sp_pos(x0, xd);
          fd = calc_energy();
          if (fd < fc) {
            // try to step past minimum
            SHIFT3(xa, xb, xc, xd);
            xd = xc + BND_MAG * (xc - xb);
            SHIFT3(fa, fb, fc, fd);
            shift_v_sp_pos(x0, xd);
            fd = calc_energy();
          }
        } else if (xd > xb + SQRT_EPS && xd < xc - SQRT_EPS) {
          // minimum falls in (ax, bx)
          shift_v_sp_pos(x0, xd);
          fd = calc_energy();
          if (fd < fc) {
            // found bracket, all done
            inv_quad_step = true;
            SHIFT2(xa, xb, xd);
            SHIFT2(fa, fb, fd);
            break;
          }
        }
      }
      if (!inv_quad_step) {
        // quadratic interpolation failed, conservative step
        xd = xc + BND_MAG * p;
        shift_v_sp_pos(x0, xd);
        fd = calc_energy();
      }
      // bookkeeping for next iteration
      SHIFT3(xa, xb, xc, xd);
      SHIFT3(fa, fb, fc, fd);
    }

    ////////////////////////////////////////////////////////////
    // Brent's method to find minimum along search direction, a translation
    // of the ALGOL 60 algorithm on page 79 of R. P. Brent, Algorithms for
    // Minimization Without Derivatives, 1973 with minor modifications. The
    // author gave permission to use this algorithm in private communication.
    ////////////////////////////////////////////////////////////
    // use values from bounding the line search
    if (fc < fa) {
      xd = xc;
      xe = xa;
      fd = fc;
      fe = fa;
    } else {
      xd = xa;
      xe = xc;
      fd = fa;
      fe = fc;
    }
    t = (xb < 0.5 * (xa + xc) ? xc : xa) - xb;
    s = PHI_SQ_INV * t;
    int c = 0;
    for (c = 0; c < MAX_LS_ITER; ++c) {
      m = 0.5 * (xa + xc);
      tol  = SQRT_EPS * (fabs(xb) + 1.);
      tol2 = 2. * tol;
      // check stopping criterion
      if (fabs(xb - m) > tol2 - 0.5 * (xc - xa)) {
        inv_quad_step = false;
        if (fabs(t) > tol) {
          // inverse quadratic interpolation
          p = (xb - xd) * (fb - fe);
          q = (xb - xe) * (fb - fd);
          r = (xb - xe) * q - (xb - xd) * p;
          q = 2. * (q - p);
          if (q > 0.) { r = -r; }
          q = fabs(q);
          SHIFT2(p, t, s);
          // mistake in ALGOL 60 routine, second condition is inverted
          if (fabs(r) < fabs(0.5 * q * p) && r > q * (xa - xb) && r < q * (xc - xb)) {
            // take inverse quadratic interpolation step
            inv_quad_step = true;
            s = r / q;
            xu = xb + s;
            // f should not be evaluated too close to xa or xc
            if (xu - xa < tol2 || xc - xu < tol2) {
              s = (xb < m ? tol : -tol);
            }
          }
        }
        if (!inv_quad_step) {
          // interpolation failed, take golden section step
          t = (xb < m ? xc : xa) - xb;
          s = PHI_SQ_INV * t;
        }

        // f should not be evaluated too close to xb
        xu = xb + (fabs(s) >= tol ? s : (s > 0. ? tol : -tol));

        shift_v_sp_pos(x0, xu);
        fu = calc_energy();
        // bookkeeping for next iteration
        if (fu <= fb) {
          if (xu < xb) { xc = xb; } else { xa = xb; }
          SHIFT3(xe, xd, xb, xu);
          SHIFT3(fe, fd, fb, fu);
        } else {
          if (xu < xb) { xa = xu; } else { xc = xu; }
          if (fu <= fd || xd == xb) {
            SHIFT2(xe, xd, xu);
            SHIFT2(fe, fd, fu);
          } else if (fu <= fe || xe == xb || xe == xd) {
            xe = xu;
            fe = fu;
          }
        }
      }
      else {
        for (int i = 0; i < v_b.size(); i++) {
          x1.at(i) = x0.at(i) + ndir.at(i) * xb;
        }

        // found minimum, apply change and update energy
        Phi1 = fb;
        break;
      }
    }
    if (c == MAX_LS_ITER)
      std::cerr << "RunSimulation: max LS iteration reached" << std::endl;

    ////////////////////////////////////////////////////////////
    // Conjugate gradient
    ////////////////////////////////////////////////////////////
    // check energy convergence
    //if (fabs(Phi1 - Phi0) < E_TOL) {
    //  x0 = x1;
    //}

    // check gradient convergence
    shift_v_sp_pos(x1, 0);

    upd_grad(g1);

    // Direction update given by Y.H Dai, C.X. Kou, SIAM Journal of
    // Optimization, v 23, p 296-320, 2013
    for (int i = 0; i < v_b.size(); i++) {
      g0.at(i) = g0.at(i) - g1.at(i);
    }
    for (int i = 0; i < v_b.size(); i++) {
      x0.at(i) = x0.at(i) - x1.at(i);
    }

    double B0 = 0;
    double B1 = 0;
    double B0_a = 0;
    double B0_b = 0;
    double B0_c = 0;

    double B1_a = 0;
    double B1_b = 0;
    for (int i = 0; i < v_b.size(); i++) {
      B0_a = B0_a + g0.at(i)*g0.at(i);
      B0_b = B0_a + g1.at(i)*x0.at(i);
      B0_c = B0_a + g0.at(i)*x0.at(i);
    }
    if(std::fabs(B0_c) < std::numeric_limits<double>::min())
      B0_c = 1;
    B0 = B0_a*B0_b/B0_c;

    for (int i = 0; i < v_b.size(); i++) {
      B0_a = B0_a + g0.at(i)*g1.at(i);
      B0_b = B0_a + dir.at(i)*g0.at(i);

      B1_a = B0_a + g1.at(i)*dir.at(i);
      B1_b = B0_a + dir.at(i)*dir.at(i);
    }
    if(std::fabs(B0_b) < std::numeric_limits<double>::min())
      B0_b = 1;

    B0 = (B0_a - B0)/B0_b;

    if(std::fabs(B1_b) < std::numeric_limits<double>::min())
      B1_b = 1;

    B1 = B1_a / (2. * B1_b);

    double mult = fmax(B0, B1);
    for (int i = 0; i < v_b.size(); i++) {
      dir.at(i) = dir.at(i) * mult - g1.at(i);
    }

    g_norm = 0;
    g_max = 0;
    for (int i = 0; i < v_b.size(); i++) {
      double g_curr = dir.at(i) * dir.at(i);
      g_norm = g_norm + g_curr;
      if(g_curr > g_max)
        g_max = g_curr;
    }
    g_norm = std::sqrt(g_norm);

    if(std::fabs(g_norm) < std::numeric_limits<double>::min())
      g_norm = 1;

    for (int i = 0; i < v_b.size(); i++) {
      ndir.at(i) = dir.at(i) / g_norm;
    }

    // prepare for next iteration
    for (int i = 0; i < v_b.size(); i++) {
      x0.at(i) = x1.at(i);
    }
    shift_v_sp_pos(x0, 0);

    Phi0 = Phi1;
    for (int i = 0; i < v_b.size(); i++) {
      g0.at(i) = g1.at(i);
    }

    inv_flag = false;
    //v_th = v_th_tol*vol_min/vol_scale;
    //vol_scale = 0;
    inv_flag = update_vol_inv_grad();

    apf::writeVtkFiles("./output/Glens_step", mesh);

    if(g_max < g_tol) {
      std::cout<< "Maximum force magnitude below threshold!" << std::endl;
      b = MAX_CG_ITER;
    }
    else
      std::cout<< "Maximum force magnitude:" << g_max << std::endl;
  }

  // End conjugate gradient method
  return g_max;
}

vd_evolve_conj_grad::~vd_evolve_conj_grad() {
}




// Update the list of interior vertices. Return the v_th_vol, ratio of the goal
// volume to the smallest volume.
double collect_relax(apf::Mesh2* m, vd_entlist* e_list, 
                                       std::vector<apf::MeshEntity*> &v_3c,
                                       std::vector<apf::MeshEntity*> &et_3c,
                                       std::map<apf::MeshEntity*, int> &e_map) {
  e_map.clear();
  int sz = e_list->e.at(3).size();
  int count = 0;
  for(int i = 0; i < sz; i++) {
    count = count +  e_list->e.at(3).at(i).at(0).size();
  }
  v_3c.clear();
  v_3c.reserve(count);

  for(int i = 0; i < sz; i++) {
    v_3c.insert(v_3c.end(), e_list->e.at(3).at(i).at(0).begin(), 
                                             e_list->e.at(3).at(i).at(0).end());
  }

  count = 0;
  for(int i = 0; i < sz; i++) {
    count = count +  e_list->e.at(3).at(i).at(3).size();
  }
  et_3c.clear();
  et_3c.reserve(count);

  for(int i = 0; i < sz; i++) {
    et_3c.insert(et_3c.end(), e_list->e.at(3).at(i).at(3).begin(), 
                                             e_list->e.at(3).at(i).at(3).end());
  }

  double v_th_tol = 1;
  // v_th = std::fabs(v_th_tol*vol_min/vol_scale);
  // v_th_tol = vol_goal/vol_min;
  if(et_3c.size() > 0) {
    double vol_min = vd_volume_tet(m, et_3c.at(0));
    double vol_tot = vol_min;
    for(int i = 1; i < et_3c.size(); i++) {
      double vol_curr = vd_volume_tet(m, et_3c.at(i));
      vol_tot = vol_tot + vol_curr;
      if(vol_curr < vol_min)
        vol_min = vol_curr;
    }
    // For feasibility the goal minimum volume should be less than the average 
    // volume of a tetrahedron. Underestimating the volume of the grain 1, divided
    // by number of tets of grain 1 achieves that.
    v_th_tol = vol_tot/15/e_list->e.at(3).at(0).at(3).size()/vol_min;
  }
  return v_th_tol;
}

void split_interior_edges(apf::Mesh2* m, field_calc* f_calc, vd_sim &sim_trial) {
  f_calc->vd_del_fields(m);

  Linear sf(m, sim_trial.get_adapt_ln());
  ma::Input* in = ma::configure(m, &sf);
  repl_sz_field(in, m, MA_SIZE_TYPE::EDGE_REFINER);
  ModelEdgeRefiner* ref = (ModelEdgeRefiner* ) in->sizeField;
  ref->split_all = true;

  in->shouldRunPreZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostParma = false;
  in->shouldRefineLayer = true;

  in->shouldCoarsen = false;
  in->shouldCoarsenLayer = false;

  in->shouldFixShape = false;
  in->maximumEdgeRatio = 100;
  in->maximumIterations = 1;

  ma::adapt(in);
  f_calc->vd_att_fields(m);

}
/*
void evolve_interior(apf::Mesh2* m, field_calc* f_calc, vd_elist* e_list) {
  int sz = e_list.e.at(3).size();
  for(int i = 0; i < sz; i++) {
    for(int j = 0; j < ; j++) {
    if(! and c_base->get_cell_ext(2, i)) {
      f_calc->upd_gam2(m, e_list->e.at(2).at(i).at(2), 0);
    }
  }
}
*/
// Return true if the relaxation was not successful.
bool evolve_interior(vd_sim &sim_trial, vd_eqn_kuprat_volm& eqn_kup, double dt_set) {
  cell_base* c_base = sim_trial.get_c_base();
  apf::Mesh2* m = sim_trial.get_mesh();
  field_calc* f_calc = sim_trial.get_f_calc();
  vd_entlist* e_list = sim_trial.get_elist();

  double sub_param = 1;
  f_calc->set_drag_glob(getAverageEntSize(m, 2));
  int iter_sub = 10;
  bool sub_suc = false;

  for (int i = 0; i < iter_sub; i++) {
    printf("Interior trial %d.\n", i);

    eqn_kup.calc_vel();

    std::stringstream ss;
    ss <<"output/before_trial" << i;

    std::string tmp = ss.str();
    const char* cstr = tmp.c_str();

    sim_trial.save_vtk_name(cstr);

    double t_sub = vd_find_min_t_all(m, false, -1)/5;
    f_calc->vdparam.adj_dt(std::min(t_sub/sub_param, dt_set));
    if(f_calc->vdparam.dt < std::numeric_limits<double>::min())
      f_calc->vdparam.adj_dt(dt_set);
    printf("dt: %.8f, dt_set: %.8f, t_sub: %.8f.\n", f_calc->vdparam.dt, 
                                                    dt_set, t_sub);

    vd_apply_vel_field(m, f_calc->vdparam.dt);

    m->verify();
    ss.str("");

    ss.clear();
    ss <<"output/after_trial" << i;
    tmp = ss.str();
    cstr = tmp.c_str();

    sim_trial.save_vtk_name(cstr);

    //PCU_Thrd_Barrier();
    bool chk_ma = sim_trial.check_ma_sgn();
    //PCU_Thrd_Barrier();

    if (PCU_Or(chk_ma)) {
      vd_apply_vel_field(m, -f_calc->vdparam.dt);
      sub_param = sub_param*2;
      //PCU_Thrd_Barrier();
    }
    else {
      sub_suc = true;
    }
  }
  if(sub_suc)
    return false;
  return true;
}

void opts2vecintpair(std::vector<std::string>& opts, 
                                      std::vector<std::pair<int,int>>& vipair) {

  vipair.clear();
  vipair.reserve(opts.size());

  for(int i = 0; i < opts.size(); i++) {
    std::string temp("");
    std::string::size_type sz;
    temp = opts.at(i);

    int c_dim = std::stoi(temp, &sz);
    temp = temp.substr(sz);
    int c_id = std::stoi(temp, &sz);
    vipair.push_back(std::make_pair(c_dim, c_id));
  }
}

///////////////////////////////////////////////////

int main(int argc, char** argv)
{

  assert(argc==6);
  const char* modelFile = "./mshfiles/troctahedra_unit.tess";
  std::vector<std::pair<int,int> > cells_ref (0, std::make_pair(0,0));
  std::vector<std::vector<std::string> > temp_str_vec(0, std::vector<std::string> (0,""));
  ReadNames("./mshfiles/ref_list_gauss.txt", ";", temp_str_vec);
  opts2vecintpair(temp_str_vec.at(0), cells_ref);


  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();

  // th to stop, FileName of the tripline mesh, and the FileSphere for the sphere
  // to use for edge length for refinement.
  const char* FileName = argv[1];
  double max_ref = atof(argv[2]);
  double mult = atof(argv[3]);
  double dt = atof(argv[4]);
  double g_tol = atof(argv[5]);

  double len_basis = 1;

  double ref_curr = 2.;

  double coarse_th = 0.9;
  double split_th = 1.1;
  int ad_nbr = 6;

  std::string temp("");
  //temp = temp + "_" + std::to_string(ref_curr);
  std::string opts("");

  std::vector<std::vector<std::string> > ext_opts(3, std::vector<std::string> (0,""));
  ext_opts.at(0).reserve(11);
  opts = "GAUSS";
  ext_opts.at(0).push_back(opts);
  opts = "./output/GAUSS" + temp + ".csv";
  ext_opts.at(0).push_back(opts);
  opts = "./output/EULERMOD" + temp + ".csv";
  ext_opts.at(0).push_back(opts);
  opts = "./output/EULERERROR" + temp + ".csv";
  ext_opts.at(0).push_back(opts);
  opts = "./output/0CDEV" + temp + ".csv";
  ext_opts.at(0).push_back(opts);
  opts = "./output/1CLEN" + temp + ".csv";
  ext_opts.at(0).push_back(opts);
  opts = "./output/3CVOL" + temp + ".csv";
  ext_opts.at(0).push_back(opts);
  opts = "./output/0C1C" + temp + ".csv";
  ext_opts.at(0).push_back(opts);
  opts = "./output/0C3C" + temp + ".csv";
  ext_opts.at(0).push_back(opts);
  opts = "./output/0CVEL" + temp + ".csv";
  ext_opts.at(0).push_back(opts);
  opts = "NEUMANN_EXT";
  ext_opts.at(0).push_back(opts);

  ext_opts.at(1).reserve(4);
  opts = "MS";
  ext_opts.at(1).push_back(opts);
  opts = "./output/MS" + temp + ".csv";
  ext_opts.at(1).push_back(opts);
  opts = "./output/ROC" + temp + ".csv";
  ext_opts.at(1).push_back(opts);
  opts = "./output/VOL" + temp + ".csv";
  ext_opts.at(1).push_back(opts);

  ext_opts.at(2).reserve(16);
  opts = "MEAS";
  ext_opts.at(2).push_back(opts);
  opts = "./output/AREA" + temp + ".csv";
  ext_opts.at(2).push_back(opts);
  for(int i = 0; i < 14; i++) {
    opts = "2 " + std::to_string(i+1);
    ext_opts.at(2).push_back(opts);
  }

  std::string ext_opts_file = std::string("./output/ext_opts_temp.txt");
  create_opts_csv(ext_opts_file, ext_opts);

  vd_sim sim_trial;
  sim_trial.set_ins_flag(false);
  sim_trial.set_col_flag(false);

  sim_trial.set_adapt(1);
  sim_trial.set_adapt_param(ref_curr);
  sim_trial.set_ad_th(coarse_th, split_th);

  //sim_trial.set_adapt_type(ADAPT_TYPE::ADAPT_STEP_1CELL);
  sim_trial.set_adapt_type(ADAPT_TYPE::ADAPT_BOUND);

  sim_trial.set_mesh(modelFile);

  sim_trial.set_field_calc(VEL_TYPE::MASON_MIR);
  //sim_trial.set_vec_sp_calc(PROJ_TYPE::FIXED);
  sim_trial.set_vec_sp_calc(PROJ_TYPE::EXT_SHELL);
  sim_trial.set_ext_gam_flag(EXT_GAM_TYPE::ZERO);
  sim_trial.set_integ_type(INTEG_TYPE::RK2);
  sim_trial.set_mov_flag(false);

  cell_base* c_base = sim_trial.get_c_base();
  apf::Mesh2* m = sim_trial.get_mesh();
  field_calc* f_calc = sim_trial.get_f_calc();
  vd_entlist* e_list = sim_trial.get_elist();


  //ent_conn e_con2c;
  //ent_conn e_con1c;
  //std::vector<int> ref_2c({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,    
  //                         20,24,25,26,27,28,29,33,34,39,40,41,42,43,44,48,
  //                         49,50,51,52,53,57,58,63,64,69,70,74,75,80,87});
  //e_con2c = ref_2c;

  for(int i = 0; i < c_base->get_sz(2); i++) {
    if(c_base->get_cell_ext(2, i))
      f_calc->b_en.c2s_zeros[i+1] = true;
  }

  ent_conn* e_con0 = new ent_conn();
  c_base->get_conn_dim_gmi(0, 3, 1, e_con0);

  f_calc->set_drag_rat(50);

  f_calc->upd_d2(m);
  f_calc->upd_gam2(m);
  f_calc->vd_att_fields(m);
  f_calc->vd_calc_vel(m);
  double v_norm_st = f_calc->vd_get_max_vel_norm_int(m);

  f_calc->vd_calc_vel(m);

  sim_trial.set_time(1.);
  sim_trial.set_extract(true, ext_opts_file);
  sim_trial.extract_data();
  sim_trial.set_extract(false, ext_opts_file);

  std::stringstream ss;
  ss << FileName << 0.;
  temp = ss.str();
  std::string tmp = temp + ".smb";

  const char* cstr = tmp.c_str();
  m->writeNative(cstr);

  tmp = temp + ".dmg";
  cstr = tmp.c_str();
  c_base->vd_write_dmg(cstr);

  sim_trial.set_time(0/10000);
  sim_trial.set_dt(dt/10000, true);
  sim_trial.set_end(50.*dt/10000);

  sim_trial.set_adapt_bound_cells(1./ref_curr, cells_ref);
  //sim_trial.set_adapt_param(ref_curr);
  for(int i = 0; i < ad_nbr; i++)
    sim_trial.adapt();

  split_interior_edges(m, f_calc, sim_trial);

  sim_trial.refresh_elist();
  sim_trial.save_vtk_name("./output/before");

  f_calc->upd_d2(m);
  f_calc->upd_gam2(m);
  f_calc->vd_att_fields(m);

  vd_eqn_kuprat_volm eqn_kup(m, c_base, f_calc, e_list);
  evolve_interior(sim_trial, eqn_kup, dt/1000.);

  f_calc->vd_calc_vel(m);

  bool adapt_flag = true;
  double dt_curr = dt;
  double t_rat = 10.;

  apf::Vector3 midpoint(0,0,0);
  double r_edge_min = 0.01;
  std::vector<apf::MeshEntity*> v_3c(0);
  std::vector<apf::MeshEntity*> et_3c(0);
  std::vector<std::vector<apf::MeshEntity*> > ent(4, 
                        std::vector<apf::MeshEntity*>(0));
  std::map<apf::MeshEntity*, int> e_map{};

  // The ratio of the smallest volume to the goal volume for the interior tets.
  double v_th_tol = collect_relax(m, e_list, v_3c, et_3c, e_map);

  while (ref_curr < max_ref) {
    double v_norm = f_calc->vd_get_max_vel_norm_int(m);
    std::cout << "v_norm " << v_norm << std::endl;

    v_th_tol = collect_relax(m, e_list, v_3c, et_3c, e_map);
    vd_att_vs_field(m, "cvx");
    vd_att_vs_field(m, "v_sp");
    vd_att_vv_field(m, "f_restore");
    vd_evolve_conj_grad rc(m, c_base, f_calc, true, true,
                                            g_tol, v_th_tol);
    f_calc->set_drag_glob(getAverageEntSize(m, 2));
    while (v_norm > g_tol) {
      double time_curr = sim_trial.get_time();

      int iter_sub = 5;
      for (int i = 0; i < iter_sub; i++) {
        printf("Trial %d.\n", i);
        v_norm = rc.relax();

        m->verify();

        ss.clear();
        ss.str("");
        ss <<"output/after_trial" << i;
        tmp = ss.str();
        const char* cstr = tmp.c_str();

        sim_trial.save_vtk_name(cstr);
      }
    }
    std::cout << "v_norm_final " << v_norm << std::endl;

    f_calc->vd_calc_vel(m);

    sim_trial.save_vtk_name("./output/evolvedad");
    double time_curr = sim_trial.get_time();
    //sim_trial.set_time(ref_curr);
    int e_nbr = e_list->e.at(1).at(0).at(1).size();
    sim_trial.set_time(e_nbr);

    sim_trial.set_extract(true, ext_opts_file);
    sim_trial.extract_data();
    sim_trial.set_extract(false, ext_opts_file);

    sim_trial.set_time(time_curr);

    ss.clear();
    ss.str("");
    ss << FileName << ref_curr;
    temp = ss.str();
    tmp = temp + ".smb";

    const char* cstr = tmp.c_str();
    m->writeNative(cstr);

    tmp = temp + ".dmg";
    cstr = tmp.c_str();
    c_base->vd_write_dmg(cstr);

    ref_curr = ref_curr + mult;

    //sim_trial.set_adapt_param(ref_curr);
    sim_trial.set_adapt_bound_cells(1./ref_curr, cells_ref);

    for(int i = 0; i < ad_nbr; i++)
      sim_trial.adapt();
    split_interior_edges(m, f_calc, sim_trial);

    sim_trial.refresh_elist();
    sim_trial.save_vtk_name("./output/before");
    //v_th_tol = collect_relax(m, e_list, v_3c, et_3c, e_map);

    f_calc->upd_d2(m);
    f_calc->upd_gam2(m);
    f_calc->vd_att_fields(m);
    f_calc->vd_calc_vel(m);

    //v_norm = f_calc->vd_get_max_vel_norm_int(m);
  }

  delete e_con0;
  sim_trial.clean_up();

  PCU_Comm_Free();
  MPI_Finalize();
}

