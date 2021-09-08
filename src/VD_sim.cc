// Adapted from ma_test.cc in SCOREC test directory.

#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>

#include <PCU.h>

#include <apfNumbering.h>
#include <apfShape.h>

int main(int argc, char** argv)
{

  assert(argc==3);
  //const char* modelFile = argv[1];
  //const char* meshFile = argv[2];

  const char* modelFile = "mshfiles/boxnbox.tess";
  const char* meshFile = "mshfiles/boxnbox.msh";

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();

  apf::Mesh2* m = apf::loadMdsFromGmsh(gmi_load(modelFile), meshFile);
  // Iterate:
/*
  apf::MeshIterator* it_e = m->begin(3);
  apf::MeshEntity* ent;
  while ((ent = m->iterate(it_e))) {
    apf::Downward down;
    m->getDownward(ent, 0, down);
  }
  m->end(it_e);
*/
  m->verify();

  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();
}

