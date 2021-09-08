#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_null.h>
#include <PCU.h>

// newdim.cc from SCORECtest

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, false);
  apf::Vector3 points[4] = {
    apf::Vector3(0,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,1,0),
    apf::Vector3(0,0,1)
  };
  apf::Vector3 param(0,0,0);
  apf::ModelEntity* surface = m->findModelEntity(2, 0);
  apf::ModelEntity* surface2 = m->findModelEntity(2, 1);
  apf::MeshEntity* verts[4];
  for (int i = 0; i < 4; ++i)
    verts[i] = m->createVertex(surface, points[i], param);
  for (int i = 0; i < 2; ++i) {
    apf::MeshEntity* tv[3];
    for (int j = 0; j < 3; ++j)
      tv[j] = verts[apf::tet_tri_verts[i][j]];
    apf::buildElement(m, surface, apf::Mesh::TRIANGLE, tv);
  }
  for (int i = 2; i < 4; ++i) {
    apf::MeshEntity* tv[3];
    for (int j = 0; j < 3; ++j)
      tv[j] = verts[apf::tet_tri_verts[i][j]];
    apf::buildElement(m, surface2, apf::Mesh::TRIANGLE, tv);
  }

  for (int i = 0; i < 3; ++i) {
    apf::MeshIterator* it = m->begin(i);
    apf::MeshEntity* ee;

    while ((ee = m->iterate(it))) {

      int e_type = m->getType(ee);
      int d = m->typeDimension[e_type];
      int type = m->getModelType(m->toModel(ee));
      int tag = m->getModelTag(m->toModel(ee));
      std::cout << d << "-ent " << ee
                << " " << type << "c" << tag << std::endl;
      apf::Downward dw;
      if(i > 0) {
        int d_count = m->getDownward(ee, i-1, dw);

        for(int d_i = 0; d_i < d_count; d_i++) {
          int e_type_d = m->getType(dw[d_i]);
          int d_d = m->typeDimension[e_type_d];
          int type_d = m->getModelType(m->toModel(dw[d_i]));
          int tag_d = m->getModelTag(m->toModel(dw[d_i]));
          std::cout << "\t" << d_d << "-ent " << dw[d_i]
                    << " " << type_d << "c" << tag_d << std::endl;
        }
      }
    }
    m->end(it);
  }

  m->acceptChanges();
  m->verify();
  apf::changeMdsDimension(m, 3);
  apf::ModelEntity* interior = m->findModelEntity(3, 0);
  apf::buildElement(m, interior, apf::Mesh::TET, verts);
  m->acceptChanges();
  m->verify();
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


