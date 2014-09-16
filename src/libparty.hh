/*
 * Copyright 2014.  Los Alamos National Security, LLC. This material was produced 
 * under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
 * Laboratory (LANL), which is operated by Los Alamos National Security, LLC 
 * for the U.S. Department of Energy. The U.S. Government has rights to use, 
 * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR 
 * LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, 
 * OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is 
 * modified to produce derivative works, such modified software should be 
 * clearly marked, so as not to confuse it with the version available from LANL.   
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not 
 * use this file except in compliance with the License. You may obtain a copy 
 * of the License at http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software distributed 
 * under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the 
 * specific language governing permissions and limitations under the License.
 * 
 *  Additionally, redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Los Alamos National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *
 * LibParty, Version 1.x, Copyright number C114122, LA-CC-14-086
 *
 * Authors:
 *    Matthew Kinsey, kinsey@lanl.gov
 *    Verinder Rana, vrana@lanl.gov
 *    Bob Robey, XCP-2, LANL, brobey@lanl.gov
 *    Jeremy Sauer, EES-16, LANL, jsauer@lanl.gov
 *    Jon Reisner, XCP-4 LANL, reisner@lanl.gov
 */ 

#ifndef VECL
  #define VECL 8
#endif
#ifndef ALIGNL 
  #define ALIGNL 32
#endif
#define alignled __declspec(align(ALIGNL))

#ifndef PARTICLE_H
#define PARTICLE_H
//#pragma offload_attribute(push target(mic))
#include <vector>
//#pragma offload_attribute(pop)
#include <unordered_map>
#include <cstdint>
//#include "util.hh"
#include <iostream>
#include "libparty.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/mpi.hpp>
#include <iostream>
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#include <bitset>
#include "grid.hh"

typedef std::vector<double> bbox;
typedef std::vector<bbox> bboxv;

struct part : CCpart
{
private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
// NOTE!: YOU HAVE TO SERIALIZE _EVERY_ VARIABLE IN THE CLASS OR IT BREAKS
    ar & pos & vel & mass & id & del & proc & nn_h & rho;
  }
};

BOOST_IS_MPI_DATATYPE(part)

struct partvec {
  __declspec(align(ALIGNL)) double x[VECL], y[VECL], z[VECL];
  __declspec(align(ALIGNL)) double vx[VECL], vy[VECL], vz[VECL];
  __declspec(align(ALIGNL)) double mass[VECL], nn_h[VECL], rho[VECL];
//  int id[VECL], proc[VECL], del[VECL];
};
class part_sys_aosoa
{
  public:
  boost::mpi::environment env();
  boost::mpi::communicator world;

  int mproc;// = world.rank();
  int nprocs;// = world.size();


  std::vector<partvec> partv;
  std::vector<std::bitset<VECL>> partvfill;
  void step( double dt );
  void build_hashtable(double bin_size);
  double current_bin_size;
  int table_size;
  tbb::concurrent_unordered_multimap<unsigned int,int> hash_table;
  tbb::concurrent_vector<tbb::concurrent_vector<int>> find_neighbors(int nmax);
  double nn_cutoff = pow(10.0,-6);
  double max_nn_h = pow((1/(4*3.1415*nn_cutoff)),(1.0/3.0));

  //IO
  int iter = 0;
  int openslot = 0;
  void output_h5part();

  void calc_density();
 
  void operator+=(const part &p) {
    if(openslot%VECL == 0) {
      partvec pv;
      pv.x[0] = p.pos[0];
      pv.y[0] = p.pos[1];
      pv.z[0] = p.pos[2];
      pv.vx[0] = p.vel[0];
      pv.vy[0] = p.vel[1];
      pv.vz[0] = p.vel[2];
      pv.mass[0] = p.mass;
      pv.nn_h[0] = p.nn_h;
      std::bitset<VECL> fill;
      fill[0]=1;
      partv.push_back(pv);
      partvfill.push_back(fill);
//      openslot++;
    } else {
//      for(int i = partv.size() - 1; i >= 0; --i) {
//        for(int b = 0; b < VECL; ++b) {
//          if(!partvfill[i][b]) {
int oi = (int)(openslot/VECL);
int ob = (int)(openslot%VECL);
            partv[oi].x[ob] = p.pos[0];
            partv[oi].y[ob] = p.pos[1];
            partv[oi].z[ob] = p.pos[2];
            partv[oi].vx[ob] = p.vel[0];
            partv[oi].vy[ob] = p.vel[1];
            partv[oi].vz[ob] = p.vel[2];
            partv[oi].mass[ob] = p.mass;
            partv[oi].nn_h[ob] = p.nn_h;
            partvfill[oi][ob] = 1;
//            return;
//          }
//        }
      }
openslot++;
}
      //we didnt find an empty spot :(
      //      partvec pv;
      //            pv.x[0] = p.pos[0];
      //                  pv.y[0] = p.pos[1];
      //                        pv.z[0] = p.pos[2];
      //                              pv.vx[0] = p.vel[0];
      //                                    pv.vy[0] = p.vel[1];
      //                                          pv.vz[0] = p.vel[2];
      //                                                pv.mass[0] = p.mass;
      //                                                      pv.nn_h[0] = p.nn_h;
      //                                                            std::bitset<VECL> fill;
      //                                                                  fill[0]=1;
      //                                                                        partv.push_back(pv);
      //                                                                              partvfill.push_back(fill);
      //                                                                                  }
      //                                                                                    }
/*
  void operator+=(const part &p) {
    if(partv.empty()) {
      partvec pv;
      pv.x[0] = p.pos[0];
      pv.y[0] = p.pos[1];
      pv.z[0] = p.pos[2];
      pv.vx[0] = p.vel[0];
      pv.vy[0] = p.vel[1];
      pv.vz[0] = p.vel[2];
      pv.mass[0] = p.mass;
      pv.nn_h[0] = p.nn_h;
      std::bitset<VECL> fill;
      fill[0]=1;
      partv.push_back(pv);
      partvfill.push_back(fill);
    } else {
      for(int i = partv.size() - 1; i >= 0; --i) {
        for(int b = 0; b < VECL; ++b) {
          if(!partvfill[i][b]) {
            partv[i].x[b] = p.pos[0];
            partv[i].y[b] = p.pos[1];
            partv[i].z[b] = p.pos[2];
            partv[i].vx[b] = p.vel[0];
            partv[i].vy[b] = p.vel[1];
            partv[i].vz[b] = p.vel[2];
            partv[i].mass[b] = p.mass;
            partv[i].nn_h[b] = p.nn_h;
            partvfill[i][b] = 1;
            return;
          }
        }
      }
      //we didnt find an empty spot :(
      partvec pv;
      pv.x[0] = p.pos[0];
      pv.y[0] = p.pos[1];
      pv.z[0] = p.pos[2];
      pv.vx[0] = p.vel[0];
      pv.vy[0] = p.vel[1];
      pv.vz[0] = p.vel[2];
      pv.mass[0] = p.mass;
      pv.nn_h[0] = p.nn_h;
      std::bitset<VECL> fill;
      fill[0]=1;
      partv.push_back(pv);
      partvfill.push_back(fill);
    }
  }
*/
  void operator-=(const int &ind) {
    int i, v;
    i = ind/VECL;
    v = ind%VECL;
    partvfill[i][v] = 0;
    if(partvfill[i].none()) {
      partvfill.erase(partvfill.begin()+i);
      partv.erase(partv.begin()+i);
    }
  }

};


class part_sys_soa
{
  public:
//TODO add pointers to these to an array for iteration
  __declspec(align(ALIGNL)) std::vector<double> x, y, z, vx, vy, vz, mass, nn_h, rho;
  __declspec(align(ALIGNL)) std::vector<double> gx, gy, gz, gvx, gvy, gvz, gmass, gnn_h, grho;
  __declspec(align(ALIGNL)) std::vector<int> id, proc;
  __declspec(align(ALIGNL)) std::vector<bool> ghost;
//  __declspec(align(ALIGNL)) std::vector<int> del;

  double * vxptr, * vyptr, * vzptr;

  bool mic_alloced = false;
  bool have_grid = false;
  bool have_bboxes = false;

  grid grid;
  void interp2dvelwithgrid(double * gvx, double * gvy, bool fort_order = false) {
    grid.interp2d(&x[0], &y[0], gvx, x.size(), &vx[0], fort_order);
    grid.interp2d(&x[0], &y[0], gvy, x.size(), &vy[0], fort_order);
  }

  void interp3dvelwithgrid(double * gvx, double * gvy, double * gvz, bool fort_order = false) {
    grid.interp3d(&x[0], &y[0], &z[0], gvx, x.size(), &vx[0], fort_order);
    grid.interp3d(&x[0], &y[0], &z[0], gvy, x.size(), &vy[0], fort_order);
    grid.interp3d(&x[0], &y[0], &z[0], gvz, x.size(), &vz[0], fort_order);
  }
  void periodic2d(bbox b);
  void initgauss();

///bool is_ghost( auto itb, auto it) {
//  int i = std::distance(itb, it)
//  return ghost[i];
  //return false;
//}
  boost::mpi::environment env();
  boost::mpi::communicator world;

  int mproc;// = world.rank();
  int nprocs;// = world.size();

  part_sys_soa( int n_part );

//  part_data operator+=(const part &p);
  void clearghosts() {
  //  x.erase(std::remove_if(partv.begin(), partv.end(), delpart), partv.end());
    std::vector<double>::iterator x_it = x.begin();
    std::vector<double>::iterator y_it = y.begin();
    std::vector<double>::iterator z_it = z.begin();
    std::vector<double>::iterator vx_it = vx.begin();
    std::vector<double>::iterator vy_it = vy.begin();
    std::vector<double>::iterator vz_it = vz.begin();
    std::vector<double>::iterator mass_it = mass.begin();
    std::vector<double>::iterator nn_h_it = nn_h.begin();
    std::vector<double>::iterator rho_it = rho.begin();
    std::vector<bool>::iterator ghost_it = ghost.begin();
    int delnum = 0;
    //std::cout << "blah " << x.size() << " " << ghost.size() << std::endl;
    while(ghost_it != ghost.end()) {
      if(*ghost_it) {
        x_it = x.erase(x_it);
        y_it = y.erase(y_it);
        z_it = z.erase(z_it);
        vx_it = vx.erase(vx_it);
        vy_it = vy.erase(vy_it);
        vz_it = vz.erase(vz_it);
        mass_it = mass.erase(mass_it);
        nn_h_it = nn_h.erase(nn_h_it);
        rho_it = rho.erase(rho_it);
        ghost_it = ghost.erase(ghost_it);
        delnum++;
      } else {
        ++x_it;
        ++y_it;
        ++z_it;
        ++vx_it;
        ++vy_it;
        ++vz_it;
        ++mass_it;
        ++nn_h_it;
        ++rho_it;
        ++ghost_it;
      }
    }
    //std::cout << "removed " << delnum << std::endl;
  }
  void pushpart(const part &p, bool ghostp = false) {
    x.push_back(p.pos[0]);
    y.push_back(p.pos[1]);
    z.push_back(p.pos[2]);
    vx.push_back(p.vel[0]);
    vy.push_back(p.vel[1]);
    vz.push_back(p.vel[2]);
    mass.push_back(p.mass);
    nn_h.push_back(p.nn_h);
    rho.push_back(p.rho);
//    del.push_back(p.del);
    ghost.push_back(ghostp);
/*    if(ghostp) {
    gx.push_back(p.pos[0]);
    gy.push_back(p.pos[1]);
    gz.push_back(p.pos[2]);
    gvx.push_back(p.vel[0]);
    gvy.push_back(p.vel[1]);
    gvz.push_back(p.vel[2]);
    gmass.push_back(p.mass);
    gnn_h.push_back(p.nn_h);
    grho.push_back(p.rho);
    }*/
  }
  void operator+=(const part &p) {
    pushpart(p);
  }

  void erasepart(const int i) {
    x.erase(x.begin()+i);
    y.erase(y.begin()+i);
    z.erase(z.begin()+i);
    vx.erase(vx.begin()+i);
    vy.erase(vy.begin()+i);
    vz.erase(vz.begin()+i);
    mass.erase(mass.begin()+i);
    nn_h.erase(nn_h.begin()+i);
    rho.erase(rho.begin()+i);
//    del.erase(del.begin()+i);
    ghost.erase(ghost.begin()+i);
  }
  void operator-=(const int i) {
    erasepart(i);
  }

  part poppart(const int i) {
    part p;
    p.pos[0] = x[i];
    p.pos[1] = y[i];
    p.pos[2] = z[i];
    p.vel[0] = vx[i];
    p.vel[1] = vy[i];
    p.vel[2] = vz[i];
    p.mass = mass[i];
    p.nn_h = nn_h[i];
    p.rho = rho[i];
//    p.del = del[i];
    return p;
  }
  part operator[](const int i) {
    return poppart(i);
  }
  void step( double dt );
  void build_hashtable(double bin_size);
  tbb::concurrent_vector<tbb::concurrent_vector<int>> find_neighbors(int nmax);
  double current_bin_size;
  int table_size;
  tbb::concurrent_unordered_multimap<unsigned int,int> hash_table;
  double nn_cutoff;// = pow(10.0,-6);
  double max_nn_h;// = pow((1/(4*3.1415*nn_cutoff)),(1.0/3.0));
  void calc_density(); 

  //MPI
  bboxv bboxes;
  void boostsend();
  void syncbboxes(bbox lbbox);
  int checkproc( int p);
//  void findprocs();
//  void setids();
  void sync_ghosts( double buffersize );

  //IO
  int iter = 0;
  void output_h5part();
//  std::ostream& operator<< (std::ostream &out, const part_data& p)
//  {
//    for( int i = 0; i < x.size(); ++i) {
//      out << "(" << p.x[i] << "," << p.y[i] << "," << p.z[i] << std::endl;
//   }
//    return out;
//  }
  friend std::ostream & operator<<(std::ostream &, const part_sys_soa &);
};


//typedef std::vector<double> bbox;
//typedef std::vector<bbox> bboxv;

typedef std::vector<part>::iterator partIter;
typedef std::vector<part>::const_iterator c_partIter;
//typedef std::vector<part*>::iterator p_partIter;
typedef tbb::concurrent_vector<int>::iterator p_partIter;

typedef std::pair<unsigned int,int> keypart_p;
typedef std::pair<unsigned int,uintptr_t> ppart_p;
typedef std::pair<unsigned int,std::pair<unsigned int,uintptr_t>> keyprocpart_p;


class part_sys
{
  public:
  boost::mpi::environment env();
  boost::mpi::communicator world;

  int mproc;// = world.rank();
  int nprocs;// = world.size();

  std::vector<part> partv;
  std::vector<part> ghostv;
    tbb::concurrent_unordered_multimap<unsigned int,int> hash_table;
  double current_bin_size;
  int table_size;
  bboxv bboxes;
//  int l_npart;

  //Basic Stuff
  part_sys( int n_part );
  void step( double dt );
  int get_npart();
  void print_sys();

  //NN stuff
  double min_bin_h = 1.0;
  double nn_cutoff;// = pow(10.0,-6);
  //TODO better PI value
  double max_nn_h;// = pow((1/(4*3.1415*nn_cutoff)),(1.0/3.0));
  void build_hashtable(double bin_size);
//  std::vector<std::vector<part*>> find_neighbors(/*std::vector<part> parts,*/ double search_dist);
//  std::vector<std::vector<part*>> find_neighbors(std::vector<double> search_distv, int nmax);
  tbb::concurrent_vector<tbb::concurrent_vector<int>> find_neighbors(std::vector<double> search_distv, int nmax);
  //MPI
  void boostsend();
  void syncbboxes(bbox lbbox);
  int checkproc( part p);
  void findprocs();
  void setids();
  void sync_ghosts( double buffersize );

  //IO
  int iter = 0;
  void output_h5part();

  void calc_density();
};
#endif
