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

#include "libparty.hh"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
//#include <boost/mpi.hpp>
//#include "libparty.h"
#include "MersenneTwister.hh"

inline double gaussian(double r, double sigma) {
  double val;
  double cons = 1.0/sigma/sqrt(2.0*3.1415);
  val = cons * exp(-(r*r)/(2*sigma*sigma));
  return val;
}
part_sys_soa::part_sys_soa( int npart )
{
  MTRand mtrand;
  double sphere_outerr = 10.0;
  for( int i = 0; i < npart; i++)
  {
    id.push_back(i);
//    double r, prob;
//    bool keep = false;
//    do {
//    x.push_back(sphere_outerr - mtrand.rand(2*sphere_outerr));
//    y.push_back(sphere_outerr - mtrand.rand(2*sphere_outerr));
//    z.push_back(sphere_outerr - mtrand.rand(2*sphere_outerr));
//      r = sqrt(p.pos[0]*p.pos[0] + p.pos[1]*p.pos[1] + p.pos[2]*p.pos[2]);
//      prob = gaussian(r,5.0);
//      if(mtrand.rand() < prob) keep = true;
//    } while (!keep);
//    } while(false);
//
//        p.mass = 1.;
//        p.nn_h = 0.5;
//        p.vel[0] = 0.2 - mtrand.rand(0.4);
//        p.vel[1] = 0.2 - mtrand.rand(0.4);
//        p.vel[2] = 0.2 - mtrand.rand(0.4);
//        p.del = false;
//        partv.push_back(p);
     }

  x.resize(npart,0.0);
  y.resize(npart,0.0);
  z.resize(npart,0.0);
  vx.resize(npart,0.0);
  vy.resize(npart,0.0);
  vz.resize(npart,0.0);
  mass.resize(npart,1.0);
  nn_h.resize(npart,1.);
  rho.resize(npart,0.0);
  ghost.resize(npart,0);

#ifndef USE_MIC
  //boost::mpi::environment env();
  if(boost::mpi::environment::initialized()){
    mproc = world.rank();
    nprocs = world.size();
  } else {
#endif
    mproc = 0;
    nprocs = 1;
#ifndef USE_MIC
  }
#endif
  nn_cutoff = pow(10.0,-6);
  max_nn_h = pow((1/(4*3.1415*nn_cutoff)),(1.0/3.0));
}

void part_sys_soa::initgauss()
{
  MTRand mtrand = MTRand(1337);
  double sphere_outerr = 10.0;
  for( int i = 0; i < x.size(); i++)
  {
    double r, prob;
    bool keep = false;
    do {
      x[i] = sphere_outerr - mtrand.rand(2*sphere_outerr);
      y[i] = sphere_outerr - mtrand.rand(2*sphere_outerr);
      z[i] = sphere_outerr - mtrand.rand(2*sphere_outerr);
      r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
      prob = gaussian(r,8.0);
      if(mtrand.rand() < prob) keep = true;
    } while (!keep);
 //    } while(false);
//    //
//    //        p.mass = 1.;
//    //        p.nn_h = 0.5;
    vx[i] = 0.2 - mtrand.rand(0.4);
    vy[i] = 0.2 - mtrand.rand(0.4);
    vz[i] = 0.2 - mtrand.rand(0.4);
//    //        p.del = false;
//    //        partv.push_back(p);
   }
}

void part_sys_soa::periodic2d(bbox b) {
  for(int i = 0; i < x.size(); ++i) {
    if(x[i] < b[0]) x[i] += b[1]-b[0];
    if(x[i] > b[1]) x[i] -= b[1]-b[0];
    if(y[i] < b[2]) y[i] += b[3]-b[2];
    if(y[i] > b[3]) y[i] -= b[3]-b[2];
  }
}
/*
void part_sys_soa::step( double dt )
{
  int numpart = x.size();
  double* cx = &x[0];
  double* cy = &y[0];
  double* cz = &z[0];
  double* cvx = &vx[0];
  double* cvy = &vy[0];
  double* cvz = &vz[0];
  short int* cg = &ghost[0];

  for( int proc = 0; proc < nprocs; ++proc) {
    if(mproc == proc && x.size() != 0) {
#pragma offload target(mic:0) in(cvx, cvy, cvz, cg : length(numpart) ) in(dt) inout(cx,cy,cz : length(numpart) alloc_if(1) free_if(1) )
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for( int i = 0; i < numpart; ++i) {
        if(cg[i]) continue;
        cx[i] += cvx[i] * dt;
        cy[i] += cvy[i] * dt;
        cz[i] += cvz[i] * dt;
        }
      }
  }
}
*/
void part_sys_soa::step( double dt )
{
  int numpart = x.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
for( int i = 0; i < numpart; ++i) {
if(ghost[i]) continue;
//#pragma omp parallel for
x[i] += vx[i] * dt;
//}
//#pragma omp parallel for
//for( int i = 0; i < numpart; ++i) {
////#pragma omp parallel for
y[i] += vy[i] * dt;
//}

//#pragma omp parallel for
//for( int i = 0; i < numpart; ++i) {
////#pragma omp parallel for
z[i] += vz[i] * dt;
}
}

std::ostream& operator<< (std::ostream &out, const part_sys_soa& p)
{
  for( int i = 0; i < p.x.size(); ++i) {
    out << "(" << p.x[i] << "," << p.y[i] << "," << p.z[i] << ")" << std::endl;
  }
  return out;
}

// Sync bboxes between all ranks
void part_sys_soa::syncbboxes(bbox lbox) {
  bboxes.resize(nprocs);
#ifndef USE_MIC
  boost::mpi::all_gather(world, lbox,bboxes);
#endif
  have_bboxes = true;
}

bool checkbbox_soa( double posx, double posy, double posz, bbox b) {
  return (posx > b[0]) && (posx < b[1]) && \
         (posy > b[2]) && (posy < b[3]) && \
         (posz > b[4]) && (posz < b[5]);
}

// return the rank that a particle p lives on (sync bboxes first)
int part_sys_soa::checkproc( int pind ) {
  if (checkbbox_soa(x[pind], y[pind], z[pind], bboxes[mproc])) {
    return mproc;
  } else {
    for(int proc=0; proc < nprocs; proc++) {
      if (checkbbox_soa(x[pind], y[pind], z[pind], bboxes[proc])) {
        return proc;
      }
    }
    return -1;
  }
}

void part_sys_soa::boostsend() {
#ifndef USE_MIC
  //Setup Send
  std::vector<std::vector<part> > sendpart;
  sendpart.resize(nprocs);

  for(int p=0; p < nprocs; p++) {
//    part temp;
    if( p == mproc ) {
      for( int i = 0; i < x.size(); ++i ) {
//        if((*it).proc != p) {
        if(ghost[i]) continue;
        int proc = checkproc(i);
        if(proc!=p) {
          if ( (proc >= 0) && (proc < nprocs)) {
            sendpart[proc].push_back(poppart(i));
          }
          //erasepart(i);
//          ghost[i] = true;
        }
      }
      //partv.erase(std::remove_if(partv.begin(), partv.end(), delpart), partv.end());
    }
  }

  // Setup recv
  std::vector<std::vector<part> > recvpart;
  recvpart.resize(nprocs);

  // let boost do its mpi magic
  boost::mpi::all_to_all(world, sendpart, recvpart);

  // Push the recved parts into partv
  int numrecv = 0;
  for(int q=0; q < nprocs; q++){
    if( recvpart[q].size() > 0){
      int recv_size = recvpart[q].size();
      for ( int i = 0; i < recv_size; ++i) {
        pushpart(recvpart[q][i]);
      }
    }
    numrecv += recvpart[q].size();
  }
  //if(verbose > 0)  
  //  std::cout << "recved: " << numrecv << std::endl;
  //}
#endif
} 

void part_sys_soa::sync_ghosts (double buff)
{
 clearghosts();
//  int gnum=0;
//  for(int i = ghost.size(); i > 0; --i) {
//    if(ghost[i]) {
//      erasepart(i);
//      gnum++;
//    }
//  }
//  for(std:vector<bool>::iterator it = ghost.begin(); it < ghost.end();) {
//  std:vector<bool>::iterator it = ghost.begin();
//  while(it != ghost.end()) {
//    if(*it) erasepart(i)
//    it++;
//  }
//  std::cout << "erased " << gnum << std::endl;
  //clearghosts();

  std::vector<double> buffsizes;

#ifndef USE_MIC
  boost::mpi::all_gather(world, buff, buffsizes);
#endif

  std::vector<bbox> buffboxes;
  buffboxes.resize(nprocs);
  std::vector<std::vector<part>> sendpart;
  sendpart.resize(nprocs);
  for(int proc=0; proc < nprocs; proc++)
  {
    if(proc == mproc) continue;
    buffboxes[proc].resize(6);
    for(int d=0; d<6; d=d+2) {
      buffboxes[proc][d] = bboxes[proc][d] - buffsizes[proc];
      buffboxes[proc][d+1] = bboxes[proc][d+1] + buffsizes[proc];
    }
    for( int p = 0; p < x.size(); ++p ) {
      if(checkbbox_soa(x[p], y[p], z[p], buffboxes[proc])) {
        sendpart[proc].push_back(poppart(p));
      }
    }
  }

  std::vector<std::vector<part> > recvpart;
  recvpart.resize(nprocs);

#ifndef USE_MIC
  boost::mpi::all_to_all(world, sendpart, recvpart);
#endif

  int numrecv = 0;
  for(int q=0; q < nprocs; q++){
    if( q == mproc) continue;
    if( recvpart[q].size() > 0){
      int recv_size = recvpart[q].size();
      for ( int i = 0; i < recv_size; ++i) {
        pushpart(recvpart[q][i],true);
      }
    }
    numrecv += recvpart[q].size();
  }
//  std::cout << "recved " << numrecv << std::endl;
}
 
