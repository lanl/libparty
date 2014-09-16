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

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
//#pragma offload_attribute(push target(mic))
#include "libparty.hh"
#include <iostream>
#include <omp.h>
#include "MersenneTwister.hh"
#include "splineutils.h"
//#include "offload.h"
//#pragma offload_attribute(pop)
inline double gaussian(double r, double sigma) {
  double val;
  double cons = 1.0/sigma/sqrt(2.0*3.1415);
  val = cons * exp(-(r*r)/(2*sigma*sigma));
  return val;
}

void part_sys_aosoa::step( double dt )
{
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for( int i = 0; i < partv.size(); ++i) {
//    double x[16],y[16],z[16],d[16];
//#pragma omp parallel for
    partv[i].x[:] += partv[i].vx[:] * dt;
//#pragma omp parallel for
    partv[i].y[:] += partv[i].vy[:] * dt;
//#pragma omp parallel for
    partv[i].z[:] += partv[i].vz[:] * dt;
  }
}
//#pragma omp parallel for
//y[i] += vy[i] * dt;
//#pragma omp parallel for
//z[i] += vz[i] * dt;
//}
/*  int numpart = x.size();
  double* cx = &x[0];
  double* cy = &y[0];
  double* cz = &z[0];
  double* cvx = &vx[0];
  double* cvy = &vy[0];
  double* cvz = &vz[0];
if(mic_alloced) {
//  #pragma offload target(mic) in(cvx : length(numpart) ) in(dt) out(cx : length(numpart) alloc_if(0) free_if(0) )
#ifdef _OPENMP
  #pragma omp parallel for 
#endif
  for( int i = 0; i < numpart; ++i) cx[i] += cvx[i] * dt;

//  #pragma offload target(mic) in(cvy : length(numpart) ) in(dt) out(cy : length(numpart) alloc_if(0) free_if(0) )
#ifdef _OPENMP
  #pragma omp parallel for 
#endif
  for( int i = 0; i < numpart; ++i) cy[i] += cvy[i] * dt;

//  #pragma offload target(mic) in(cvz : length(numpart) ) in(dt) out(cz : length(numpart) alloc_if(0) free_if(0) )  
#ifdef _OPENMP
  #pragma omp parallel for 
#endif
  for( int i = 0; i < numpart; ++i) cz[i] += cvz[i] * dt;

}else{

//#pragma offload target(mic) in(cvx : length(numpart) ) in(dt) inout(cx : length(numpart) alloc_if(1) free_if(0) )
#ifdef _OPENMP
  #pragma omp parallel for 
#endif
  for( int i = 0; i < numpart; ++i) cx[i] += cvx[i] * dt;

//#pragma offload target(mic) in(cvy : length(numpart) ) in(dt) inout(cy : length(numpart) alloc_if(1) free_if(0) )
#ifdef _OPENMP
  #pragma omp parallel for 
#endif
  for( int i = 0; i < numpart; ++i) cy[i] += cvy[i] * dt;

//#pragma offload target(mic) in(cvz : length(numpart) ) in(dt) inout(cz : length(numpart) alloc_if(1) free_if(0) )  
#ifdef _OPENMP
  #pragma omp parallel for 
#endif
  for( int i = 0; i < numpart; ++i) cz[i] += cvz[i] * dt;

mic_alloced = true;
}*/
  //(&z[0])[0:size:] += (&vz[0])[:] * dt;
//    y[i] += vy[i] * dt;
//    z[i] += vz[i] * dt;
//  }
//}


part_sys::part_sys( int npart )
{
  MTRand mtrand;
  double sphere_outerr = 10.0;
  for( int i = 0; i < npart; i++)
  {
    part p;
    p.id = i;
    double r, prob;
    bool keep = false;
    do {
      p.pos[0] = sphere_outerr - mtrand.rand(2*sphere_outerr);
      p.pos[1] = sphere_outerr - mtrand.rand(2*sphere_outerr);
      p.pos[2] = sphere_outerr - mtrand.rand(2*sphere_outerr);
      r = sqrt(p.pos[0]*p.pos[0] + p.pos[1]*p.pos[1] + p.pos[2]*p.pos[2]);
      prob = gaussian(r,5.0);
      if(mtrand.rand() < prob) keep = true;
    } while (!keep);
//    } while(false);

    p.mass = 1.;
    p.nn_h = 0.5;
    p.vel[0] = 0.2 - mtrand.rand(0.4);
    p.vel[1] = 0.2 - mtrand.rand(0.4);
    p.vel[2] = 0.2 - mtrand.rand(0.4);
    p.del = false;
    partv.push_back(p);
  }
//  l_npart = partv.size();
  mproc = world.rank();
  nprocs = world.size();
  nn_cutoff = pow(10.0,-6);
  //TODO better PI value
  max_nn_h = pow((1/(4*3.1415*nn_cutoff)),(1.0/3.0));
  
}

void part_sys::setids() {

  int localsize = partv.size();
  std::vector<int> numpart_per_proc;

#ifndef USE_MIC
  boost::mpi::all_gather(world, localsize, numpart_per_proc);
#endif
  int accum_part = 0;
  for (int proc = 0; proc < world.size(); ++proc) {
    if(proc==mproc) {
      for( partIter it = partv.begin(); it != partv.end(); ++it ) {
        (*it).id = accum_part + std::distance(partv.begin(),it);
      }
    } else {
      accum_part += numpart_per_proc[proc];
    }
  }
}

int part_sys::get_npart() { return partv.size(); }

void part_sys::step( double dt )
{
  int pimax = partv.size();
#ifdef _OPENMP
  #pragma omp parallel for collapse(2) 
#endif
  for( int i = 0; i < pimax; ++i) {
    for(int d = 0; d<3; ++d) partv[i].pos[d] += partv[i].vel[d] * dt;
  }
//  for( int i = 0; i < pimax; ++i) {
//partv[i].pos[1] += partv[i].vel[1] * dt;
//}
//  for( int i = 0; i < pimax; ++i) {
///partv[i].pos[2] += partv[i].vel[2] * dt;
//}
//  iter++;
}

void part_sys::print_sys( )
{
  //std::cout << "data @ " << &partv[0] << std::endl;
  for( partIter p = partv.begin(); p != partv.end(); ++p ) {
    std::cout << "Particle[" << p->id << "] @ (" << p->pos[0] << "," << p->pos[1] << "," << p->pos[2]
         << ") w/ vel=(" << p->vel[0] << "," << p->vel[1] << "," << p->vel[2] << ")" << " on rank " << mproc <<std::endl;
  }
}

bool delpart( const part p) {
  return p.del;
  //return false;
}

void part_sys::boostsend() {
#ifndef USE_MIC
  // Setup Send
  std::vector<std::vector<part> > sendpart;
  sendpart.resize(nprocs);

  for(int p=0; p < nprocs; p++) {
    part temp;
    if( p == mproc ) {
      for( partIter it = partv.begin(); it != partv.end(); ++it ) {
        if((*it).proc != p) {
          if ( ((*it).proc >= 0) && ((*it).proc < nprocs)) {
            sendpart[(*it).proc].push_back((*it));
          }
          (*it).del = true;
        }
      }
      partv.erase(std::remove_if(partv.begin(), partv.end(), delpart), partv.end());
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
        partv.push_back(recvpart[q][i]);
      }
    }
    numrecv += recvpart[q].size();
  }
  //if(verbose > 0)  
  //  std::cout << "recved: " << numrecv << std::endl;
  //}
#endif
}

// Check if a pos[3] is contained in a bbox
bool checkbbox( double x[3], bbox b) {
  return (x[0] > b[0]) && (x[0] < b[1]) && (x[1] > b[2]) && (x[1] < b[3]) && (x[2] > b[4]) && (x[2] < b[5]);
}

void part_sys::sync_ghosts (double buff)
{
  ghostv.clear();

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
    for( partIter p = partv.begin(); p != partv.end(); ++p ) {
      if(checkbbox(p->pos, buffboxes[proc])) {
        sendpart[proc].push_back(*p);
      }
    }
  }

  // Setup recv
  std::vector<std::vector<part> > recvpart;
  recvpart.resize(nprocs);

#ifndef USE_MIC
  // let boost do its mpi magic
  boost::mpi::all_to_all(world, sendpart, recvpart);
#endif

  // Push the recved parts into partv
  int numrecv = 0;
  for(int q=0; q < nprocs; q++){
    if( recvpart[q].size() > 0){
      int recv_size = recvpart[q].size();
      for ( int i = 0; i < recv_size; ++i) {
        ghostv.push_back(recvpart[q][i]);
      }
    }
    numrecv += recvpart[q].size();
  }
}

// Sync bboxes between all ranks
void part_sys::syncbboxes(bbox lbox) {
  bboxes.resize(nprocs);
#ifndef USE_MIC
  boost::mpi::all_gather(world, lbox,bboxes);
#endif
}

// return the rank that a particle p lives on (sync bboxes first)
int part_sys::checkproc( part p) {
  if (checkbbox(p.pos, bboxes[mproc])) {
    return mproc;
  } else {
    for(int proc=0; proc < nprocs; proc++) {
      if (checkbbox(p.pos, bboxes[proc])) {
        return proc;
      } 
    }
    //when we cant find a bbox, stage the particle to be deleted
//    p.del = true;
    return -1;
  }
}

// check and set the proc each local particle should be on
void part_sys::findprocs() {
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for( partIter p = partv.begin(); p != partv.end(); ++p ) {
    p->proc = checkproc(*p);
  }  
}

//TODO change var names from N to rho
void part_sys::calc_density() {

  tbb::concurrent_vector<tbb::concurrent_vector<int>> neighbors;
  std::vector<double> search_distv;

  for( partIter p = partv.begin(); p != partv.end(); ++p ) {
    if(p->nn_h < max_nn_h) search_distv.push_back(p->nn_h);
    else search_distv.push_back(max_nn_h);
  }

  double min_nn_h, max_nn_h;
  if(partv.empty()) {
    min_nn_h = 0.0;
    max_nn_h = 0.0;
  } else {
    auto minmax = std::minmax_element(search_distv.begin(), search_distv.end());
    min_nn_h = *minmax.first;
    max_nn_h = *minmax.second;
  }
//  double global_ghost_size;
//  boost::mpi::all_reduce(world, *minmax.second, global_ghost_size, boost::mpi::maximum<double>());

//  sync_ghosts(global_ghost_size);
  sync_ghosts(max_nn_h);
//ghostv.clear();
  build_hashtable(min_nn_h);

  neighbors = find_neighbors(search_distv,32);

  int pi = 0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for( partIter p = partv.begin(); p != partv.end(); ++p ) {
    double distance, weight, Dhweight;
    double newN, dNdh, h_err;
    double eta = 1.2;
    double h_guess = p->nn_h;
    double h_new;
    bool keep_niter;
    int niter, maxiter;
//    double SORfac = 0.8;
    maxiter = 200;
    keep_niter = true;
    niter = 0;
    while(keep_niter) {
      newN = 0.0;
      dNdh = 0.0;

//      for(p_partIter n = neighbors[pi].begin(); n != neighbors[pi].end(); ++n) {
      for(int ni = 0; ni < neighbors[pi].size(); ++ni) {
        int n = neighbors[pi][ni];
        double dist2 = pow(partv[n].pos[0] - p->pos[0],2.0)
                     + pow(partv[n].pos[1] - p->pos[1],2.0)
                     + pow(partv[n].pos[2] - p->pos[2],2.0);
        distance = sqrt(dist2);
        weight = cubic_spline(distance,h_guess);
        Dhweight = Dhcubic_spline(distance,h_guess);
        if (weight > nn_cutoff) {
          newN += partv[n].mass * weight;
          dNdh += partv[n].mass * Dhweight;
        }
      }

      h_new = eta/pow(newN,(1.0/3.0));
//      if(h_new > max_nn_h) h_new = max_nn_h;
      h_err = h_new - h_guess;
//      double dFdh = -1.0 - (eta/3.0)*dNdh/(pow(newN,(4.0/3.0)));
//      double new_h_guess = (1.0 - SORfac)*h_guess + SORfac*(h_guess - h_err/dFdh);
      
      if(abs(h_err) < pow(10.0,-8.0)) {
        p->rho = newN;
        p->nn_h = h_new;
        keep_niter = false;
//        g->Omega = 1.0 - (eta/3.0)*dNdh/(pow(newN,(4.0/3.0)));
//        ^^ thats one of the 'gradh' terms for sph
      } else {
//        h_guess = new_h_guess;
//        double SORfac = 0.8;
        double dFdh = -1.0 - (eta/3.0)*dNdh/(pow(newN,(4.0/3.0)));
//        h_guess = (1.0 - SORfac)*h_guess + SORfac*(h_guess - h_err/dFdh);
        h_guess = (h_guess - h_err/dFdh);
      }
      if(niter > maxiter || !isfinite(h_new)) {
        keep_niter = false;
        //partv.erase(p);
        p->del = true;
        std::cout << "failed" << " " << h_guess<< " "<< h_err << " n nn: " << neighbors[pi].size() << std::endl;
      }
      niter++;
    }
    pi++;
  }
std::cout << partv[0].pos[0] << " " << partv[0].rho << std::endl;
}
