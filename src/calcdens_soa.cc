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
#include "splineutils.h"

void part_sys_soa::calc_density() {

//  std::vector<std::vector<part*>> neighbors;
  tbb::concurrent_vector<tbb::concurrent_vector<int>> neighbors;
  std::vector<double> search_distv;

  for( int i = 0; i < x.size(); ++i ) {
      if(ghost[i]) continue;
      if(nn_h[i] < max_nn_h) search_distv.push_back(nn_h[i]);
      else search_distv.push_back(max_nn_h);
  }

  double min_nn_h, max_nn_h;
  if(search_distv.empty()) {
    min_nn_h = 0.0;
    max_nn_h = 0.0;
  } else {
    auto minmax = std::minmax_element(search_distv.begin(), search_distv.end());
    min_nn_h = *minmax.first;
    max_nn_h = *minmax.second;
  }

  if(boost::mpi::environment::initialized()) sync_ghosts(max_nn_h);

  build_hashtable(min_nn_h);

  neighbors = find_neighbors(32);

  int pi = 0;

//  for( partIter p = partv.begin(); p != partv.end(); ++p ) {
//  #pragma omp parallel for
  for(int i = 0; i < x.size(); ++i) {
    if(ghost[i]) continue;
    double distance, weight, Dhweight;
    double newN, dNdh, h_err;
    double eta = 1.2;
    double h_guess = nn_h[i];
    double h_new;
    bool keep_niter;
    int niter, maxiter;
    maxiter = 200;
    keep_niter = true;
    niter = 0;
    while(keep_niter) {
      newN = 0.0;
      dNdh = 0.0;

//      for(p_partIter n = neighbors[pi].begin(); n != neighbors[pi].end(); ++n) {
      for(int n = 0; n < neighbors[i].size(); ++n) {
        int ni;
        ni = (int)neighbors[i][n];
//        nv = neighbors[i*VECL+v][n]%VECL;

        double dist2 = pow(x[ni] - x[i],2.0)
                     + pow(y[ni] - y[i],2.0)
                     + pow(z[ni] - z[i],2.0);

        distance = sqrt(dist2);
        weight = cubic_spline(distance,h_guess);
        Dhweight = Dhcubic_spline(distance,h_guess);
        if (weight > nn_cutoff) {
          newN += mass[ni] * weight;
          dNdh += mass[ni] * Dhweight;
        }
      }

      h_new = eta/pow(newN,(1.0/3.0));
    //  std::cout << h_new << std::endl;
      h_err = h_new - h_guess;
      if(abs(h_err) < pow(10.0,-8.0)) {
        rho[i] = newN;
        nn_h[i] = h_new;
        keep_niter = false;
      } else {
        double dFdh = -1.0 - (eta/3.0)*dNdh/(pow(newN,(4.0/3.0)));
        h_guess = h_guess - h_err/dFdh;
      }
      if(niter > maxiter || !isfinite(h_new)) {
        keep_niter = false;
//        this -= ind;
//TODO DELTE!!!!        
//partvfill[i][v] = 0;
//        erasepart(i);;
//        del[i] = true;
//        ghost[i] = true;
        rho[i] = 0;
        std::cout << "failed" << " " << h_guess<< " "<< h_err << "n nn: " << neighbors[i].size() << std::endl;
      }
      niter++;
    }
    pi++;
  }
//std::cout << x[0] << " " << rho[0] << std::endl;
}
