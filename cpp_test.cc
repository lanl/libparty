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
#include <iostream>
#include <cstdio>
#include <ctime>
#include <math.h>
#include <iomanip>
#include <chrono>
#include <vector>

int main() {
  std::clock_t start;
  float duration;
//  start = std::clock();
//  double start = omp_get_wtime();
//  part_sys parts = part_sys(1000000);

  double t, tmax, dt;
  t = 0.0;
  tmax = 2.0;
  dt = 0.2;
  int nx = 1000;
  int ny = 1000;
  double xmin = -30;
  double ymin = -30;
  double xs = 2.*abs(xmin)/nx;
  double ys = 2.*abs(ymin)/ny;

  int nump = 100000;
  part_sys_soa parts = part_sys_soa(nump);
  parts.initgauss();
  parts.grid.xmin = xmin;
  parts.grid.ymin = ymin;
  parts.grid.zmin = 0;
  parts.grid.nx = nx;
  parts.grid.ny = ny;
  parts.grid.nz = 0;
  parts.grid.xstride = xs;
  parts.grid.ystride = ys;

  std::vector<double> vx, vy;
  for(int j = 0; j < ny; j++) {
    for(int i = 0; i < nx; i++) {
      double x = xmin + i * xs;
      double y = ymin + j * ys;
      double r = sqrt(x*x + y*y)+1.;
//      vx.push_back(-1./(1.)*y/abs(r));
//      vy.push_back(-1./(1.)*x/abs(r));
      vx.push_back(y/abs(r*r)+0.1);
      vy.push_back(-x/abs(r*r));
   
    }
  }

  parts.calc_density();
  parts.calc_density();
  parts.calc_density();

  bbox b;
  b.resize(6);
  b[0] = -10;
  b[1] = 10;
  b[2] = -10;
  b[3] = 10;
  b[4] = 0.;
  b[5] = 0.;

  for(; t < tmax; t+=dt) {
    parts.periodic2d(b);
    parts.interp2dvelwithgrid(&vx[0], &vy[0]); 
    parts.step(dt);
    parts.calc_density();
    parts.output_h5part();
    std::cout << "done step at time " << t << std::endl;
  }
/*
//  while(t<tmax)
for(int n = 1; n < 30; ++n)
  {
    int nump = pow(10.0,(double)n/5.0);
    part_sys parts = part_sys(nump);
//  start = std::clock();
    auto begin = std::chrono::high_resolution_clock::now();
    double bin = 10*0.5/nump;
    parts.build_hashtable(bin);
//    parts.find_neighbors(bin);
    parts.step(dt);
//  duration = ( std::clock() - start ) ;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout<<std::setprecision(9)<< nump << ","<< std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1000. << std::endl;
    t += dt;
  }
*/
//  parts.print_sys();

//  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
//  std::cout<<"Timer: "<< duration << std::endl;
                                 
}                                  
