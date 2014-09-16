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

// icc gridtest.cc -I./src/ -o grid -g -vec-report3 -vec -openmp -std=c++11

#include <grid.hh>
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>

int main(int argc, char* argv[])
{
  double xmin = 10.;
  double ymin = 10.;
  double zmin = 10.;
  int nx = 100;
  int ny = 100;
  int nz = 100;
  double xs = 2.;
  double ys = 2.;
  double zs = 2.;
  std::vector<double> x,y,z,v;
  std::vector<double> px,py,pz;

  for(int k = 0; k < nz; ++k) {
    for(int j = 0; j < ny; ++j) {
      for(int i = 0; i < nx; ++i) {
        for(int d = 0; d<100; ++d) {
          px.push_back(xmin+(i-d/100.)*xs);
          py.push_back(ymin+(j-d/100.)*ys);
          pz.push_back(zmin+(k-d/100.)*zs);
        }
        x.push_back(xmin+i*xs);
        y.push_back(ymin+j*ys);
        z.push_back(zmin+k*zs);
        v.push_back(i*i + j*j + k*k);
      }
    }
  }

  grid g;
  g.xmin = xmin;
  g.ymin = ymin;
  g.zmin = zmin;
  g.nx = nx;
  g.ny = ny;
  g.nz = nz;
//  g.xmax = xmin+(nx-1)*xs;
//  g.ymax = ymin+(ny-1)*ys;
//  g.zmax = zmin+(nz-1)*zs;
  g.xstride = xs;
  g.ystride = ys;
  g.zstride = zs;

//  for(int i = 0; i < x.size(); ++i) {
//    std::cout << "(" << x[i] << "," << y[i] << "," << z[i] << "): " << v[i] << std::endl;
//  }
int nstep = 10;
  double interped[px.size()];
  std::cout << "interp" << std::endl;
auto begin = std::chrono::high_resolution_clock::now();
  for(int i=0; i< nstep; ++i) g.interp3d(&px[0],&py[0],&pz[0],&v[0],px.size(),interped);
auto end = std::chrono::high_resolution_clock::now();
std::cout << std::setprecision(9)<<  std::chrono::duration_cast<std::chrono::nanoseconds>( end - begin).count()/1000000000./nstep << std::endl;
//  for(int i = 0; i < x.size(); ++i) std::cout << v[i] << " " << interped[i] << std::endl;
}
