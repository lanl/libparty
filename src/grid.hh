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

#ifndef GRID_HH
#define GRID_HH

#include <iostream>

class grid {
  public:
  double xmin, ymin, zmin;
//  double xmax, ymax, zmax;
  int nx, ny, nz;
  double xstride, ystride, zstride;

//  void interp(double * x, double * y, double * z, double * val, int size, double * ret);
//  int ind1d(int i, int j, int k);
//};

int ind1d(int i, int j, int k, bool fort_order = false) {
//  int nx, ny;
//  nx = (int)(( xmax - xmin ) / xstride) + 1;
//  ny = (int)(( ymax - ymin ) / ystride) + 1;
//  nz = (int)(( zmax - zmin ) / zstride);
//  if(fort_order) return k + nz*j + nz*ny*i;
  return i + nx*j + nx*ny*k;
}

void interp3d(double * x, double * y, double * z, double * val, int size, double * ret, bool fort_order = false) {

//  int i[size], j[size], k[size];
//  double x0[size], y0[size], z0[size];
//  double xd[size], yd[size], zd[size];
//  double c00[size], c10[size], c01[size], c11[size], c0[size], c1[size];
//  for(int i = 0; i < nx; ++i) std::cout << va

  #pragma ivdep
  #pragma omp parallel for
  for(int n = 0; n < size; ++n) {
    int i = (int)(( x[n] - xmin ) / xstride);
    int j = (int)(( y[n] - ymin ) / ystride);
    int k = (int)(( z[n] - zmin ) / zstride);

    double x0 = (i) * xstride + xmin;
    double y0 = (j) * ystride + ymin;
    double z0 = (k) * zstride + zmin;

    double xd = ( x[n] - x0 ) / xstride;
    double yd = ( y[n] - y0 ) / ystride;
    double zd = ( z[n] - z0 ) / zstride;

    double c00 = val[ind1d(i,j,k,fort_order)]       * ( 1. - xd )
               + val[ind1d(i+1,j,k,fort_order)]     * xd;
    double c10 = val[ind1d(i,j+1,k,fort_order)]     * ( 1. - xd )
               + val[ind1d(i+1,j+1,k,fort_order)]   * xd;
    double c01 = val[ind1d(i,j,k+1,fort_order)]     * ( 1. - xd )
               + val[ind1d(i+1,j,k+1,fort_order)]   * xd;
    double c11 = val[ind1d(i,j+1,k+1,fort_order)]   * ( 1. - xd )
               + val[ind1d(i+1,j+1,k+1,fort_order)] * xd;

    double c0 = c00 * ( 1. - yd ) + c10 * yd;
    double c1 = c01 * ( 1. - yd ) + c11 * yd;

    ret[n] = c0 * ( 1. - zd ) + c1 * zd;
  }
}

void interp2d(double * x, double * y, double * val, int size, double * ret, bool fort_order = false) {
//  int i[size], j[size];
//  double x0[size], y0[size];
//  double xd[size], yd[size];
//  double c0[size], c1[size];
  #pragma ivdep
  #pragma omp parallel for
  for(int n = 0; n < size; ++n) {
    int i = (int)(( x[n] - xmin ) / xstride);
    int j = (int)(( y[n] - ymin ) / ystride);

    double x0 = i * xstride + xmin;
    double y0 = j * ystride + ymin;

    double xd = ( x[n] - x0 ) / xstride;
    double yd = ( y[n] - y0 ) / ystride;

    double c0 = val[ind1d(i,j,0,fort_order)]     * (1. - xd) 
              + val[ind1d(i+1,j,0,fort_order)]   * xd;
    double c1 = val[ind1d(i,j+1,0,fort_order)]   * (1. - xd) 
              + val[ind1d(i+1,j+1,0,fort_order)] * xd;

    ret[n] = c0 * (1. - yd) + c1 * yd;
  }
//  std::cout << ret[0] << std::endl;
}
};
/*
    i[:] = (int)(( x[0:size] - xmin ) / xstride);
    j[:] = (int)(( y[0:size] - ymin ) / ystride);
    k[:] = (int)(( z[0:size] - zmin ) / zstride);

    x0[:] = i[:] * xstride + xmin;
    y0[:] = j[:] * ystride + ymin;
    z0[:] = k[:] * zstride + zmin;

    xd[:] = ( x[0:size] - x0[:] ) / xstride;
    yd[:] = ( y[0:size] - y0[:] ) / ystride;
    zd[:] = ( z[0:size] - z0[:] ) / zstride;

    c00[:] = val[ind1d(i[:],j[:],k[:])]       * ( 1. - xd[:] )
           + val[ind1d(i[:]+1,j[:],k[:])]     * xd[:];
    c10[:] = val[ind1d(i[:],j[:]+1,k[:])]     * ( 1. - xd[:] )
           + val[ind1d(i[:]+1,j[:]+1,k[:])]   * xd[:];
    c01[:] = val[ind1d(i[:],j[:],k[:]+1)]     * ( 1. - xd[:] )
           + val[ind1d(i[:]+1,j[:],k[:]+1)]   * xd[:];
    c11[:] = val[ind1d(i[:],j[:]+1,k[:]+1)]   * ( 1. - xd[:] )
           + val[ind1d(i[:]+1,j[:]+1,k[:]+1)] * xd[:];

    c0[:] = c00[:] * ( 1. - yd[:] ) + c10[:] * yd[:];
    c1[:] = c01[:] * ( 1. - yd[:] ) + c11[:] * yd[:];

    ret[0:size] = c0[:] * ( 1. - zd[:] ) + c1[:] * zd[:];

}*/

#endif
