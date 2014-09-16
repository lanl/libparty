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

#include <math.h>

//TODO needs more accurate PI
inline double cubic_spline(double r, double h) {
  double rnorm = r/h;
  if (rnorm > 2.0) return 0.0;
  double c_n = 1/3.1415/h/h/h;
  if (rnorm > 1.0) return c_n*0.25*pow((2-rnorm),3);
  return c_n*(1.0 - 1.5*rnorm*rnorm + 0.75*rnorm*rnorm*rnorm);
}

inline double Dcubic_spline(double r, double h) {
  double rnorm = r/h;
  if (rnorm > 2.0) return 0.0;
  double c_n = 1/3.1415/h/h/h;
  if (rnorm > 1.0) return c_n*-3.0*0.25 * pow((rnorm-2),2);
  return c_n*0.25*(-12.0*rnorm + 9.0*rnorm*rnorm);
}

inline double Dhcubic_spline(double r, double h) {
  double rnorm = r/h;
  if (rnorm > 2.0) return 0.0;
  double c_n = 1/3.1415/h/h/h;
  if (rnorm > 1.0) return c_n*3.0*0.25 *rnorm / h * pow(2.0-rnorm,2.0) + 3.0*c_n/h*cubic_spline(r,h);
  return c_n*(3.0*rnorm*rnorm/h-9.0*0.25*rnorm*rnorm*rnorm/h) + 3.0*c_n/h*cubic_spline(r,h);
}

