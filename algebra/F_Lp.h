/*
 *  Copyright (c) 2010, INRIA, Project ALICE
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://alice.loria.fr
 *
 *     ALICE Project
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef __GEEX_LPP_FAST__
#define __GEEX_LPP_FAST__

#include "../common/types.h"
#include <vector>

namespace Geex {

    class Mesh ; 

    double compute_F_Lp(
        bool volumic,                  // IN: false for surface meshing, true for volume meshing
        unsigned int p,                // IN: Lp norm to be used
        Mesh* mesh,                    // IN: the PLC (required to resolve configurations with coplanar
                                       //   facets, see VoroFunc::set_intersection).
        const std::vector<int>& sym,   // IN: 10 integers per integration simplex:
                                       //   - index of x0
                                       //   - symbolic representation of C1,C2,C3 (3 integers each)
        const std::vector<vec3>& C,    // IN: C vertices of integration simplices (3 per integration simplex)
                                       //   (since they are computed during the combinatorial phase, they are
                                       //    kept and re-used here).
        const std::vector<vec3>& X,    // IN: vertices
        const std::vector<plane3>& Q,  // IN: boundary facets
        const std::vector<mat3>& M,    // IN: anisotropy matrix (one per integration simplex)
        std::vector<double>& g         // OUT: gradient
    ) ;

}

#endif
