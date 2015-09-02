/*
 
Copyright (c) 2015, M. Kerber, D. Morozov, A. Nigmetov
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
(Enhancements) to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to copyright holder,
without imposing a separate written license agreement for such Enhancements,
then you hereby grant the following license: a  non-exclusive, royalty-free
perpetual license to install, use, modify, prepare derivative works, incorporate
into other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.

  */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <random>

#include "basic_defs.h"
#include "wasserstein.h"

bool readDiagramPointSets(const char* fnameA,
                          const char* fnameB,
                          DiagramPointSet& A,
                          DiagramPointSet& B)
{
    std::ifstream fA(fnameA);
    if (!fA.good()) {
        std::cerr << "Cannot open file " << fnameA << std::endl;
        return false;
    }
    std::ifstream fB(fnameB);
    if (!fB.good()) {
        std::cerr << "Cannot open file " << fnameB << std::endl;
        return false;
    }
    A.clear();
    B.clear();
    Point p;
    size_t uniqueId {MIN_VALID_ID};
    size_t projId;
    while(fA >> p.x >> p.y) {
        // normal point, its projection will be added next, so +1
        projId = uniqueId+1;
        DiagramPoint dpA {p.x, p.y, DiagramPoint::NORMAL, uniqueId++, projId};
        // diagonal point, its parent has been added, so -1
        projId = uniqueId-1;
        DiagramPoint dpB {p.x, p.y, DiagramPoint::DIAG, uniqueId++, projId};
        A.insert(dpA);
        B.insert(dpB);
    }
    fA.close();
    while(fB >> p.x >> p.y) {
        // normal point, its projection will be added next, so +1
        projId = uniqueId+1;
        DiagramPoint dpB {p.x, p.y, DiagramPoint::NORMAL, uniqueId++, projId};
        // diagonal point, its parent has been added, so -1
        projId = uniqueId-1;
        DiagramPoint dpA {p.x, p.y, DiagramPoint::DIAG, uniqueId++, projId};
        B.insert(dpB);
        A.insert(dpA);
    }
    fB.close();
    return true;
}

int main(int argc, char* argv[])
{
    DiagramPointSet A, B;
    if (argc < 3 ) {
        std::cerr << "Usage: " << argv[0] << " file1 file2 [wasserstein_power] [relative_error]. By default power is 1.0, relative error is 0.01" << std::endl;
        return 1;
    }
    if (!readDiagramPointSets(argv[1], argv[2], A, B)) {
        std::exit(1);
    }

    double wasserPower = (4 <= argc) ? atof(argv[3]) : 1.0;
    //default : 1% error
    double delta = (5 <= argc) ? atof(argv[4]) : 0.01;
    if ( delta <= 0.0) {
        std::cerr << "Relative error must be positive" << std::endl;
        return 0;
    }
    double res = wassersteinDistSlow(A, B, wasserPower, delta);
    
    std::cout << std::setprecision(15) << res << std::endl;
    return 0;
}
