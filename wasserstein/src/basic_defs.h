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

#ifndef BASIC_DEFS_H
#define BASIC_DEFS_H

#include <vector>
#include <math.h>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <string>
#include <assert.h>

#include "def_debug.h"

#define MIN_VALID_ID 10

struct Point {
    double x, y;
    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;
    Point(double ax, double ay) : x(ax), y(ay) {}
    Point() : x(0.0), y(0.0) {}
    friend std::ostream& operator<<(std::ostream& output, const Point p);
};

struct DiagramPoint 
{
    // Points above the diagonal have type NORMAL
    // Projections onto the diagonal have type DIAG
    // for DIAG points only x-coordinate is relevant
    // to-do: add getters/setters, checks in constructors, etc
    enum Type { NORMAL, DIAG};
    // data members
    double x, y;
    Type type;
    size_t id;
    // for diagonal points projId is the id of the normal point whose
    // projection they are
    // for normal points projId is the id of the projection 
    size_t projId;
    // operators, constructors
    bool operator==(const DiagramPoint& other) const;
    bool operator!=(const DiagramPoint& other) const;
    DiagramPoint(double xx, double yy, Type ttype, size_t uid, size_t projId);
    DiagramPoint(double xx, double yy, Type ttype, size_t uid);
    //DiagramPoint(double xx, double yy); // type is NORMAL 
    //DiagramPoint() : x(0.0), y(0.0), type(DIAG) {}
    bool isDiagonal(void) const { return type == DIAG; }
    bool isNormal(void) const { return type == NORMAL; }
    Point getOrdinaryPoint() const; // for diagonal points return the coords or projection
    double getRealX() const; // return the x-coord
    double getRealY() const; // return the y-coord
    friend std::ostream& operator<<(std::ostream& output, const DiagramPoint p);

    struct LexicographicCmp
    {
        bool    operator()(const DiagramPoint& p1, const DiagramPoint& p2) const
        { return p1.type < p2.type || (p1.type == p2.type && (p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y))); }
    };
};

struct PointHash {
    size_t operator()(const Point& p) const{
        return std::hash<double>()(p.x)^std::hash<double>()(p.y);
    }
};

struct DiagramPointHash {
    size_t operator()(const DiagramPoint& p) const{
        //return std::hash<double>()(p.x)^std::hash<double>()(p.y)^std::hash<bool>()(p.type == DiagramPoint::NORMAL);
        assert(p.id >= MIN_VALID_ID);
        return std::hash<int>()(p.id);
    }
};

double sqrDist(const Point& a, const Point& b);
double dist(const Point& a, const Point& b);
double distLInf(const DiagramPoint& a, const DiagramPoint& b);
double distLp(const DiagramPoint& a, const DiagramPoint& b, const double p);

typedef std::unordered_set<Point, PointHash> PointSet;

class DiagramPointSet {
public:
    void insert(const DiagramPoint p);
    void erase(const DiagramPoint& p, bool doCheck = true); // if doCheck, erasing non-existing elements causes assert
    void erase(const std::unordered_set<DiagramPoint, DiagramPointHash>::const_iterator it);
    size_t size() const;
    void reserve(const size_t newSize);
    void clear();
    bool empty() const;
    bool hasElement(const DiagramPoint& p) const;
    bool operator==(const DiagramPointSet& Other) const;
    std::unordered_set<DiagramPoint, DiagramPointHash>::iterator find(const DiagramPoint& p) { return points.find(p); };
    std::unordered_set<DiagramPoint, DiagramPointHash>::const_iterator find(const DiagramPoint& p) const { return points.find(p); };
    std::unordered_set<DiagramPoint, DiagramPointHash>::iterator begin() { return points.begin(); };
    std::unordered_set<DiagramPoint, DiagramPointHash>::iterator end() { return points.end(); }
    std::unordered_set<DiagramPoint, DiagramPointHash>::const_iterator cbegin() const { return points.cbegin(); }
    std::unordered_set<DiagramPoint, DiagramPointHash>::const_iterator cend() const { return points.cend(); }
    friend std::ostream& operator<<(std::ostream& output, const DiagramPointSet& ps);
private:
    std::unordered_set<DiagramPoint, DiagramPointHash> points;
};

template<typename DiagPointContainer>
double getFurthestDistance3Approx(DiagPointContainer& A, DiagPointContainer& B)
{
    double result { 0.0 };
    DiagramPoint begA = *(A.begin());
    DiagramPoint optB = *(B.begin());
    for(const auto& pointB : B) {
        if (distLInf(begA, pointB) > result) {
            result = distLInf(begA, pointB);
            optB = pointB;
        }
    }
    for(const auto& pointA : A) {
        if (distLInf(pointA, optB) > result) {
            result = distLInf(pointA, optB);
        }
    }
    return result;
}


#endif
