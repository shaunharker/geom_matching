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

#include <algorithm>
#include <cfloat>
#include <set>
#include "basic_defs.h"

// Point

bool Point::operator==(const Point& other) const
{
    return ((this->x == other.x) and (this->y == other.y));
}

bool Point::operator!=(const Point& other) const
{
    return !(*this == other);
}

std::ostream& operator<<(std::ostream& output, const Point p)
{
    output << "(" << p.x << ", " << p.y << ")";
    return output;
}

std::ostream& operator<<(std::ostream& output, const PointSet& ps)
{
    output << "{ ";
    for(auto& p : ps) {
        output << p << ", ";
    }
    output << "\b\b }";
    return output;
}

double sqrDist(const Point& a, const Point& b)
{
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

double dist(const Point& a, const Point& b)
{
    return sqrt(sqrDist(a, b));
}

// DiagramPoint

// compute l-inf distance between two diagram points
double distLInf(const DiagramPoint& a, const DiagramPoint& b)
{
    if ( DiagramPoint::NORMAL == a.type &&
         DiagramPoint::NORMAL == b.type ) {
        // distance between normal points is a usual l-inf distance
        return std::max(fabs(a.x - b.x), fabs(a.y - b.y));
    } else if ( DiagramPoint::DIAG == a.type &&
                DiagramPoint::DIAG == b.type ) {
        // distance between points on the diagonal is 0
        return 0;
    } else {
        if ( DiagramPoint::DIAG == a.type ) {
            // real coordinates of a are ( (a.x + a.y)/2, (a.x + a.y)/2 ) 
            return std::max(fabs(0.5* (a.x + a.y) - b.x), 
                            fabs(0.5 * (a.x + a.y) - b.y));
        } else {
            // real coordinates of b are ( (b.x + b.y)/2, (b.x + b.y)/2 ) 
            return std::max(fabs(a.x - 0.5 * (b.x + b.y) ), 
                            fabs(a.y - 0.5 * ( b.y + b.x)));
        }
    }
}

bool DiagramPoint::operator==(const DiagramPoint& other) const
{
    assert(this->id >= MIN_VALID_ID);
    assert(other.id >= MIN_VALID_ID);
    bool areEqual{ this->id == other.id };
    assert(!areEqual or  ((this->x == other.x) and (this->y == other.y) and (this->type == other.type)));
    return areEqual;
}

bool DiagramPoint::operator!=(const DiagramPoint& other) const
{
    return !(*this == other);
}

std::ostream& operator<<(std::ostream& output, const DiagramPoint p)
{
    if ( p.type == DiagramPoint::DIAG ) {
        output << "(" << p.x << ", " << p.y << ", " <<  0.5 * (p.x + p.y) << ", "  << p.id <<", " << p.projId << " DIAG )";
    } else {
        output << "(" << p.x << ", " << p.y << ", " << p.id <<", " << p.projId << " NORMAL)";
    }
    return output;
}

std::ostream& operator<<(std::ostream& output, const DiagramPointSet& ps)
{
    output << "{ ";
    for(auto pit = ps.cbegin(); pit != ps.cend(); ++pit) {
        output << *pit << ", ";
    }
    output << "\b\b }";
    return output;
}

DiagramPoint::DiagramPoint(double xx, double yy, Type ttype, size_t uid, size_t _projId) : 
    x(xx),
    y(yy),
    type(ttype),
    id(uid),
    projId(_projId)
{
    if ( yy < xx )
        throw "Point is below the diagonal";
    if ( yy == xx and ttype != DiagramPoint::DIAG)
        throw "Point on the main diagonal must have DIAG type";

}


DiagramPoint::DiagramPoint(double xx, double yy, Type ttype, size_t uid) : 
    x(xx),
    y(yy),
    type(ttype),
    id(uid),
    projId(MIN_VALID_ID)
{
    if ( yy < xx )
        throw "Point is below the diagonal";
    if ( yy == xx and ttype != DiagramPoint::DIAG)
        throw "Point on the main diagonal must have DIAG type";

}

/*
DiagramPoint::DiagramPoint(double xx, double yy) : 
    x(xx),
    y(yy), 
    type(NORMAL)
{
    if ( xx < 0 )
        throw "Negative x coordinate";
    if ( yy < 0)
        throw "Negative y coordinate";
    if ( yy <= xx )
        throw "Point is below the diagonal";
}
*/

void DiagramPointSet::insert(const DiagramPoint p)
{
    points.insert(p);
}

Point DiagramPoint::getOrdinaryPoint() const
{
    if (DiagramPoint::NORMAL == type)
        return Point(x, y);
    else if (DiagramPoint::DIAG == type)
        return Point(0.5 * ( x + y), 0.5 * (x+y) );
    else
        assert(false);
}

double DiagramPoint::getRealX() const
{
    if (DiagramPoint::NORMAL == type)
        return x;
    else 
        return 0.5 * ( x + y);
}

double DiagramPoint::getRealY() const
{
    if (DiagramPoint::NORMAL == type)
        return y;
    else 
        return 0.5 * ( x + y);
}

// erase should be called only for the element of the set
void DiagramPointSet::erase(const DiagramPoint& p, bool doCheck)
{
    auto it = points.find(p);
    if (it != points.end()) {
        points.erase(it);
    } else {
        assert(!doCheck);
    }
}

void DiagramPointSet::reserve(const size_t newSize)
{
    points.reserve(newSize);
}


void DiagramPointSet::erase(const std::unordered_set<DiagramPoint, DiagramPointHash>::const_iterator it)
{
    points.erase(it);
}

void DiagramPointSet::clear()
{
    points.clear();
}

size_t DiagramPointSet::size() const
{
    return points.size();
}

bool DiagramPointSet::empty() const
{
    return points.empty();
}

bool DiagramPointSet::hasElement(const DiagramPoint& p) const
{
    return points.find(p) != points.end();
}

bool DiagramPointSet::operator==(const DiagramPointSet& other) const
{
    if (size() != other.size())
        return false;

    std::set<DiagramPoint, DiagramPoint::LexicographicCmp>      dgm1(this->cbegin(), this->cend()),
                                                                dgm2(other.cbegin(), other.cend());

    for (auto& p : other.points)
        if (dgm1.find(p) == dgm1.end())
            return false;

    for (auto& p : points)
        if (dgm2.find(p) == dgm2.end())
            return false;

    return true;
}
