#ifndef XDIC_COMMON_HPP
#define XDIC_COMMON_HPP
#include "mex.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace xdic {
inline void require(bool ok,const char* message) {
    if(!ok) throw std::runtime_error(message);
}
inline const mxArray* field(const mxArray* s,const char* name) {
    require(mxIsStruct(s) && mxGetNumberOfElements(s)==1,"Expected a scalar preparation structure.");
    const mxArray* a=mxGetField(s,0,name);
    require(a!=nullptr,"Missing preparation field.");
    return a;
}
inline const double* doubles(const mxArray* a) {
    require(mxIsDouble(a) && !mxIsComplex(a) && !mxIsSparse(a),"Expected a real full double array.");
    return mxGetDoubles(a);
}
inline double scalar(const mxArray* a) {
    const double* p=doubles(a);
    require(mxGetNumberOfElements(a)==1 && std::isfinite(p[0]),"Expected a finite scalar.");
    return p[0];
}
inline mwSize integer(double v,mwSize maximum) {
    require(std::isfinite(v) && v>=0 && v==std::floor(v) && v<=double(maximum),"Invalid integer index or dimension.");
    return mwSize(v);
}
inline mxArray* column(const std::vector<double>& v) {
    mxArray* a=mxCreateDoubleMatrix(v.size(),1,mxREAL);
    std::copy(v.begin(),v.end(),mxGetDoubles(a));
    return a;
}
inline void set(mxArray* s,const char* name,mxArray* a) {
    mxSetField(s,0,name,a);
}
struct Point { double x,y; };
inline Point subtract(Point a,Point b) { return {a.x-b.x,a.y-b.y}; }
inline double cross(Point a,Point b) { return a.x*b.y-a.y*b.x; }
inline double dot(Point a,Point b) { return a.x*b.x+a.y*b.y; }
inline double orient(Point a,Point b,Point c) { return cross(subtract(b,a),subtract(c,a)); }
inline std::array<double,3> bary(Point p,const std::array<Point,3>& t) {
    double den=orient(t[0],t[1],t[2]);
    double n1=orient(p,t[1],t[2])/den;
    double n2=orient(t[0],p,t[2])/den;
    return {{n1,n2,1-n1-n2}};
}
inline bool inside(Point p,const std::array<Point,3>& t) {
    auto n=bary(p,t);
    return *std::min_element(n.begin(),n.end())>=-1e-9;
}
inline bool segmentTouchesTriangle(Point a,Point b,const std::array<Point,3>& t) {
    auto na=bary(a,t);
    auto nb=bary(b,t);
    double lo=0,hi=1;
    for(int j=0;j<3;++j) {
        double delta=nb[j]-na[j];
        if(std::abs(delta)<1e-14) {
            if(na[j]<-1e-10) return false;
        } else {
            double limit=(-1e-10-na[j])/delta;
            if(delta>0) lo=std::max(lo,limit);
            else hi=std::min(hi,limit);
        }
    }
    return hi>lo && (hi-lo)*std::sqrt(dot(subtract(b,a),subtract(b,a)))>1e-8;
}
inline bool onSegment(Point a,Point b,Point p) {
    return std::abs(orient(a,b,p))<1e-9 &&
        p.x>=std::min(a.x,b.x)-1e-9 && p.x<=std::max(a.x,b.x)+1e-9 &&
        p.y>=std::min(a.y,b.y)-1e-9 && p.y<=std::max(a.y,b.y)+1e-9;
}
inline bool intersect(Point a,Point b,Point c,Point d) {
    double p=orient(a,b,c),q=orient(a,b,d),r=orient(c,d,a),s=orient(c,d,b);
    return ((p>0)!=(q>0) && (r>0)!=(s>0)) ||
        onSegment(a,b,c) || onSegment(a,b,d) || onSegment(c,d,a) || onSegment(c,d,b);
}
inline double heaviside(Point p,const std::vector<Point>& crack) {
    double best=std::numeric_limits<double>::infinity(),side=1;
    for(mwSize k=1;k<crack.size();++k) {
        Point v=subtract(crack[k],crack[k-1]);
        Point w=subtract(p,crack[k-1]);
        double t=std::max(0.0,std::min(1.0,dot(w,v)/dot(v,v)));
        Point e={w.x-t*v.x,w.y-t*v.y};
        double distance=dot(e,e);
        if(distance<best) {
            best=distance;
            side=cross(v,w)>=0 ? 1.0 : 0.0;
        }
    }
    return side;
}
struct Entry { mwIndex row,col; double value; };
inline mxArray* sparse(mwSize rows,mwSize cols,std::vector<Entry>& entries) {
    std::sort(entries.begin(),entries.end(),[](const Entry& a,const Entry& b) {
        return a.col<b.col || (a.col==b.col && a.row<b.row);
    });
    mxArray* out=mxCreateSparse(rows,cols,entries.size(),mxREAL);
    mwIndex* ir=mxGetIr(out);
    mwIndex* jc=mxGetJc(out);
    double* values=mxGetDoubles(out);
    mwSize offset=0;
    for(mwSize c=0;c<cols;++c) {
        jc[c]=offset;
        while(offset<entries.size() && entries[offset].col==c) {
            ir[offset]=entries[offset].row;
            values[offset]=entries[offset].value;
            ++offset;
        }
    }
    jc[cols]=offset;
    return out;
}
}
#endif
