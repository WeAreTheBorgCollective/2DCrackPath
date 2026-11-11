#include "xdic_common.hpp"
using namespace xdic;
void prepare(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    require(nrhs==4 && nlhs==1,"prep=xdicPrepMex(nodesXY,triangles,imageSize,crackXY).");
    const double* pos=doubles(prhs[0]);
    const double* con=doubles(prhs[1]);
    const double* size=doubles(prhs[2]);
    const double* path=doubles(prhs[3]);
    mwSize nn=mxGetM(prhs[0]),ne=mxGetM(prhs[1]);
    require(nn>=3 && mxGetN(prhs[0])==2 && ne>0 && mxGetN(prhs[1])==3,"Expected 2D nodes and T3 connectivity.");
    require(mxGetNumberOfElements(prhs[2])==2,"imageSize must be [rows columns].");
    mwSize height=integer(size[0],100000),width=integer(size[1],100000);
    require(height>=5 && width>=5 && height*width<=100000000,"Invalid image size.");
    mwSize np=height*width,nc=mxGetM(prhs[3]);
    if(mxIsEmpty(prhs[3])) nc=0;
    else require(mxGetN(prhs[3])==2 && nc>=2,"Crack must be an open ordered Nx2 polyline or empty.");
    std::vector<Point> nodes(nn),crack(nc);
    for(mwSize i=0;i<nn;++i) {
        nodes[i]={pos[i],pos[i+nn]};
        require(std::isfinite(nodes[i].x) && std::isfinite(nodes[i].y),"Nonfinite mesh node.");
    }
    for(mwSize i=0;i<nc;++i) {
        crack[i]={path[i],path[i+nc]};
        require(std::isfinite(crack[i].x) && std::isfinite(crack[i].y),"Nonfinite crack point.");
        if(i) {
            double l2=dot(subtract(crack[i],crack[i-1]),subtract(crack[i],crack[i-1]));
            require(std::isfinite(l2) && l2>1e-16,"Collapsed crack segment.");
        }
    }
    for(mwSize i=1;i<nc;++i)
        for(mwSize j=i+2;j<nc;++j)
            require(!intersect(crack[i-1],crack[i],crack[j-1],crack[j]),"Self-intersecting crack paths are unsupported.");
    std::vector<std::array<mwIndex,3>> triangles(ne);
    std::vector<bool> candidate(nn,false),tip(nn,false),cut(ne,false);
    std::vector<double> cutList,tipList,nodeH(nn,0);
    for(mwSize i=0;i<nn;++i) if(nc) nodeH[i]=heaviside(nodes[i],crack);
    for(mwSize e=0;e<ne;++e) {
        std::array<Point,3> t;
        for(int j=0;j<3;++j) {
            mwSize id=integer(con[e+ne*j],nn);
            require(id>0,"Connectivity must be one based.");
            triangles[e][j]=id-1;
            t[j]=nodes[id-1];
        }
        double area2=orient(t[0],t[1],t[2]);
        require(std::isfinite(area2) && std::abs(area2)>1e-10,"Degenerate triangle.");
        if(nc) {
            bool containsTip=inside(crack.front(),t) || inside(crack.back(),t);
            for(mwSize k=1;k<nc;++k)
                if(segmentTouchesTriangle(crack[k-1],crack[k],t)) { cut[e]=true; break; }
            for(int j=0;j<3;++j) {
                if(containsTip) tip[triangles[e][j]]=true;
                if(cut[e]) candidate[triangles[e][j]]=true;
            }
            if(cut[e]) cutList.push_back(double(e+1));
        }
    }
    // Unique ownership of shared-edge pixels; no duplicate quadrature.
    std::vector<int> owner(np,-1);
    std::vector<double> pixel,element,pixelH,plus(nn,0),minus(nn,0);
    std::vector<std::array<double,3>> shape;
    for(mwSize e=0;e<ne;++e) {
        std::array<Point,3> t={{nodes[triangles[e][0]],nodes[triangles[e][1]],nodes[triangles[e][2]]}};
        double xmin=std::min({t[0].x,t[1].x,t[2].x}),xmax=std::max({t[0].x,t[1].x,t[2].x});
        double ymin=std::min({t[0].y,t[1].y,t[2].y}),ymax=std::max({t[0].y,t[1].y,t[2].y});
        double lx=std::max(2.0,std::ceil(xmin-1e-9)),hx=std::min(double(width-1),std::floor(xmax+1e-9));
        double ly=std::max(2.0,std::ceil(ymin-1e-9)),hy=std::min(double(height-1),std::floor(ymax+1e-9));
        if(lx>hx || ly>hy) continue;
        for(mwSize x=mwSize(lx);x<=mwSize(hx);++x) {
            for(mwSize y=mwSize(ly);y<=mwSize(hy);++y) {
                mwSize index=(x-1)*height+y-1;
                if(owner[index]>=0) continue;
                auto n=bary({double(x),double(y)},t);
                if(*std::min_element(n.begin(),n.end()) < -1e-9) continue;
                double sum=0;
                for(double& v:n) { v=std::max(0.0,v); sum+=v; }
                for(double& v:n) v/=sum;
                owner[index]=int(e);
                double H=nc ? heaviside({double(x),double(y)},crack) : 0;
                pixel.push_back(double(index+1));
                element.push_back(double(e+1));
                pixelH.push_back(H);
                shape.push_back(n);
                for(int j=0;j<3;++j) {
                    mwIndex node=triangles[e][j];
                    (H ? plus[node] : minus[node])+=n[j]*n[j];
                }
            }
        }
    }
    require(!pixel.empty(),"Mesh contains no eligible interior image pixels.");
    std::vector<double> enriched,excluded,geometric,index(nn,0);
    for(mwSize i=0;i<nn;++i) {
        if(tip[i]) tipList.push_back(double(i+1));
        if(candidate[i] && !tip[i]) {
            geometric.push_back(double(i+1));
            if(plus[i]>1e-10 && minus[i]>1e-10) {
                enriched.push_back(double(i+1));
                index[i]=double(enriched.size());
            } else excluded.push_back(double(i+1));
        }
    }
    mwSize count=pixel.size(),nb=nn+enriched.size();
    std::vector<Entry> entries;
    entries.reserve(count*6);
    mxArray* shapeOut=mxCreateDoubleMatrix(count,3,mxREAL);
    double* shapeData=mxGetDoubles(shapeOut);
    for(mwSize p=0;p<count;++p) {
        auto t=triangles[mwSize(element[p])-1];
        for(int j=0;j<3;++j) {
            shapeData[p+count*j]=shape[p][j];
            mwSize node=t[j];
            if(shape[p][j]!=0) entries.push_back({p,node,shape[p][j]});
            double shifted=shape[p][j]*(pixelH[p]-nodeH[node]);
            if(index[node]>0 && shifted!=0)
                entries.push_back({p,nn+mwSize(index[node])-1,shifted});
        }
    }
    const char* names[]={"basis","pixelIndex","elementIndex","shape","imageSize","nNodes",
        "enrichedNodes","tipSupportNodes","cutElements","fullyCutCandidates","unresolvedNodes",
        "pixelHeaviside","nodeHeaviside","crackPathXY","version"};
    plhs[0]=mxCreateStructMatrix(1,1,15,names);
    set(plhs[0],"basis",sparse(count,nb,entries));
    set(plhs[0],"pixelIndex",column(pixel));
    set(plhs[0],"elementIndex",column(element));
    set(plhs[0],"shape",shapeOut);
    set(plhs[0],"imageSize",mxDuplicateArray(prhs[2]));
    set(plhs[0],"nNodes",mxCreateDoubleScalar(double(nn)));
    set(plhs[0],"enrichedNodes",column(enriched));
    set(plhs[0],"tipSupportNodes",column(tipList));
    set(plhs[0],"cutElements",column(cutList));
    set(plhs[0],"fullyCutCandidates",column(geometric));
    set(plhs[0],"unresolvedNodes",column(excluded));
    set(plhs[0],"pixelHeaviside",column(pixelH));
    set(plhs[0],"nodeHeaviside",column(nodeH));
    set(plhs[0],"crackPathXY",mxDuplicateArray(prhs[3]));
    set(plhs[0],"version",mxCreateDoubleScalar(1));
}
void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    try { prepare(nlhs,plhs,nrhs,prhs); }
    catch(const std::exception& e) { mexErrMsgIdAndTxt("xdic:Input","%s",e.what()); }
}
