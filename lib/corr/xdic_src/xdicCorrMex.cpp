#include "xdic_common.hpp"
#include <unordered_map>
using namespace xdic;
struct Image {
    const void* data;
    bool single;
    mwSize height,width;
    double at(mwSize i) const {
        return single ? double(static_cast<const float*>(data)[i]) : static_cast<const double*>(data)[i];
    }
    static double cubic(double vm,double v0,double v1,double v2,double t) {
        return vm*(-t)*(1-t)*(2-t)/6+v0*(1+t)*(1-t)*(2-t)/2+
            v1*(1+t)*t*(2-t)/2+v2*(1+t)*t*(t-1)/6;
    }
    double sample(double x,double y,bool useCubic) const {
        mwSize ix=mwSize(std::floor(x)),iy=mwSize(std::floor(y));
        double sx=x-ix,sy=y-iy;
        if(!useCubic) {
            mwSize o=ix*height+iy;
            return (1-sx)*(1-sy)*at(o)+sx*(1-sy)*at(o+height)+
                (1-sx)*sy*at(o+1)+sx*sy*at(o+height+1);
        }
        double row[4];
        for(int j=0;j<4;++j) {
            mwSize o=(ix-1)*height+(iy-1+j);
            row[j]=cubic(at(o),at(o+height),at(o+2*height),at(o+3*height),sx);
        }
        return cubic(row[0],row[1],row[2],row[3],sy);
    }
};
Image image(const mxArray* a,mwSize h,mwSize w) {
    require(!mxIsComplex(a) && !mxIsSparse(a) && mxGetNumberOfDimensions(a)==2 &&
        (mxIsSingle(a)||mxIsDouble(a)) && mxGetM(a)==h && mxGetN(a)==w,
        "Images must be real full single or double arrays of the prepared size.");
    return {mxGetData(a),mxIsSingle(a),h,w};
}
void correlateKernel(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    require(nrhs==6 && nlhs==1,"state=xdicCorrMex(prep,f,g,a,interpolation,assemble).");
    const mxArray* prep=prhs[0];
    require(scalar(field(prep,"version"))==1,"Unsupported preparation version.");
    const mxArray* sz=field(prep,"imageSize");
    const double* dims=doubles(sz);
    require(mxGetNumberOfElements(sz)==2,"Invalid prepared image size.");
    mwSize h=integer(dims[0],100000),w=integer(dims[1],100000);
    require(h>=5 && w>=5 && h*w<=100000000,"Invalid prepared image dimensions.");
    Image reference=image(prhs[1],h,w),deformed=image(prhs[2],h,w);
    const mxArray* basis=field(prep,"basis");
    require(mxIsSparse(basis) && mxIsDouble(basis) && !mxIsComplex(basis),"Preparation basis must be real double sparse.");
    mwSize P=mxGetM(basis),nb=mxGetN(basis),ndof=2*nb;
    require(nb>0 && nb<100000000,"Invalid basis dimension.");
    mwSize nn=integer(scalar(field(prep,"nNodes")),nb);
    require(nn>0,"Invalid node count.");
    const mxArray* pixels=field(prep,"pixelIndex");
    const double* pixel=doubles(pixels);
    require(mxGetNumberOfElements(pixels)==P,"Pixel list and basis rows disagree.");
    const double* a=doubles(prhs[3]);
    require(mxGetNumberOfElements(prhs[3])==ndof,"Expected interleaved ux,uy coefficients for every basis function.");
    for(mwSize i=0;i<ndof;++i) require(std::isfinite(a[i]),"Nonfinite displacement coefficient.");
    require(mxIsChar(prhs[4]),"Interpolation must be 'linear' or 'cubic'.");
    char* text=mxArrayToString(prhs[4]);
    require(text!=nullptr,"Could not read interpolation.");
    std::string interpolation(text);
    mxFree(text);
    require(interpolation=="linear" || interpolation=="cubic","Interpolation must be 'linear' or 'cubic'.");
    bool cubic=interpolation=="cubic";
    require(mxIsLogicalScalar(prhs[5]),"assemble must be a scalar logical.");
    bool assemble=mxIsLogicalScalarTrue(prhs[5]);
    const mwIndex* ir=mxGetIr(basis);
    const mwIndex* jc=mxGetJc(basis);
    const double* values=mxGetDoubles(basis);
    mwSize nz=jc[nb];
    std::vector<mwSize> offset(P+1,0),cursor;
    for(mwSize k=0;k<nz;++k) {
        require(ir[k]<P && std::isfinite(values[k]),"Invalid sparse basis.");
        ++offset[ir[k]+1];
    }
    for(mwSize p=0;p<P;++p) offset[p+1]+=offset[p];
    cursor=offset;
    std::vector<mwIndex> columns(nz);
    std::vector<double> rowValues(nz);
    for(mwSize c=0;c<nb;++c)
        for(mwSize k=jc[c];k<jc[c+1];++k) {
            mwSize dest=cursor[ir[k]]++;
            columns[dest]=c;
            rowValues[dest]=values[k];
        }
    std::vector<double> rhs(ndof,0),nodal(nn,0);
    mxArray* residual=mxCreateDoubleMatrix(h,w,mxREAL);
    double* R=mxGetDoubles(residual);
    std::fill(R,R+h*w,mxGetNaN());
    mwSize shape[3]={h,w,2};
    mxArray* displacement=mxCreateNumericArray(3,shape,mxDOUBLE_CLASS,mxREAL);
    double* U=mxGetDoubles(displacement);
    std::fill(U,U+2*h*w,mxGetNaN());
    mxArray* valid=mxCreateLogicalMatrix(h,w);
    mxLogical* validData=mxGetLogicals(valid);
    std::unordered_map<uint64_t,double> matrix;
    if(assemble) matrix.reserve(nb*60);
    double squareResidual=0;
    mwSize used=0,gradientUsed=0;
    for(mwSize p=0;p<P;++p) {
        mwSize index1=integer(pixel[p],h*w);
        require(index1>0,"Pixel indices must be one based.");
        mwSize index=index1-1,y=index%h,x=index/h;
        require(x>0 && x+1<w && y>0 && y+1<h,"Prepared pixel violates reference-gradient boundary margin.");
        double ux=0,uy=0;
        for(mwSize k=offset[p];k<offset[p+1];++k) {
            ux+=rowValues[k]*a[2*columns[k]];
            uy+=rowValues[k]*a[2*columns[k]+1];
        }
        U[index]=ux;
        U[index+h*w]=uy;
        double qx=double(x)+ux,qy=double(y)+uy,margin=cubic ? 3.0 : 2.0;
        if(!(qx>=margin && qy>=margin && qx<double(w)-margin && qy<double(h)-margin)) continue;
        double r=reference.at(index)-deformed.sample(qx,qy,cubic);
        if(!std::isfinite(r)) continue;
        R[index]=r;
        validData[index]=true;
        squareResidual+=r*r;
        ++used;
        for(mwSize k=offset[p];k<offset[p+1];++k)
            if(columns[k]<nn) nodal[columns[k]]+=rowValues[k]*r*r;
        if(!assemble) continue;
        double gx=0.5*(reference.at(index+h)-reference.at(index-h));
        double gy=0.5*(reference.at(index+1)-reference.at(index-1));
        if(!std::isfinite(gx) || !std::isfinite(gy)) continue;
        ++gradientUsed;
        std::vector<mwIndex> dofs;
        std::vector<double> J;
        dofs.reserve(2*(offset[p+1]-offset[p]));
        J.reserve(dofs.capacity());
        for(mwSize k=offset[p];k<offset[p+1];++k) {
            dofs.push_back(2*columns[k]);
            J.push_back(rowValues[k]*gx);
            dofs.push_back(2*columns[k]+1);
            J.push_back(rowValues[k]*gy);
        }
        for(mwSize i=0;i<dofs.size();++i) {
            rhs[dofs[i]]+=J[i]*r;
            for(mwSize j=0;j<=i;++j) {
                double v=J[i]*J[j];
                if(v!=0) {
                    mwSize row=std::min(dofs[i],dofs[j]),col=std::max(dofs[i],dofs[j]);
                    matrix[uint64_t(col)*ndof+row]+=v;
                }
            }
        }
    }
    std::vector<Entry> entries;
    entries.reserve(matrix.size()*2);
    for(const auto& kv:matrix) {
        mwSize row=mwSize(kv.first%ndof),col=mwSize(kv.first/ndof);
        entries.push_back({row,col,kv.second});
        if(row!=col) entries.push_back({col,row,kv.second});
    }
    const char* names[]={"H","b","R","r","nodalResidual","U","validMask","validCount","gradientCount"};
    plhs[0]=mxCreateStructMatrix(1,1,9,names);
    set(plhs[0],"H",sparse(ndof,ndof,entries));
    set(plhs[0],"b",column(rhs));
    set(plhs[0],"R",residual);
    set(plhs[0],"r",mxCreateDoubleScalar(used ? std::sqrt(squareResidual/used) : mxGetNaN()));
    set(plhs[0],"nodalResidual",column(nodal));
    set(plhs[0],"U",displacement);
    set(plhs[0],"validMask",valid);
    set(plhs[0],"validCount",mxCreateDoubleScalar(double(used)));
    set(plhs[0],"gradientCount",mxCreateDoubleScalar(double(gradientUsed)));
}
void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    try { correlateKernel(nlhs,plhs,nrhs,prhs); }
    catch(const std::exception& e) { mexErrMsgIdAndTxt("xdic:Input","%s",e.what()); }
}
