#include "morpher.hh"
#include <iostream>
#include <igl/readOBJ.h>
#include <igl/mat_max.h>
#include <igl/dot.h>
#include <igl/per_face_normals.h>
#include <igl/barycentric_coordinates.h>

void morpher::CalculateBaryCoord(){
    baryCoords.clear();
    baryCoords = GenerateBarycentricCoord(Vf,Tf,V);
    assert(baryCoords.size()==V.rows());
}

void morpher::SetTGF(MatrixXd C, MatrixXi BE){
    orientation.clear();origin.clear();
    vector<Vector3d> cVec;
    for(int i=0;i<C.rows();i++)
        cVec.push_back(Vector3d(C(i,0),C(i,1),C(i,2)));
    for(int i=0;i<BE.rows();i++){
        origin.push_back(cVec[BE(i,0)]);
        orientation.push_back((cVec[BE(i,1)]-cVec[BE(i,0)]).normalized());
    }
}

#include <igl/ray_mesh_intersect.h>
//void morpher::StartMorph(string target){
//    cout<<"<Start Morph>"<<endl;cout<<"Reading "+target+"..."<<flush;

//    cout<<"\tread "<<target<<"..."<<flush;
//    MatrixXd Vt;
//    MatrixXi Ft;
//    igl::readPLY(target,Vt,Ft);
//    cout<<"done ("<<Vt.rows()<<"v, "<<Ft.rows()<<"f)"<<endl;

//    cout<<"move vertices..."<<flush;
//    U = V;
//    //for(int i:selected){
////    for(int i=0;i<V1.rows();i++){
////        if(closestBone(i)<0) continue;
////        int bone = closestBone(i);
////        //armpit chk
////        bool armpit(false);
////        if(closestBone(i)==10||closestBone(i)==2){
////            armpit=true;
////            //bone = 7;
////        }

////        Vector3d s(V1(i,0),V1(i,1),V1(i,2));
////        Vector3d p0 = origin[bone]+orientation[bone]*(s-origin[bone]).dot(orientation[bone]);
////        Vector3d dir = s-p0; dir.normalize();
////        igl::Hit hit, hit1;
////        if(!igl::ray_mesh_intersect(s,dir,Vt,Ft,hit)) hit.t=100;
////        dir = -dir;
////        if(!igl::ray_mesh_intersect(s,dir,Vt,Ft,hit1)) hit1.t=100;
////        double dist;
////        if(hit.t<hit1.t) {
////            dist = hit.t; dir = -dir;
////        }else dist = hit1.t;

////        if(limited.find(i)!=limited.end()){
////            if(limited[i]>0) dist*=limited[i];
////            else if(dist>(-limited[i])) dist = 0;
////        }
////        if(armpit && dist>1) dist=0;
////        if(dist>10) dist=0;

////        Vector3d v = s + dir*dist;
////        V2(i,0) = v(0);
////        V2(i,1) = v(1);
////        V2(i,2) = v(2);
////    }
//    for(auto w:W_morph){
//      //  if(closestBone(w.first)<0) continue;
//        Vector3d dir(0,0,0);
//        Vector3d s(V(w.first,0),V(w.first,1),V(w.first,2));
//        for(auto b:w.second){
//            Vector3d p0 = origin[b.first]+orientation[b.first]*(s-origin[b.first]).dot(orientation[b.first]);
//            dir += (s-p0).normalized() * b.second;
//        }

//        igl::Hit hit, hit1;
//        if(!igl::ray_mesh_intersect(s,dir,Vt,Ft,hit)) hit.t=100;
//        dir = -dir;
//        if(!igl::ray_mesh_intersect(s,dir,Vt,Ft,hit1)) hit1.t=100;
//        double dist;
//        if(hit.t<hit1.t) {
//            dist = hit.t; dir = -dir;
//        }else dist = hit1.t;

//        if(limited.find(w.first)!=limited.end()){
//            dist*=limited[w.first];
//        }
//        //if(dist>10) dist=0;

//        Vector3d v = s + dir*dist;
//        U(w.first,0) = v(0);
//        U(w.first,1) = v(1);
//        U(w.first,2) = v(2);
//    }
//    igl::writePLY("out.ply",V2,F1,0);
//    cout<<"done"<<endl;
//}

void morpher::CalculateMorphWeight(MatrixXd &W_b){
//    cout<<"<CalculateClosestBone>"<<endl;

//    cout<<"\textract weights..."<<flush;
//    MatrixXd W_b1;
//    W_b1.resize(V.rows(),W_b.cols());
//    for(size_t i=0;i<wID.size();i++)
//        W_b1.row(i) = W_b.row(wID[i]);
//    cout<<"done"<<endl;

//    cout<<"\torganize weights..."<<flush;
//    W_morph.clear();
//    double epsl(1e-6);
//    int count(0);
//    for(int i=0;i<W_b1.rows();i++){
//        double sum(0);
//        for(int j=0;j<W_b1.cols();j++){
//            if(W_b1(i,j)<epsl) continue;
//            W_morph[i][j]=W_b1(i,j);
//            sum += W_b1(i,j);void Morph(){
//            count++;
//        }
//        for(auto &w:W_morph[i]) w.second/=sum;
//    }
//    cout<<"done ("<<count<<" w for "<<W_morph.size()<<" v)"<<endl;

//    cout<<"\tset extrimities vertices..."<<flush;
//    Eigen::VectorXd maxVal1;
//    Eigen::VectorXi maxIdx1;
//    igl::mat_max(W_b1, 2, maxVal1, maxIdx1);
//    set<int> extrimities = {6,8,9,14,18,22};
//    set<int> fixIdx;
//    for(int i=0;i<maxIdx1.rows();i++)
//        if(extrimities.find(maxIdx1(i))!=extrimities.end()) fixIdx.insert(i);
//    cout<<"done ("<<fixIdx.size()<<")"<<endl;
//    MatchVertices("fixed.ply",fixIdx);
//    cout<<"\tfix vertices..."<<flush;
//    for(auto id:fixIdx) W_morph.erase(id);
//    cout<<"done ("<<fixIdx.size()<<")"<<endl;

//    cout<<"\tset limits..."<<flush;
//    for(int i=0;i<F.rows();i++){
//        vector<int> fixed;
//        vector<int> notFixed;
//        for(int j=0;j<2;j++){
//            if(fixIdx.find(F(i,j))!=fixIdx.end()) fixed.push_back(F1(i,j));
//            else notFixed.push_back(F1(i,j));
//        }
//        if(fixed.size()==0 || fixed.size()==3) continue;
//        for(auto id:fixed) limited[id]=0.3;
//        for(auto id:notFixed) limited[id]=0.7;
//    }
//    count = 0;
//    for(auto iter:limited) if(iter.second<0.5) count++;
//    cout<<"done (0.3: "<<count<<" / 0.5: "<<limited.size()-count<<")"<<endl;
}

map<int, MatrixXd> morpher::GenerateBarycentricCoord(MatrixXd V_f, MatrixXi T_f, MatrixXd V){
    map<tuple<int,int,int>,vector<int>> grid = GenerateGrid(V);
    map<int, MatrixXd> _bary;
    for(int n=0;n<T_f.rows();n++){
        MatrixXd tet(4,3);
        for(int e=0;e<4;e++) tet.row(e) = V_f.row(T_f(n,e));
        MatrixXd max = tet.colwise().maxCoeff();
        MatrixXd min = tet.colwise().minCoeff();
        int i_max = floor(max(0,0)+0.5);int i_min = floor(min(0,0)+0.5);
        int j_max = floor(max(0,1)+0.5);int j_min = floor(min(0,1)+0.5);
        int k_max = floor(max(0,2)+0.5);int k_min = floor(min(0,2)+0.5);
        for(int i=i_min-1;i<i_max+2;i++){
            for(int j=j_min-1;j<j_max+2;j++){
                for(int k=k_min-1;k<k_max+2;k++){
                    auto key = make_tuple(i,j,k);
                    for(int idx:grid[key]){
                        bool done = _bary.find(idx) != _bary.end();
                        if(done && _bary[idx].minCoeff()>0) continue;
                        MatrixXd bary;
                        igl::barycentric_coordinates(V.row(idx),tet.row(0),tet.row(1),tet.row(2),tet.row(3),bary);
                        MatrixXd bary1(1,5); bary1 << bary, n;
                        if(!done) _bary[idx] = bary1;
                        else if(bary.minCoeff() > _bary[idx].minCoeff())
                            _bary[idx] = bary1;
                    }
                }
            }
        }
    }
    cout<<"done ("<<_bary.size()<<"/"<<V.rows()<<")"<<endl;
    return _bary;
}

map<tuple<int,int,int>,vector<int>> morpher::GenerateGrid(MatrixXd V){
    map<tuple<int,int,int>,vector<int>> grid;
    for(int i=0;i<V.rows();i++){
        int x = floor(V(i,0)+0.5);
        int y = floor(V(i,1)+0.5);
        int z = floor(V(i,2)+0.5);
        auto key = make_tuple(x,y,z);
        if(grid.find(key)==grid.end()) grid[key]={};
        grid[key].push_back(i);
    }
    return grid;
}
