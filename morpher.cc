#include "morpher.hh"
#include <iostream>
#include <igl/readOBJ.h>
#include <igl/dot.h>
#include <igl/per_face_normals.h>
#include <igl/barycentric_coordinates.h>

void morpher::CalculateBaryCoord(double _offset){
    offset = _offset;
    map<int, map<int, double>> baryCoords;
    map<tuple<int,int,int>,vector<int>> grid = GenerateGrid(U);
    //double epsl(1e-5);
    MatrixXd Usbary=Us;
    for(auto i:outer){
        Vector3d dir(0,0,0);
        for(auto f:adjacent[i]){
            vector<Vector3d> tri;
            for(int e=0;e<3;e++) tri.push_back(Us.row(Fs(f,e)).transpose());
            Vector3d normal = (tri[1]-tri[0]).cross(tri[2]-tri[0]);
            dir += normal;
        }
        //if(fix.find(i)!=fix.end())
            Usbary.row(i) = Us.row(i) + offset*dir.normalized().transpose();
        //else
           // Usbary.row(i) = Us.row(i) + 0.02*dir.normalized().transpose();
    }
    for(int n=0;n<Ts.rows();n++){
        vector<Vector3d> tet;
        Vector3d max(-1e10,-1e10,-1e10), min(1e10,1e10,1e10);
        for(int e=0;e<4;e++){
            tet.push_back(Usbary.row(Ts(n,e)).transpose());
            max(0) = max(0)>tet[e](0)? max(0):tet[e](0);
            max(1) = max(1)>tet[e](1)? max(1):tet[e](1);
            max(2) = max(2)>tet[e](2)? max(2):tet[e](2);
            min(0) = min(0)<tet[e](0)? min(0):tet[e](0);
            min(1) = min(1)<tet[e](1)? min(1):tet[e](1);
            min(2) = min(2)<tet[e](2)? min(2):tet[e](2);
        }
        double invVol6 = 1./(tet[1]-tet[0]).cross(tet[2]-tet[0]).dot(tet[3]-tet[0]);
        int i_max = floor(max(0)+0.5);int i_min = floor(min(0)+0.5);
        int j_max = floor(max(1)+0.5);int j_min = floor(min(1)+0.5);
        int k_max = floor(max(2)+0.5);int k_min = floor(min(2)+0.5);
        for(int i=i_min;i<i_max+1;i++){
            for(int j=j_min;j<j_max+1;j++){
                for(int k=k_min;k<k_max+1;k++){
                    auto key = make_tuple(i,j,k);
                    for(int idx:grid[key]){
                        if(baryCoords.find(idx)!=baryCoords.end()) continue;
                        Vector3d v = U.row(idx).transpose();
                        double b0 = (tet[1]-v).cross(tet[2]-v).dot(tet[3]-v)*invVol6;
                        if(b0<0) continue;
                        double b1 = (v-tet[0]).cross(tet[2]-tet[0]).dot(tet[3]-tet[0])*invVol6;
                        if(b1<0) continue;
                        double b2 = (tet[1]-tet[0]).cross(v-tet[0]).dot(tet[3]-tet[0])*invVol6;
                        if(b2<0) continue;
                        double b3 = (tet[1]-tet[0]).cross(tet[2]-tet[0]).dot(v-tet[0])*invVol6;
                        if(b3<0) continue;
                        map<int, double> bary;
                        bary[Ts(n,0)] = b0; bary[Ts(n,1)] = b1; bary[Ts(n,2)] = b2; bary[Ts(n,3)] = b3;
                        baryCoords[idx] = bary;
                    }
                }
            }
        }
        cout<<"\rGenerating barycentric coord..."<<baryCoords.size()<<"/"<<V.rows()<<flush;
    }
     map<int, map<int,double>> baryCoordsOut;

    map<int, pair<int, double>> triGap;
    map<int, double> triDist;
    map<int, map<int,double>> baryCoordsTri;
    for(int n=0;n<Fs.rows();n++){
        vector<Vector3d> tri, tet;
        Vector3d max(-1e10,-1e10,-1e10), min(1e10,1e10,1e10);
        for(int e=0;e<3;e++){
            tri.push_back(Usbary.row(Fs(n,e)).transpose());
            max(0) = max(0)>tri[e](0)? max(0):tri[e](0);
            max(1) = max(1)>tri[e](1)? max(1):tri[e](1);
            max(2) = max(2)>tri[e](2)? max(2):tri[e](2);
            min(0) = min(0)<tri[e](0)? min(0):tri[e](0);
            min(1) = min(1)<tri[e](1)? min(1):tri[e](1);
            min(2) = min(2)<tri[e](2)? min(2):tri[e](2);
        }

        Vector3d normal = (tri[1]-tri[0]).cross(tri[2]-tri[0]).normalized();
        int i_max = floor(max(0)+0.5);int i_min = floor(min(0)+0.5);
        int j_max = floor(max(1)+0.5);int j_min = floor(min(1)+0.5);
        int k_max = floor(max(2)+0.5);int k_min = floor(min(2)+0.5);
        for(int i=i_min-1;i<i_max+2;i++){
            for(int j=j_min-1;j<j_max+2;j++){
                for(int k=k_min-1;k<k_max+2;k++){
                    auto key = make_tuple(i,j,k);
                    for(int idx:grid[key]){
                        if(baryCoords.find(idx)!=baryCoords.end()) continue;
                        Vector3d v = U.row(idx).transpose();
                        double invVol6 = 1./(tri[2]-v).cross(tri[1]-v).dot(tri[0]-v);
                        if(invVol6<0) continue;
                        double dist=normal.dot(v-tri[0]);
                        Vector3d v_proj = v-normal*dist;
                        double b0 = (tri[2]-v).cross(tri[1]-v).dot(v_proj-v)*invVol6; //if(b0<-0.3) continue;
                        double b1 = (tri[2]-v).cross(v_proj-v).dot(tri[0]-v)*invVol6; //if(b1<-0.3) continue;
                        double b2 = (v_proj-v).cross(tri[1]-v).dot(tri[0]-v)*invVol6; //if(b2<-0.3) continue;

                        double minDist;
                        if(b0<0 && b1<0) minDist = (v-tri[2]).norm();
                        else if(b2<0 && b1<0) minDist = (v-tri[0]).norm();
                        else if(b0<0 && b2<0) minDist = (v-tri[1]).norm();
                        else if(b0<0) minDist = (v-tri[1]).cross(tri[2]-tri[1]).norm()/(tri[2]-tri[1]).norm();
                        else if(b1<0) minDist = (v-tri[0]).cross(tri[2]-tri[0]).norm()/(tri[2]-tri[0]).norm();
                        else if(b2<0) minDist = (v-tri[1]).cross(tri[0]-tri[1]).norm()/(tri[0]-tri[1]).norm();
                        else minDist =dist;

                        if(baryCoordsOut.find(idx)!=baryCoordsOut.end()){
                            if(minDist>triDist[idx]) continue;
                        }
                        map<int, double> bary;
                        bary[Fs(n,0)] = b0; bary[Fs(n,1)] = b1; bary[Fs(n,2)] = b2;
                        baryCoordsOut[idx] = bary;
                        triDist[idx] = minDist;
                    }
                }
            }
        }
        cout<<"\rGenerating barycentric coord..."<<baryCoords.size()+baryCoordsOut.size()<<"/"<<V.rows()<<flush;
    }
    for(int i=0;i<V.rows();i++)
        if(baryCoords.find(i)==baryCoords.end()) baryCoords[i]=baryCoordsOut[i];

    if(baryCoords.size()!=V.rows()){
        cout<<"Check if all the vertices are in frame model!!"<<endl; exit(100);
    }
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    for(int i=0;i<baryCoords.size();i++){
        for(auto w:baryCoords[i])
            triplets.push_back(T(i,w.first,w.second));
    }
    SparseMatrix<double> mat(baryCoords.size(),Vs.rows());
    mat.setFromTriplets(triplets.begin(),triplets.end());
    baryMat = mat;

    cout<<endl<<"Generated barycentric coord...(tet "<<baryCoords.size()-baryCoordsOut.size()<<"/ tri "<<baryCoordsOut.size()<<")"<<endl;
}

void morpher::SetTGF(MatrixXd _C, MatrixXi _BE){
    C = _C; BE = _BE;
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
void morpher::StartMorph(MatrixXd Vt, MatrixXi Ft){
    cout<<"<Start Morph>"<<endl;
    MatrixXd Us_0 = Us;
    set<int> outer_set(outer.begin(),outer.end());
    int count(0);
    for(auto i:outer_set){
        if(W_morph.find(i)==W_morph.end()) continue;
        Vector3d s(Us(i,0),Us(i,1),Us(i,2));
        Vector3d dir(0,0,0);
        for(auto b:W_morph[i]){
            Vector3d p0 = origin[b.first]+orientation[b.first]*(s-origin[b.first]).dot(orientation[b.first]);
            dir += (s-p0).normalized() * b.second;
        }
        igl::Hit hit, hit1;
        if(!igl::ray_mesh_intersect(s,dir,Vt,Ft,hit)) hit.t=100;
        dir = -dir;
        if(!igl::ray_mesh_intersect(s,dir,Vt,Ft,hit1)) hit1.t=100;
        double dist;
        if(hit.t<hit1.t) {
            dist = hit.t; dir = -dir;
        }else dist = hit1.t;

        if(limited.find(i)!=limited.end()){
            dist*=limited[i];
        }
        //if(dist>10) dist=0;

        Vector3d v = s + dir*dist;
        Us(i,0) = v(0);
        Us(i,1) = v(1);
        Us(i,2) = v(2);
        cout<<"\rmove vertices..."<<++count<<"/"<<W_morph.size()<<flush;
    }
    for(auto w:W_morph){
        if(outer_set.find(w.first)!=outer_set.end()) continue;
        Vector3d s(Us(w.first,0),Us(w.first,1),Us(w.first,2));
        Vector3d move(0,0,0);
        for(auto b:w.second){
            Vector3d p0 = origin[b.first]+orientation[b.first]*(s-origin[b.first]).dot(orientation[b.first]);
            Vector3d dir = (s-p0).normalized();
            igl::Hit hit0;
            igl::ray_mesh_intersect(s,dir,Us_0,Fs,hit0);
            Vector3d s1=s+hit0.t*dir;
            igl::Hit hit, hit1;
            if(!igl::ray_mesh_intersect(s1,dir,Us,Fs,hit)) hit.t=100;
            dir = -dir;
            if(!igl::ray_mesh_intersect(s1,dir,Us,Fs,hit1)) hit1.t=100;
            double dist;
            if(hit.t<hit1.t) {
                dist = hit.t; dir = -dir;
            }else dist = hit1.t;
            dist *= b.second*(s-p0).norm()/(s1-p0).norm();
            move += dir*dist;
        }
//                    if(limited.find(w.first)!=limited.end()){
//                        move*=limited[w.first];
//                    }
        Vector3d v = s + move;
        Us(w.first,0) = v(0);
        Us(w.first,1) = v(1);
        Us(w.first,2) = v(2);
        // }
        cout<<"\rmove vertices..."<<++count<<"/"<<W_morph.size()<<flush;
    }
    //igl::writePLY("out.ply",V2,F1,0);
    cout<<"\rmove vertices...done - "<<W_morph.size()<<endl;
}

void morpher::CalculateMorphWeight(){
    cout<<"<CalculateMorphWeight>"<<endl;
    cout<<"\torganize weights..."<<flush;
    W_morph.clear();
    double epsl(1e-6);
    int count(0);
    for(int i=0;i<W_sb.rows();i++){
        double sum(0);
        for(int j=0;j<W_sb.cols();j++){
            if(W_sb(i,j)<epsl) continue;
            if(j==15 || j==19) continue;//ignore horizontal skeleton nearby pelvis
            W_morph[i][j]=W_sb(i,j);
            sum += W_sb(i,j);
            count++;
        }
        for(auto &w:W_morph[i]) w.second/=sum;
    }
    cout<<"done ("<<count<<" w for "<<W_morph.size()<<" v)"<<endl;

    cout<<"\tset extrimities vertices..."<<flush;
    Eigen::VectorXd maxVal1;
    Eigen::VectorXi maxIdx1;
    igl::mat_max(W_sb, 2, maxVal1, maxIdx1);
    set<int> extrimities = {6,8,9,14,18,22};
    //set<int> fixIdx;
    for(int i=0;i<maxIdx1.rows();i++)
        if(extrimities.find(maxIdx1(i))!=extrimities.end()) fix.insert(i);
    cout<<"done ("<<fix.size()<<")"<<endl;
    cout<<"\tfix vertices..."<<flush;
    for(auto id:fix) W_morph.erase(id);
    cout<<"done ("<<fix.size()<<")"<<endl;

    cout<<"\tset limits..."<<flush;
    for(int i=0;i<Fs.rows();i++){
        vector<int> fixed;
        vector<int> notFixed;
        for(int j=0;j<2;j++){
            if(fix.find(Fs(i,j))!=fix.end()) fixed.push_back(Fs(i,j));
            else notFixed.push_back(Fs(i,j));
        }
        if(fixed.size()==0 || fixed.size()==3) continue;
        for(auto id:fixed) limited[id]=0.3;
        for(auto id:notFixed) limited[id]=0.7;
    }

    set<int> outer_set(outer.begin(),outer.end());
    double epsilon(1e-6);
    vector<int> boneV;
    for(auto w:W_morph){
        if(outer_set.find(w.first)!=outer_set.end()) continue;
        Vector3d s(Vs(w.first,0),Vs(w.first,1),Vs(w.first,2));
        int idMax; double wMax(-1);
        for(auto b:w.second){
            if(wMax<b.second){wMax = b.second; idMax = b.first;}
        }
        Vector3d p0 = origin[idMax]+orientation[idMax]*(s-origin[idMax]).dot(orientation[idMax]);
        if((s-p0).squaredNorm()<epsilon) boneV.push_back(w.first);
    }
    for(auto i:boneV) W_morph.erase(i);

    count = 0;
    for(auto iter:limited) if(iter.second==0.7) count++;
    cout<<"done (0.3: "<<limited.size()-count<<" / 0.7: "<<count<<" / b "<<boneV.size()<<")"<<endl;
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

void morpher::MatchJoints(MatrixXd C_target){
    //morphing - use different weights for scaling and rotation
    pair<vector<Vec3>,vector<Transfo>> scaleRot = CalculateScaleRot(C, C_target, BE);
    Us = Vs;
    for(int i=0;i<Vs.rows();i++){
        Vec3 scale_v = scaleRot.first[0]*W_sj(i, 0);
        for(int j=1;j<C.rows();j++)
            scale_v = scale_v + scaleRot.first[j]*W_sj(i, j);

        Transfo rot_tf = scaleRot.second[0]*W_sb(i, 0);
        for(int j=1;j<BE.rows();j++)
            rot_tf = rot_tf + scaleRot.second[j]*W_sb(i, j);

        Point3 p1 = rot_tf*(Point3(Us(i,0),Us(i,1),Us(i,2))+scale_v);
        Us(i, 0) = p1.x; Us(i, 1) = p1.y; Us(i, 2) = p1.z;
    }

    U = V;
    for(int i=0;i<V.rows();i++){
        Vec3 scale_v = scaleRot.first[0]*W_j(i, 0);
        for(int j=1;j<C.rows();j++)
            scale_v = scale_v + scaleRot.first[j]*W_j(i, j);

        Transfo rot_tf = scaleRot.second[0]*W_b(i, 0);
        for(int j=1;j<BE.rows();j++)
            rot_tf = rot_tf + scaleRot.second[j]*W_b(i, j);

        Point3 p1 = rot_tf*(Point3(U(i,0),U(i,1),U(i,2))+scale_v);
        U(i, 0) = p1.x; U(i, 1) = p1.y; U(i, 2) = p1.z;
    }
    SetTGF(C_target,BE);
}

pair<vector<Vec3>,vector<Transfo>> morpher::CalculateScaleRot(MatrixXd C, MatrixXd C_target, MatrixXi BE){
    vector<Vec3> scales; scales.resize(C.rows());
    vector<Vec3> newCen; newCen.resize(C.rows());
    newCen[0] = Vec3(C_target(0,0),C_target(0,1),C_target(0,2));

    for(int i=0;i<BE.rows();i++){
        Vec3 model_c  = Vec3(C(BE(i, 0),0),C(BE(i, 0),1),C(BE(i, 0),2));
        Vec3 target_c = Vec3(C_target(BE(i, 0),0),C_target(BE(i, 0),1),C_target(BE(i, 0),2));
        Vec3 model_o  = Vec3(C(BE(i, 1),0),C(BE(i, 1),1),C(BE(i, 1),2)) - model_c;
        Vec3 target_o = Vec3(C_target(BE(i, 1),0),C_target(BE(i, 1),1),C_target(BE(i, 1),2)) - target_c;
        newCen[BE(i,1)] = newCen[BE(i,0)] + model_o.normalized()*target_o.norm();
    }

    for(int i=0;i<C.rows();i++){
        Vec3 model_c  = Vec3(C(i,0),C(i,1),C(i,2));
        scales[i] = newCen[i] - model_c;
    }

    //rotation
    vector<Transfo> rotVec;
    for(int i=0;i<BE.rows();i++){
        //translation to mother joint
        Vec3 model_c  = Vec3(C(BE(i, 0),0),C(BE(i, 0),1),C(BE(i, 0),2));
        Vec3 target_c = Vec3(C_target(BE(i, 0),0),C_target(BE(i, 0),1),C_target(BE(i, 0),2));

        //rotation
        Vec3 model_o  = Vec3(C(BE(i, 1),0),C(BE(i, 1),1),C(BE(i, 1),2)) - model_c;
        Vec3 target_o = Vec3(C_target(BE(i, 1),0),C_target(BE(i, 1),1),C_target(BE(i, 1),2)) - target_c;
        Vec3 axis = model_o.cross(target_o); axis.normalize();
        double angle = model_o.normalized().dot(target_o.normalized());
        Transfo rot;
        if(angle > 0.9999999) rot = Transfo::translate(target_c-newCen[BE(i,0)]);
        else rot = Transfo::translate(target_c)*Transfo::rotate(axis, acos(angle))*Transfo::translate(-newCen[BE(i,0)]);

        rotVec.push_back(rot);
    }
    return make_pair(scales,rotVec);
}


