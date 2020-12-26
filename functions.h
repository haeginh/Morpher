#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cctype>
#include <algorithm>
#include "dual_quat_cu.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/bbw.h>
#include <igl/boundary_conditions.h>
#include <igl/readMESH.h>
#include <igl/writeMESH.h>
#include <igl/writePLY.h>
#include <igl/opengl/glfw/Viewer.h>
#include "G4Tet.hh"
using namespace  std;
using namespace  Tbx;
using namespace  Eigen;

bool is_number(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(),
        s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

void ReadTargetOBJ(string target, MatrixXd &V, MatrixXi &F, MatrixXd &C){
    ifstream ifs(target);
    if(!ifs.is_open()){
        cout<<target<<" is not open..Abort program"<<endl;
        exit(1);
    }

    string aLine;
    vector<Vec3> verts, verts_temp, vertsPLY;
    map<int, vector<int>> shellV;
    vector<vector<int>> facePLY;
    double x,y,z;
    int centerID(-1), a, b, c, prevNum;
    map<int, Vec3> centers;
    while(getline(ifs, aLine)){
        stringstream ss(aLine);
        string dump;
        ss>>dump;
        if(dump=="v"){
            ss>>x>>y>>z;
            verts_temp.push_back(Vec3(x,y,z));
        }
        else if(dump=="g"){
            ss>>dump;
            if(!is_number(dump)){
                centerID = -1;
                vertsPLY = verts_temp;
                prevNum = verts.size();
            }
            else {
                centerID=atoi(dump.c_str());
                shellV[centerID]={};
            }
            verts.insert(verts.end(),verts_temp.begin(),verts_temp.end());
            verts_temp.clear();
        }
        else if(dump=="f"){
            ss>>a>>b>>c;
            if(centerID<0){
                facePLY.push_back({a-1-prevNum,b-1-prevNum,c-1-prevNum});
                continue;
            }
            shellV[centerID].push_back(a-1);
            shellV[centerID].push_back(b-1);
            shellV[centerID].push_back(c-1);
        }
    }ifs.close();
    C.resize(shellV.size(), 3);
    for(auto &iter:shellV){
        if(iter.first<0) continue;
        centers[iter.first] = Vec3(0,0,0);
        sort(iter.second.begin(), iter.second.end());
        iter.second.erase(unique(iter.second.begin(), iter.second.end()), iter.second.end());

        for(auto v:iter.second){
            centers[iter.first] += verts[v];
        }
        centers[iter.first] /= (double)iter.second.size();
        C(iter.first-1, 0) = centers[iter.first].x;
        C(iter.first-1, 1) = centers[iter.first].y;
        C(iter.first-1, 2) = centers[iter.first].z;
    }
    //MatrixXd V; MatrixXi F;
    V.resize(vertsPLY.size(),3);
    for(int i=0;i<vertsPLY.size();i++){
        V(i,0) = vertsPLY[i].x;
        V(i,1) = vertsPLY[i].y;
        V(i,2) = vertsPLY[i].z;
    }
    F.resize(facePLY.size(),3);
    for(int i=0;i<facePLY.size();i++){
        F(i,0) = facePLY[i][0];
        F(i,1) = facePLY[i][1];
        F(i,2) = facePLY[i][2];
    }
}

vector<Transfo> CalculateTransfo(MatrixXd C, MatrixXd C_target, MatrixXi BE){
    vector<Transfo> transformations;
    for(int i=0;i<BE.rows();i++){
        //translation to mother joint
        Vec3 model_c  = Vec3(C(BE(i, 0),0),C(BE(i, 0),1),C(BE(i, 0),2));
        Vec3 target_c = Vec3(C_target(BE(i, 0),0),C_target(BE(i, 0),1),C_target(BE(i, 0),2));
        //Transfo trans = target_c - model_c;

        //rotation
        Vec3 model_o  = Vec3(C(BE(i, 1),0),C(BE(i, 1),1),C(BE(i, 1),2)) - model_c;
        Vec3 target_o = Vec3(C_target(BE(i, 1),0),C_target(BE(i, 1),1),C_target(BE(i, 1),2)) - target_c;
        Vec3 axis = model_o.cross(target_o); axis.normalize();
        double angle = model_o.normalized().dot(target_o.normalized());
        Transfo rot;
        if(angle > 0.9999999) rot = Transfo::translate(target_c-model_c);
        else rot = Transfo::translate(target_c)*Transfo::rotate(axis, acos(angle))*Transfo::translate(-model_c);

        //scaling
        double r = (target_o.norm()-model_o.norm()) - 1.;
        Transfo scaling = target_o.normalized() * r;

        transformations.push_back(scaling*rot);
    }
    return transformations;
}

pair<vector<Vec3>,vector<Transfo>> CalculateScaleRot(MatrixXd C, MatrixXd C_target, MatrixXi BE){
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

map<tuple<int,int,int>,vector<int>> GenerateGrid(MatrixXd V, double size=1.){
    map<tuple<int,int,int>,vector<int>> grid;
    for(int i=0;i<V.rows();i++){
        int x = floor(V(i,0)*size+0.5);
        int y = floor(V(i,1)*size+0.5);
        int z = floor(V(i,2)*size+0.5);
        auto key = make_tuple(x,y,z);
        if(grid.find(key)==grid.end()) grid[key]={};
        grid[key].push_back(i);
    }
    return grid;
}

void CalculateScalingWeights(MatrixXd &C, MatrixXd &V,MatrixXi &T,MatrixXd &W){
    cout<<endl<<"<CalculateScailingWeights>"<<endl;
    VectorXi b; b.resize(C.rows());
    MatrixXd bc = MatrixXd::Zero(C.rows(),C.rows());
    map<tuple<int,int,int>,vector<int>> grid = GenerateGrid(V);
    for(int n=0;n<C.rows();n++){
        int x = floor(C(n,0)+0.5);
        int y = floor(C(n,1)+0.5);
        int z = floor(C(n,2)+0.5);
        auto key = make_tuple(x,y,z);
        for(int i:grid[key]){
            if(fabs(C(n,0)-V(i,0))>0.01) continue;
            if(fabs(C(n,1)-V(i,1))>0.01) continue;
            if(fabs(C(n,2)-V(i,2))>0.01) continue;
            b(n) = i;
            bc(n,n) = 1;
            break;
        }
    }
    // compute BBW weights matrix
    igl::BBWData bbw_data;
    // only a few iterations for sake of demo
    bbw_data.active_set_params.max_iter = 15;
    bbw_data.verbosity = 2;
    igl::bbw(V,T,b,bc,bbw_data,W);
}

void EraseInnerTri(string meshN){
    ifstream ifs(meshN);
    ofstream ofs(meshN.substr(0, meshN.size()-4)+"2.mesh");
    int numTri; string dump, str; int a, b, c, id;

    while(getline(ifs, dump)){
        stringstream ss(dump);
        ofs<<dump<<endl;
        ss>>str;
        if(str=="Triangles") break;
    }
    ifs>>numTri;
    vector<vector<int>> triVec;
    for(int i=0;i<numTri;i++){
        ifs>>a>>b>>c>>id;
        if(id==-1) triVec.push_back({b,a,c});
    }
    cout<<"Surface Triangles: "<<triVec.size()<<"/"<<numTri<<endl;
    ofs<<triVec.size()<<endl;
    for(auto tri:triVec){
        ofs<<tri[0]<<" "<<tri[1]<<" "<<tri[2]<<" 0"<<endl;
    }

    while(getline(ifs, dump)){
        ofs<<dump<<endl;
    }
    ifs.close();
    ofs.close();
}

void GenerateMorphingWeights(string fileN, string inputMesh, MatrixXd &C, MatrixXd &V,MatrixXi &T,MatrixXi &BE, MatrixXi &F){
    //read mesh
    MatrixXd V1; MatrixXi T1, F1;
    cout<<"read "+inputMesh<<"..."<<flush;
    igl::readMESH(inputMesh,V1,T1,F1);
    ifstream ifs(inputMesh); string dump; int vNum;
    while(ifs>>dump){
        if(dump=="Vertices") {
            ifs>>vNum;
            break;
        }
    }
    double x,y,z; int idx;
    vector<int> inner;
    for(int i=0;i<vNum;i++){
        ifs>>x>>y>>z>>idx;
        if(idx!=0) inner.push_back(i);
    }ifs.close();
    cout<<"done (v "<<V1.rows()<<", i "<<inner.size()<<")"<<endl;

    //combine two meshes
    cout<<"comparing vertices..."<<flush;
    map<tuple<int, int, int>, int> v10;
    map<int, int> i2w;
    for(int i=0;i<V.rows();i++)
        v10[make_tuple(floor(V(i,0)*10+0.5),floor(V(i,1)*10+0.5),floor(V(i,2)*10+0.5))] = i;
    for(int i=0;i<V1.rows();i++){
        if(find(inner.begin(),inner.end(),i)!=inner.end()) continue;
        auto key = make_tuple(floor(V1(i,0)*10+0.5),floor(V1(i,1)*10+0.5),floor(V1(i,2)*10+0.5));
        auto iter = v10.find(key);
        if(iter == v10.end()) continue;
        i2w[i] = iter->second;
    }cout<<i2w.size()<<" v found"<<endl;
    MatrixXd V2; V2.resize(V.rows()+V1.rows()-i2w.size(),3);
    V2.block(0,0,V.rows(),3) = V;
    int count(0);
    for(int i=0;i<V1.rows();i++){
        if(i2w.find(i)!=i2w.end()) continue;
        i2w[i]=V.rows()+count;
        V2.row(V.rows()+count) = V1.row(i);
        count++;
    }
    MatrixXi T2; T2.resize(T.rows()+T1.rows(),4);
    T2.block(0,0,T.rows(),4) = T;
    for(int i=0;i<T1.rows();i++){
        T2(T.rows()+i,0) = i2w[T1(i,0)];
        T2(T.rows()+i,1) = i2w[T1(i,1)];
        T2(T.rows()+i,2) = i2w[T1(i,2)];
        T2(T.rows()+i,3) = i2w[T1(i,3)];
    }
    MatrixXi F2; F2.resize(F.rows()+F1.rows(),3);
    F2.block(0,0,F.rows(),3) = F;
    for(int i=0;i<F1.rows();i++){
        F2(F.rows()+i,0) = i2w[F1(i,0)];
        F2(F.rows()+i,1) = i2w[F1(i,1)];
        F2(F.rows()+i,2) = i2w[F1(i,2)];
    }
//    igl::writePLY("merged.ply", V2, F2, 0);
//    igl::writeMESH("merged.mesh", V2, T2, F2);
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V2, F2);
    viewer.launch();
}

MatrixXd GenerateBonePoints(MatrixXd C, MatrixXi BE, double interval){
    MatrixXd output=C;
    for(int i=0;i<BE.rows();i++){
        Vector3d c0 = C.row(BE(i,0)).transpose();
        Vector3d c1 = C.row(BE(i,1)).transpose();
        Vector3d dir = (c1-c0).normalized();
        double l = (c1-c0).norm();
        for(double i=interval;i<l;i+=interval){
            output.conservativeResize(output.rows()+1, NoChange);
            output.row(output.rows()-1) = (c0+dir*i).transpose();
        }
    }
    return output;
}

void readMESH(string fileName,MatrixXd &V,MatrixXi &T,MatrixXi &F){
    ifstream ifs(fileName);
    if(!ifs.is_open())
        cerr<<fileName<<" is not open!!"<<endl;
    string dump;
    int num;
    while(getline(ifs,dump)){
        stringstream ss(dump);
        ss>>dump;
        if(dump=="Vertices"){
            ifs>>num;
            break;
        }
    }
    V.resize(num, 3);
    double x,y,z;
    int chk;
    for(int i=0;i<num;i++){
        ifs>>x>>y>>z>>chk;
        V(i,0)=x;V(i,1)=y;V(i,2)=z;
    }

    while(getline(ifs,dump)){
        stringstream ss(dump);
        ss>>dump;
        if(dump=="Triangles"){
            ifs>>num;
            break;
        }
    }
    int a, b, c, d;
    vector<vector<int>> faces;
    for(int i=0;i<num;i++){
        ifs>>a>>b>>c>>chk;
//        cout<<a<<" "<<b<<" "<<c<<" "<<chk<<flush;
//        getchar();
        if(chk==-1) faces.push_back({b-1,a-1,c-1});
    }
    F.resize(faces.size(),3);
    for(int i=0;i<faces.size();i++){
        for(int j=0;j<3;j++)
            F(i,j)=faces[i][j];
    }

    while(getline(ifs,dump)){
        stringstream ss(dump);
        ss>>dump;
        if(dump=="Tetrahedra"){
            ifs>>num;
            break;
        }
    }
    T.resize(num,4);
    for(int i=0;i<num;i++){
        ifs>>a>>b>>c>>d>>chk;
//        cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<chk<<flush;
//        getchar();
        T(i,0)=a-1;T(i,1)=b-1;T(i,2)=c-1;T(i,3)=d-1;
    }ifs.close();
}

void CleanWeights(MatrixXd &W){
    double epsilon (1e-6);
    for(int i=0;i<W.rows();i++){
        double sum(0);
        for(int j=0;j<W.cols();j++){
            if(W(i,j)<epsilon) W(i,j)=0;
            else sum += W(i,j);
        }
        for(int j=0;j<W.cols();j++) W(i,j)/=sum;
    }
}

//#include "igl/biharmonic_coordinates.h"
//#include "igl/point_mesh_squared_distance.h"
#include "igl/point_simplex_squared_distance.h"
#include "igl/barycentric_coordinates.h"

map<int, MatrixXd> GenerateBarycentricCoord(MatrixXd V_f, MatrixXi T_f, MatrixXd V){
    map<tuple<int,int,int>,vector<int>> grid = GenerateGrid(V);
    map<int, MatrixXd> baryCoords;
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
                        bool done = baryCoords.find(idx) != baryCoords.end();
                        if(done && baryCoords[idx].minCoeff()>0) continue;
                        MatrixXd bary;
                        igl::barycentric_coordinates(V.row(idx),tet.row(0),tet.row(1),tet.row(2),tet.row(3),bary);
                        MatrixXd bary1(1,5); bary1 << bary, n;
                        if(!done) baryCoords[idx] = bary1;
                        else if(bary.minCoeff() > baryCoords[idx].minCoeff())
                            baryCoords[idx] = bary1;
                    }
                }
            }
        }
    }
    cout<<"("<<baryCoords.size()<<"/"<<V.rows()<<")"<<endl;
    return baryCoords;
}
void PrintBaryCoords(string file, map<int, MatrixXd> &baryCoords){
    ofstream ofs(file); ofs<<baryCoords.size()<<endl;
    for(int i=0;i<baryCoords.size();i++){
        ofs<<baryCoords[i](0,4)<<"\t"
           <<baryCoords[i](0,0)<<"\t"
           <<baryCoords[i](0,1)<<"\t"
           <<baryCoords[i](0,2)<<"\t"
           <<baryCoords[i](0,3)<<endl;
    }ofs.close();
}
map<int, MatrixXd> ReadBaryFile(string file){
    map<int, MatrixXd> baryCoords;
    ifstream ifs(file);
    if(!ifs.is_open()) {
        cout<<file+" is not open!!!"<<endl;
        return baryCoords;
    }
    int num, tetID; ifs>>num;
    double a,b,c,d;
    for(int i=0;i<num;i++){
        ifs>>tetID>>a>>b>>c>>d;
        MatrixXd bary(1,5);
        bary<<a,b,c,d,tetID;
        baryCoords[i] = bary;
    }ifs.close();
    return baryCoords;
}
SparseMatrix<double> GenerateBarySparse(map<int, MatrixXd> &baryCoords, int v_size, MatrixXi tet){
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    for(auto iter:baryCoords){
        for(int i=0;i<4;i++)
            triplets.push_back(T(iter.first,tet(iter.second(0,4),i),iter.second(0,i)));
    }
    SparseMatrix<double> mat(baryCoords.size(),v_size);
    mat.setFromTriplets(triplets.begin(),triplets.end());
    return mat;
}
string readOBJ(string fileN,MatrixXd &V){
    ifstream ifs(fileN);
    if(!ifs.is_open())
        cout<<fileN<<" is not open!!"<<endl;
    stringstream ssOut;
    vector<Vector3d> vertices;
    while(!ifs.eof()){
        string dump;
        getline(ifs,dump);
        stringstream ss(dump);
        string first;
        ss>>first;
        double x,y,z;
        if(first=="v"){
            ss>>x>>y>>z;
            vertices.push_back(Vector3d(x,y,z));
        }
        else if(first=="mtllib") continue;
        else ssOut<<dump<<endl;
    }ifs.close();
    V.resize(vertices.size(),3);
    for(int i=0;i<vertices.size();i++){
        V.row(i)=vertices[i].transpose();
    }
    return ssOut.str();
}

map<int,int> CompareVertices(MatrixXd V1,MatrixXd V2, double epsil = 0.0001){
    map<tuple<int,int,int>,vector<int>> grid = GenerateGrid(V1);
    map<int,int> conv;
    for(int i=0;i<V2.rows();i++){
        int x = floor(V2(i,0)+0.5);
        int y = floor(V2(i,1)+0.5);
        int z = floor(V2(i,2)+0.5);
        auto key = make_tuple(x,y,z);
        bool chk(false);
        for(int idx:grid[key]){
            if(fabs(V1(idx,0)-V2(i,0))>epsil) continue;
            if(fabs(V1(idx,1)-V2(i,1))>epsil) continue;
            if(fabs(V1(idx,2)-V2(i,2))>epsil) continue;
            conv[idx] = i;
        }
    }
    return conv;
}

MatrixXd ReadShellOBJ(string shellOBJ, MatrixXd &V_s, MatrixXi &F_s){
    ifstream ifs(shellOBJ);
    if(!ifs.is_open())
        cerr<<shellOBJ + " is not open!!"<<endl;

    cout<<"read "+shellOBJ+"..."<<flush;
    int idx(0), prev_idx(0),a,b,c;
    double x,y,z;
    string shellN;
    vector<Vector3d> verts_temp;
    map<string,vector<Vector3d>> verts;
    map<string,vector<Vector3i>> faces;
    while (!ifs.eof()) {
        string dump, first;
        getline(ifs,dump);
        stringstream ss(dump);
        ss>>first;
        if(first=="v"){
            ss>>x>>y>>z;
            verts_temp.push_back(Vector3d(x,y,z));
        }else if(first=="g"){
            idx += prev_idx;
            ss>>dump; shellN = dump;
            verts[shellN] = verts_temp;
            faces[shellN] = {};
            prev_idx = verts_temp.size();
            verts_temp.clear();
        }else if(first=="f"){
            ss>>a>>b>>c;
            faces[shellN].push_back(Vector3i(a-1-idx,b-1-idx,c-1-idx));
        }
    }ifs.close();
    cout<<"done"<<endl;

    V_s.resize(verts["1"].size(),3);
    F_s.resize(faces["1"].size(),3);
    for(int i=0;i<verts["1"].size();i++)
        V_s.row(i) = verts["1"][i].transpose();
    for(int i=0;i<faces["1"].size();i++)
        F_s.row(i) = faces["1"][i].transpose();

    MatrixXd V_fix;
    V_fix.resize(verts["-1"].size(),3);
    for(int i=0;i<verts["-1"].size();i++)
        V_fix.row(i) = verts["-1"][i].transpose();

    return V_fix;
}

#endif // FUNCTIONS_H
