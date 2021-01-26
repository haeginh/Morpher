#ifndef MORPHER_HH
#define MORPHER_HH

#include <Eigen/Dense>
#include <igl/writePLY.h>
#include <map>
#include <set>
#include <tuple>
#include <Eigen/Sparse>
#include <igl/writeDMAT.h>
#include <igl/mat_max.h>
#include <dual_quat_cu.hpp>

using namespace  Eigen;
using namespace  std;
using namespace  Tbx;

class morpher
{
public:
    morpher(){}
    void SetModelV(MatrixXd _V){V=_V; U=V;}
    int SetFrame(MatrixXd _V, MatrixXi _T, MatrixXi _F){
        Vs=_V; Us=_V; Ts=_T; Fs=_F;
        outer.clear();
        adjacent.clear();
        for(int i=0;i<Fs.rows();i++){
            for(int e=0;e<3;e++){
                outer.push_back(Fs(i,e));
                if(adjacent.find(Fs(i,e))==adjacent.end()) adjacent[Fs(i,e)]={i};
                else adjacent[Fs(i,e)].push_back(i);
            }
        }
        sort(outer.begin(),outer.end());
        outer.erase(unique(outer.begin(),outer.end()),outer.end());
        return outer.size();
    }
    void CalculateBaryCoord(double _offset=0.8);
    void SetTGF(MatrixXd _C, MatrixXi _BE);
    void SetShellVFix(set<int> _fix){fix=_fix;}
    void SetShellWeights(MatrixXd _Wb, MatrixXd _Wj){W_sb = _Wb; W_sj = _Wj;}
    void SetWholeWeights(MatrixXd _Wb, MatrixXd _Wj){
        W_b = _Wb; W_j = _Wj;
        Eigen::VectorXd maxVal1;
        Eigen::VectorXi maxIdx1;
        igl::mat_max(W_b, 2, maxVal1, maxIdx1);
        set<int> extrimities = {6,8,9,14,18,22};
        for(int i=0;i<maxIdx1.rows();i++)
            if(extrimities.find(maxIdx1(i))!=extrimities.end()) extrmtV[i] = maxVal1(i);
    }
    void MatchJoints(MatrixXd C_target);
    void CalculateMorphWeight();
    map<int, map<int, double>> GetMorphWeights(){return W_morph;}
    void StartMorph(MatrixXd _V, MatrixXi _F);
    void GetMorphedVertices(MatrixXd &_V){
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
        _V = baryMat*Usbary;//+triGapMat;
        for(auto i:extrmtV) _V.row(i.first)= _V.row(i.first)*(1.-i.second) + U.row(i.first)*i.second;
    }
    MatrixXd GetMorphedShellVerts(){return Us;}
    MatrixXd GetWholeVerts(){return U;}
    void PrintPLY(string file){
        //igl::writePLY(file, V1, F1, 0);
    }
private:
    map<int, MatrixXd> GenerateBarycentricCoord(MatrixXd V_f, MatrixXi T_f, MatrixXd V);
    map<tuple<int,int,int>,vector<int>> GenerateGrid(MatrixXd V);
    pair<vector<Vec3>,vector<Transfo>> CalculateScaleRot(MatrixXd C, MatrixXd C_target, MatrixXi BE);

private:
    MatrixXd Vs, Us, V, U, W_b, W_j, W_sb, W_sj, C;
    MatrixXi Ts, BE, Fs;
    set<int> fix;
    SparseMatrix<double> baryMat;
    map<int, pair<int, double>> triGap;
    vector<int> outer;
    //vector<int> outerTet;
    map<int,vector<int>> adjacent;
    map<int, double> extrmtV;
    double offset;
    vector<Vector3d> origin, orientation;
    map<int, double> limited;
    map<int, map<int, double>> W_morph;
};

#endif // MORPHER_HH
