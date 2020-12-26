#ifndef MORPHER_HH
#define MORPHER_HH

#include <Eigen/Dense>
#include <igl/writePLY.h>
#include <map>
#include <set>
#include <tuple>
#include <Eigen/Sparse>
#include <igl/writeDMAT.h>

using namespace  Eigen;
using namespace  std;

class morpher
{
public:
    morpher(){}
    void SetModelV(MatrixXd _V){V=_V;}
    void SetFrame(MatrixXd _V, MatrixXi _T){Vf=_V; Tf=_T;}
    void CalculateBaryCoord();
    void PrintBaryCoord(string file){
        ofstream ofs(file); ofs<<baryCoords.size()<<endl;
        for(int i=0;i<baryCoords.size();i++){
            ofs<<baryCoords[i](0,4)<<"\t"
               <<baryCoords[i](0,0)<<"\t"
               <<baryCoords[i](0,1)<<"\t"
               <<baryCoords[i](0,2)<<"\t"
               <<baryCoords[i](0,3)<<endl;
        }ofs.close();
    }
    void SetBaryCoord(map<int, MatrixXd> _bary) {baryCoords=_bary;}
    void SetTGF(MatrixXd C, MatrixXi BE);
    void SetShellVFix(set<int> _fix){fix=_fix;}
    void CalculateMorphWeight(MatrixXd &W_b);
    void StartMorph(MatrixXd _V, MatrixXi _F){}
    void GetMorphedVertices(MatrixXd &V){}//V=V2;}
    void PrintPLY(string file){
        //igl::writePLY(file, V1, F1, 0);
    }
private:
    map<int, MatrixXd> GenerateBarycentricCoord(MatrixXd V_f, MatrixXi T_f, MatrixXd V);
    map<tuple<int,int,int>,vector<int>> GenerateGrid(MatrixXd V);

private:
    MatrixXd Vf, Uf, V, U;
    MatrixXi Tf;
    set<int> fix;
    map<tuple<int,int,int>, vector<int>> grider;
    map<int, MatrixXd> baryCoords;
    vector<int> wID; //morph to whole
    vector<int> selected;
    vector<Vector3d> origin, orientation;
    map<int, double> limited;
    map<int, map<int, double>> W_morph;
};

#endif // MORPHER_HH
