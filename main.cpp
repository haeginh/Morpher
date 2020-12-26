#include <igl/readTGF.h>
#include <igl/writeTGF.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/readPLY.h>
#include <igl/directed_edge_parents.h>
#include <vector>
#include <algorithm>
#include <iostream>

#include "functions.h"
#include "morpher.hh"
#define PI       3.14159265358979323846

using namespace Eigen;
using namespace std;
const RowVector3d sea_green(70./255.,252./255.,167./255.);
const RowVector3d red(1.,0.,0.);
const RowVector3d blue(0.,0.,1.);
MatrixXd V,U,C, C_widen;
MatrixXd V_target,C_target;
MatrixXd W, W_s, W_shell, W_s_shell, W_dp;
MatrixXd V_shell;
MatrixXi T,F,BE,F_target, F_shell;
Eigen::VectorXi P;
map<int, int> m2p, o2m;
int selected = 0;
bool toggle(false);
double degree_arm(0),degree_leg(0);
void PrintPly(string);
void PrintOBJ(string fileN, string mtlN);
string objStr;
void Morph();
morpher* morph;
void LimbWiden_arm(double degree, MatrixXd &U, MatrixXd W);
void LimbWiden_arm(double degree, MatrixXd &C);
void LimbWiden_leg(double degree, MatrixXd &U, MatrixXd W);
void LimbWiden_leg(double degree, MatrixXd &C);
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case '0':
      viewer.data().clear();
      viewer.data().set_mesh(V,F);
      //morph->UpdateFrame(V);
      viewer.data().set_edges(C,BE,sea_green);
      break;
  case '1':
      toggle = !toggle;
    viewer.data().clear();
    if(toggle)  {
        viewer.data().set_mesh(V_shell,F_shell);
        W_dp = W_shell;
    }
    else{
        viewer.data().set_mesh(U,F);
        W_dp = W;
    }
    //morph->UpdateFrame(V);
    viewer.data().set_edges(C_widen,BE,blue);
    break;
  case '2':
    viewer.data().clear();
    viewer.data().set_mesh(V_target,F_target);
    viewer.data().set_edges(C_target,BE,red);
    break;
  case '3':
      toggle = !toggle;
    viewer.data().clear();
    if(toggle)  viewer.data().set_mesh(V_shell,F_shell);
    else viewer.data().set_mesh(U,F);
    //morph->UpdateFrame(V);
    viewer.data().set_edges(C_target,BE,red);
    break;
//  case '=':
//  {
//      cout<<"target obj? "<<flush;
//      string targetJ;
//      cin>>targetJ;
//      ReadTargetOBJ(targetJ, C_target);
//      Morph();
//      //morph->UpdateFrame(U);
//      MatrixXd V_t; MatrixXi F_t;
//      igl::readPLY(targetJ.substr(0,targetJ.size()-3)+"ply", V_t,F_t);
//      V_target = V_t; F_target = F_t;
//      viewer.data().clear();
//      viewer.data().set_mesh(V_target,F_target);
//      viewer.data().set_edges(C_target,BE,red);
//  }
    break;
//  case '2':
//    morph->GetMorphedVertices(V1);
//    viewer.data().clear();
//    viewer.data().set_mesh(V1,F1);
//    viewer.data().set_edges(C_target,BE,red);
//    break;
  case ']':
      degree_arm += 5;
      LimbWiden_arm(5,U,W);
      LimbWiden_arm(5,V_shell,W_shell);
      LimbWiden_arm(5,C_widen);
      cout<<"arm degree: "<<degree_arm<<endl;
     // morph->UpdateFrame(U);
      if(toggle)viewer.data().set_mesh(V_shell,F_shell);
      else viewer.data().set_mesh(U,F);
      viewer.data().set_edges(C_widen,BE,blue);
      break;
  case '[':
      degree_arm -= 5;
      LimbWiden_arm(-5,U,W);
      LimbWiden_arm(-5,V_shell,W_shell);
      LimbWiden_arm(-5,C_widen);
      cout<<"arm degree: "<<degree_arm<<endl;
     // morph->UpdateFrame(U);
      if(toggle)viewer.data().set_mesh(V_shell,F_shell);
      else viewer.data().set_mesh(U,F);
      viewer.data().set_edges(C_widen,BE,blue);
      break;
  case '\'':
      degree_leg += 5;
      LimbWiden_leg(5,U,W);
      LimbWiden_leg(5,V_shell,W_shell);
      LimbWiden_leg(5,C_widen);
      cout<<"leg degree: "<<degree_leg<<endl;
     // morph->UpdateFrame(U);
      if(toggle)viewer.data().set_mesh(V_shell,F_shell);
      else viewer.data().set_mesh(U,F);
      viewer.data().set_edges(C_widen,BE,blue);
      break;
  case ';':
      degree_leg -= 5;
      LimbWiden_leg(-5,U,W);
      LimbWiden_leg(-5,V_shell,W_shell);
      LimbWiden_leg(-5,C_widen);
      cout<<"leg degree: "<<degree_leg<<endl;
     // morph->UpdateFrame(U);
      if(toggle)viewer.data().set_mesh(V_shell,F_shell);
      else viewer.data().set_mesh(U,F);
      viewer.data().set_edges(C_widen,BE,blue);
      break;
  case ',':
      selected++;
      selected = std::min(std::max(selected,0),(int)W_dp.cols()-1);
      viewer.data().set_data(W_dp.col(selected));
      break;
  case '.':
      selected--;
      selected = std::min(std::max(selected,0),(int)W_dp.cols()-1);
      viewer.data().set_data(W_dp.col(selected));
      break;
  case '-':
      toggle = !toggle;
      if(toggle) W_dp = W_s;
      else W_dp = W;
      selected = std::min(std::max(selected,0),(int)W_dp.cols()-1);
      viewer.data().set_data(W_dp.col(selected));
      break;
  case '9':
     {string fileN;
      cout<<"ply file name? "<<flush;
      cin>>fileN;
      igl::writePLY(fileN,V_shell,F_shell);
//      PrintPly(fileN + ".ply");
//      morph->PrintPLY(fileN + "_f.ply");
      }
      break;
  case '/':
     {string fileN;
      cout<<"obj file name? "<<flush;
      cin>>fileN;
      PrintOBJ(fileN + ".obj", "mrcp.mtl");
      }
      break;
//  case ';':
//     {string fileN;
//      cout<<"mesh/tgf file name? "<<flush;
//      cin>>fileN;
//      igl::writeMESH(fileN+".mesh", U, T, F);
//      igl::writeTGF(fileN+".tgf",C_widen,BE);
//      }
      break;
  }
  return true;
}

void Morph(){
    //morphing - use different weights for scaling and rotation
    pair<vector<Vec3>,vector<Transfo>> scaleRot = CalculateScaleRot(C_widen, C_target, BE);
    for(int i=0;i<U.rows();i++){
        Vec3 scale_v = scaleRot.first[0]*W_s(i, 0);
        for(int j=1;j<C.rows();j++)
            scale_v = scale_v + scaleRot.first[j]*W_s(i, j);

        Transfo rot_tf = scaleRot.second[0]*W(i, 0);
        for(int j=1;j<BE.rows();j++)
            rot_tf = rot_tf + scaleRot.second[j]*W(i, j);

        Point3 p1 = rot_tf*(Point3(U(i,0),U(i,1),U(i,2))+scale_v);
        U(i, 0) = p1.x; U(i, 1) = p1.y; U(i, 2) = p1.z;
    }
    for(int i=0;i<V_shell.rows();i++){
        Vec3 scale_v = scaleRot.first[0]*W_s_shell(i, 0);
        for(int j=1;j<C.rows();j++)
            scale_v = scale_v + scaleRot.first[j]*W_s_shell(i, j);

        Transfo rot_tf = scaleRot.second[0]*W_shell(i, 0);
        for(int j=1;j<BE.rows();j++)
            rot_tf = rot_tf + scaleRot.second[j]*W_shell(i, j);

        Point3 p1 = rot_tf*(Point3(V_shell(i,0),V_shell(i,1),V_shell(i,2))+scale_v);
        V_shell(i, 0) = p1.x; V_shell(i, 1) = p1.y; V_shell(i, 2) = p1.z;
    }
}

void Morph_scale(){
    pair<vector<Vec3>,vector<Transfo>> scaleRot = CalculateScaleRot(C, C_target, BE);
    for(int i=0;i<U.rows();i++){
        Vec3 scale_v = scaleRot.first[0]*W_s(i, 0);
        for(int j=1;j<C.rows();j++)
            scale_v = scale_v + scaleRot.first[j]*W_s(i, j);
        Point3 p1 = Point3(V(i,0),V(i,1),V(i,2))+scale_v;
        U(i, 0) = p1.x; U(i, 1) = p1.y; U(i, 2) = p1.z;
    }
}

void PrintPly(string fileN){
    ofstream ofs(fileN);
    ofs<<"ply"<<endl;
    ofs<<"format ascii 1.0"<<endl;
    ofs<<"comment Exported in Morpher"<<endl;
    ofs<<"element vertex "<<m2p.size()<<endl;
    ofs<<"property float x"<<endl;
    ofs<<"property float y"<<endl;
    ofs<<"property float z"<<endl;
    ofs<<"element face "<<F.rows()<<endl;
    ofs<<"property list uchar int vertex_index"<<endl;
    ofs<<"end_header"<<endl;
    for(auto v:m2p)
        ofs<<U(v.first,0)<<" "<<U(v.first,1)<<" "<<U(v.first,2)<<endl;
    for(int i=0;i<F.rows();i++)
        ofs<<"3 "<<m2p[F(i,0)]<<" "<<m2p[F(i,1)]<<" "<<m2p[F(i,2)]<<endl;
    ofs.close();
}

void PrintOBJ(string fileN, string mtlN){
     ofstream ofs(fileN);
     ofs<<"mtllib "+mtlN<<endl<<endl;

     for(auto iter:o2m){
         int i=iter.second;
         ofs<<"v "<<U(i,0)<<" "<<U(i,1)<<" "<<U(i,2)<<endl;
     }
     ofs<<endl<<objStr;
     ofs.close();
}
void LimbWiden_arm(double degree, MatrixXd &_U, MatrixXd _W){
    double rad = degree * PI / 180.;
    vector<Transfo> tfVec; vector<Vec3> cenVec;
    for(int i=0;i<C.rows();i++)  cenVec.push_back(Vec3(C_widen(i,0),C_widen(i,1),C_widen(i,2)));
    for(int i=0;i<BE.rows();i++) tfVec.push_back(Transfo::identity());
    //C_widen = C;
    for(int i=0;i<BE.rows();i++){
        if(i==4)       tfVec[i] = Transfo::rotate(cenVec[BE(i,0)],Vec3(0,-1.,0),rad);
        else if(i==12) tfVec[i] = Transfo::rotate(cenVec[BE(i,0)],Vec3(0,1.,0),rad);
        else if(P(i)>0)  tfVec[i] = tfVec[P(i)];
    }

    for(int i=0;i<_U.rows();i++){
        Transfo tf = tfVec[0]*_W(i, 0);
        for(int j=1;j<BE.rows();j++)
            tf = tf + tfVec[j]*_W(i, j);
        Point3 p1 = tf*Point3(_U(i,0),_U(i,1),_U(i,2));
        _U(i, 0) = p1.x; _U(i, 1) = p1.y; _U(i, 2) = p1.z;
    }
}
void LimbWiden_leg(double degree, MatrixXd &_U, MatrixXd _W){
    double rad = degree * PI / 180.;
    vector<Transfo> tfVec; vector<Vec3> cenVec;
    for(int i=0;i<C.rows();i++)  cenVec.push_back(Vec3(C_widen(i,0),C_widen(i,1),C_widen(i,2)));
    for(int i=0;i<BE.rows();i++) tfVec.push_back(Transfo::identity());
    //C_widen = C;
    for(int i=0;i<BE.rows();i++){
        if(i==16) tfVec[i] = Transfo::rotate(cenVec[BE(i,0)],Vec3(0,-1.,0),rad*0.5);
        else if(i==20) tfVec[i] = Transfo::rotate(cenVec[BE(i,0)],Vec3(0,1.,0),rad*0.5);
        else if(P(i)>0)  tfVec[i] = tfVec[P(i)];
    }

    for(int i=0;i<_U.rows();i++){
        Transfo tf = tfVec[0]*_W(i, 0);
        for(int j=1;j<BE.rows();j++)
            tf = tf + tfVec[j]*_W(i, j);
        Point3 p1 = tf*Point3(_U(i,0),_U(i,1),_U(i,2));
        _U(i, 0) = p1.x; _U(i, 1) = p1.y; _U(i, 2) = p1.z;
    }
}
void LimbWiden_arm(double degree, MatrixXd &_C){
    double rad = degree * PI / 180.;
    vector<Transfo> tfVec; vector<Vec3> cenVec;
    for(int i=0;i<_C.rows();i++)  cenVec.push_back(Vec3(_C(i,0),_C(i,1),_C(i,2)));
    for(int i=0;i<BE.rows();i++) tfVec.push_back(Transfo::identity());
    //C_widen = C;
    for(int i=0;i<BE.rows();i++){
        if(i==4)       tfVec[i] = Transfo::rotate(cenVec[BE(i,0)],Vec3(0,-1.,0),rad);
        else if(i==12) tfVec[i] = Transfo::rotate(cenVec[BE(i,0)],Vec3(0,1.,0),rad);
        else if(P(i)>0)  tfVec[i] = tfVec[P(i)];
        Vec3 cen = tfVec[i]*Point3(_C(BE(i,1),0),_C(BE(i,1),1),_C(BE(i,1),2));
        _C(BE(i,1),0) = cen.x; _C(BE(i,1),1) = cen.y; _C(BE(i,1),2) = cen.z;
    }
}
void LimbWiden_leg(double degree, MatrixXd &_C){
    double rad = degree * PI / 180.;
    vector<Transfo> tfVec; vector<Vec3> cenVec;
    for(int i=0;i<_C.rows();i++)  cenVec.push_back(Vec3(_C(i,0),_C(i,1),_C(i,2)));
    for(int i=0;i<BE.rows();i++) tfVec.push_back(Transfo::identity());
    //C_widen = C;
    for(int i=0;i<BE.rows();i++){
        if(i==16) tfVec[i] = Transfo::rotate(cenVec[BE(i,0)],Vec3(0,-1.,0),rad*0.5);
        else if(i==20) tfVec[i] = Transfo::rotate(cenVec[BE(i,0)],Vec3(0,1.,0),rad*0.5);
        else if(P(i)>0)  tfVec[i] = tfVec[P(i)];
        Vec3 cen = tfVec[i]*Point3(_C(BE(i,1),0),_C(BE(i,1),1),_C(BE(i,1),2));
        _C(BE(i,1),0) = cen.x; _C(BE(i,1),1) = cen.y; _C(BE(i,1),2) = cen.z;
    }
}

int BBW_option(int argc, char *argv[]);
int BARY_option(int argc, char *argv[]);
int MORPH_option(int argc, char *argv[]);
int main(int argc, char *argv[])
{
    //arguments
    if(string(argv[1])=="-bbw")
        return BBW_option(argc, argv);
    else if(string(argv[1])=="-bary")
        return BARY_option(argc, argv);
    else
        return MORPH_option(argc, argv);
//    //read model files
// //   igl::readMESH("IR2_f.1.2.mesh",V,T,F);
//    igl::readMESH("widen15.mesh",V,T,F);
//    U=V;
// //   igl::readTGF("MRCP_AM.tgf",C,BE);
//    igl::readTGF("widen15.tgf",C,BE);
//    igl::directed_edge_parents(BE,P);


//    //GenerateMorphingWeights("morphing//  case '1':
    //    Morph();
    //    morph->UpdateFrame(U);
    //    viewer.data().clear();
    //    viewer.data().set_mesh(U,F);
    //    viewer.data().set_edges(C_target,BE,red);
    //    break;.weight", "Default.2.2.mesh", C, V, T, BE, F);
//    //CalculateScalingWeights("Scaling.weight", C, V, T);
//    //extract polygon
//    vector<int> extV;
//    for(int i=0;i<F.rows();i++){
//        extV.push_back(F(i,0));
//        extV.push_back(F(i,1));
//        extV.push_back(F(i,2));
//    }
//    sort(extV.begin(), extV.end());
//    extV.erase(unique(extV.begin(), extV.end()), extV.end());
//     for(int i=0;i<extV.size();i++) m2p[extV[i]] = i;

//    //argument
//    string targetJ = string(argv[1]);
//    ReadTargetOBJ(targetJ, C_target);

//    //Read weights (posture, scaling)
//    cout<<"Read posture weight..."<<flu//  case '1':
    //    Morph();
    //    morph->UpdateFrame(U);
    //    viewer.data().clear();
    //    viewer.data().set_mesh(U,F);
    //    viewer.data().set_edges(C_target,BE,red);
    //    break;sh;
//    ifstream ifs("IR2_f.1.m.weight");
//    double w;
//    W.resize(V.rows(), BE.rows());
//    for(int i=0;i<V.rows();i++){
//        for(int j=0;j<BE.rows();j++){
//            ifs>>w;
//            W(i,j) = w;
//        }
//    }ifs.close(); cout<<"done"<<endl;

//    ifstream ifs_s("Scaling.weight");
//    cout<<"Read scaling weight..."<<flush;
//    W_s.resize(V.rows(), C.rows());
//    for(int i=0;i<V.rows();i++){
//        for(int j=0;j<C.rows();j++){
//            ifs_s>>w;
//            W_s(i,j) = w;
//        }
//    }ifs_s.close(); cout<<"done"<<endl;

//    //morphing test
//    morph = new morpher(V,"downsampled3.ply");
//    morph->CalculateMorphWeight(W);
//    Morph2();
//    morph->SetTGF(C_target,BE);
//    morph->UpdateFrame(U);
//    string targetS = string(argv[2]);
//    morph->StartMorph(targetS);

//    igl::readPLY("downsampled3.ply",V1, F1);

//    // Plot the mesh with pseudocolors
//    igl::opengl::glfw::Viewer viewer;
//    viewer.data().set_mesh(V1, F1);
//    viewer.data().set_edges(C,BE,sea_green);
//    viewer.data().show_lines = false;
//    viewer.data().show_overlay_depth = false;
//    viewer.data().line_width = 1;
//    viewer.callback_key_down = &key_down;
//    viewer.launch();

  return EXIT_SUCCESS;
}
int BBW_option(int argc, char *argv[]){
    if(argc!=5 && argc!=6) return EXIT_FAILURE;
    string filePLY(argv[3]);
    string fileMESH(argv[2]);
    MatrixXd V1, V3, C, W_b, W_j, bc,V_w;
    MatrixXi F1, F2, T1, BE,F_w,T_w;
    VectorXi b;

    readMESH(fileMESH,V_w,T_w,F_w);
    igl::readPLY(filePLY,V1,F1);
    igl::readTGF(string(argv[4]),C,BE);

    double interval(1.);
    if(argc==6) interval = atof(argv[5]);
    MatrixXd boneP = GenerateBonePoints(C,BE,1);
    MatrixXd V2(V1.rows()+boneP.rows(),3);
    V2<<V1, boneP;
    string prefix = fileMESH.substr(0, fileMESH.size()-5);
    igl::writePLY(prefix+"_bbwF.ply",V2,F1, 0);
    if(system(("tetgen -pYqgNEF "+prefix+"_bbwF.ply").c_str())) return EXIT_FAILURE;
    if(system(("rm "+prefix+"_bbwF.1.smesh").c_str())) return EXIT_FAILURE;
    if(system(("mv "+prefix+"_bbwF.1.mesh "+prefix+"_bbwF.mesh").c_str())) return EXIT_FAILURE;

    cout<<"read "+prefix+"_bbwF.mesh..."<<flush;
    readMESH(prefix+"_bbwF.mesh",V3,T1,F2);
    cout<<"done ("<<V3.rows()<<"v, "<<T1.rows()<<"e)"<<endl;

    cout<<"Generate bary. coord. for whole phantom in frame..."<<flush;
    map<int, MatrixXd> baryCoord = GenerateBarycentricCoord(V3,T1,V_w);
    assert(baryCoord.size()==V_w.rows());
    SparseMatrix<double> bary = GenerateBarySparse(baryCoord,V3.rows(),T1);
    cout<<"done (printed "+prefix+".fbary)"<<endl<<"Printing "+prefix+".fbary..."<<flush;
    PrintBaryCoords(prefix+".fbary",baryCoord); cout<<"done"<<endl;

    igl::boundary_conditions(V3,T1,C,VectorXi(),BE,MatrixXi(),b,bc);
    igl::BBWData bbw_data;
    bbw_data.active_set_params.max_iter = 15;
    bbw_data.verbosity = 2;
    if(!igl::bbw(V3,T1,b,bc,bbw_data,W_b))  return EXIT_FAILURE;
    CleanWeights(W_b);
    igl::writeDMAT(prefix+"_fb.dmat",W_b,0);
    cout<<"Printed "+prefix+"_fb.dmat"<<endl;

    CalculateScalingWeights(C,V3,T1,W_j);
    CleanWeights(W_j);
    igl::writeDMAT(prefix+"_fj.dmat",W_j,0);
    cout<<"Printed "+prefix+"_fj.dmat"<<endl;

    MatrixXd W_w_j = bary*W_j;
    MatrixXd W_w_b = bary*W_b;

    igl::writeTGF(prefix+".tgf",C,BE);

    // Plot the mesh with pseudocolors
    cout<<endl;
    cout<<"<INFO> Check weights given to whole phantom by toggling with , and . "<<endl;
    cout<<"<INFO> Toggle between joint/bone weights by with - "<<endl<<endl;
    W = W_w_b; W_s = W_w_j;
    V=V_w; F=F_w; W_dp=W;
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_edges(C,BE,sea_green);
    viewer.data().show_lines = false;
    viewer.data().show_overlay_depth = false;
    viewer.data().line_width = 1;
    viewer.callback_key_down = &key_down;
    viewer.launch();

    return EXIT_SUCCESS;
}
int BARY_option(int argc, char *argv[]){
    string prefix(argv[2]);
    string shellOBJ(argv[3]);

    cout<<"<Tetrahedralize shell PLY>"<<endl;
    cout<<"read "+shellOBJ+"..."<<flush;
    MatrixXd V_s; MatrixXi F_s;
    MatrixXd V_fix = ReadShellOBJ(shellOBJ,V_s,F_s);
    cout<<"done ("<<V_s.rows()<<" v/"<<V_fix.rows()<<" fix)"<<endl;
    igl::readTGF(prefix+".tgf",C,BE);
    MatrixXd boneP = GenerateBonePoints(C,BE,1);
    cout<<"added "<<boneP.rows()<<" bone points"<<endl;
    MatrixXd V2(V_s.rows()+boneP.rows(),3);
    V2<<V_s, boneP;
    igl::writePLY(prefix+"_s.ply",V2,F_s, 0);
    if(system(("tetgen -pYqgNEF "+prefix+"_s.ply").c_str())) return EXIT_FAILURE;
    if(system(("rm "+prefix+"_s.1.smesh").c_str())) return EXIT_FAILURE;
    if(system(("mv "+prefix+"_s.1.mesh "+prefix+"_s.mesh").c_str())) return EXIT_FAILURE;

    string shellMESH = prefix+"_s.mesh";
    MatrixXi T_s;
    cout<<"read "+shellMESH+"..."<<flush;
    readMESH(shellMESH,V_shell,T_s,F_shell);
    cout<<"done ("<<V_shell.rows()<<"v, "<<T_s.rows()<<"t, "<<flush;
    map<int,int> f2w = CompareVertices(V_fix,V_shell);
    cout<<f2w.size()<<"f)"<<endl;
    ofstream ofs(prefix+"_s.fix"); ofs<<f2w.size()<<endl;
    for(auto iter:f2w) ofs<<iter.second<<endl; ofs.close();
    cout<<"printed "+prefix+"_s.fix";

    cout<<"read "+prefix+".mesh..."<<flush;
    readMESH(prefix+".mesh",V,T,F);
    cout<<"done ("<<V.rows()<<"v, "<<T.rows()<<"t)"<<endl;
    cout<<"read weight files "+prefix+"_fb/fj.dmat..."<<flush;
    MatrixXd Wf_b,Wf_j;
    igl::readDMAT(prefix+"_fb.dmat",Wf_b);
    igl::readDMAT(prefix+"_fj.dmat",Wf_j);
    cout<<"done"<<endl;
    cout<<"read "+prefix+"_bbwF.mesh..."<<flush;
    MatrixXd V_f; MatrixXi T_f, F_f;
    readMESH(prefix+"_bbwF.mesh",V_f,T_f,F_f);
    cout<<"done ("<<V_f.rows()<<" v/"<<T_f.rows()<<" t)"<<endl;
    cout<<"read "+prefix+".fbary..."<<flush;
    map<int, MatrixXd> baryMap = ReadBaryFile(prefix+".fbary");
    assert(baryMap.size()==V.rows());
    SparseMatrix<double> bary = GenerateBarySparse(baryMap,V_f.rows(),T_f);
    cout<<"done ("<<baryMap.size()<<")"<<endl;
    W = bary*Wf_b; W_s = bary*Wf_j;
    cout<<"done"<<endl;

    cout<<"Generate barycen. coord. for downsampled shell..."<<flush;
    map<int, MatrixXd> baryShell = GenerateBarycentricCoord(V_f,T_f,V_shell);
    PrintBaryCoords(prefix+"_s.fbary",baryShell);
    SparseMatrix<double> baryCoordShell = GenerateBarySparse(baryShell,V_f.rows(),T_f);
    cout<<"calculate shell weights..."<<flush;
    W_shell = baryCoordShell*Wf_b;
    W_s_shell = baryCoordShell*Wf_j;
    cout<<"done"<<endl;

    cout<<"Generate bary. coord. for whole phantom in shell..."<<flush;
    morph = new morpher;
    morph->SetModelV(V);
    morph->SetFrame(V_shell,T_s);
    morph->CalculateBaryCoord();
    morph->PrintBaryCoord(prefix + ".bary");
    cout<<"Printed "+prefix + ".bary"<<endl;

    cout<<"Read "+prefix+".obj..."<<flush;
    MatrixXd V1;
    objStr = readOBJ(prefix+".obj",V1);
    o2m = CompareVertices(V1,V);
    cout<<"done ("<<o2m.size()<<"/"<<V1.rows()<<")"<<endl;

    U=V;
    C_widen = C;
    W_dp = W;
    igl::directed_edge_parents(BE,P);
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_edges(C,BE,sea_green);
    viewer.data().show_lines = false;
    viewer.data().show_overlay_depth = false;
    viewer.data().line_width = 1;
    viewer.callback_key_down = &key_down;
    viewer.launch();

    delete morph;
}
int MORPH_option(int argc, char *argv[]){
    string prefix(argv[1]);
    MatrixXi T_s;
    cout<<"Read preprocessed files..."<<endl;
    igl::readTGF(prefix+".tgf",C,BE); cout<<"\t"+prefix+".tgf (C"<<C.rows()<<",BE"<<BE.rows()<<")"<<endl;
    readMESH(prefix+".mesh",V,T,F);cout<<"\t"+prefix+".mesh (V"<<V.rows()<<",T"<<T.rows()<<")"<<endl;
    readMESH(prefix+"_s.mesh", V_shell, T_s,F_shell);cout<<"\t"+prefix+"_s.mesh (V"<<V_shell.rows()<<",T"<<T_s.rows()<<flush;
    ifstream ifs(prefix+"_s.fix"); int num, fixID; set<int> fix; ifs>>num;
    for(int i=0;i<num;i++) {ifs>>fixID; fix.insert(fixID);} ifs.close();cout<<",F"<<fix.size()<<")"<<endl;
    map<int, MatrixXd> bary = ReadBaryFile(prefix+".bary");cout<<"\t"+prefix+".bary ("<<bary.size()<<")"<<endl<<endl;

    cout<<"Read target files..."<<flush;
    string target(argv[2]);
     ReadTargetOBJ(target, V_target, F_target, C_target);
    cout<<"done (V"<<V_target.rows()<<",F"<<F_target.rows()<<",C"<<C_target.rows()<<")"<<endl;
    assert(C.rows()==C_target.rows());

    cout<<"Start morphing to "+target<<endl;
    morph = new morpher;
    morph->SetModelV(V);
    cout<<"\tset "<<V.rows()<<" model vertices"<<endl;
    morph->SetFrame(V_shell,T_s);
    cout<<"\tset "<<V_shell.rows()<<"v/"<<T_s.rows()<<"t shell model vertices"<<endl;
    morph->SetShellVFix(fix);
    cout<<"\tset "<<fix.size()<<" shell vertices to fix"<<endl;
    morph->SetBaryCoord(bary);
    cout<<"\tset "<<bary.size()<<" bary. coords. for the model"<<endl;

    U=V;
    W_dp = W;
    igl::directed_edge_parents(BE,P);
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_shell, F_shell);
    viewer.data().set_edges(C,BE,sea_green);
    viewer.data().show_lines = false;
    viewer.data().show_overlay_depth = false;
    viewer.data().line_width = 1;
    viewer.callback_key_down = &key_down;
    viewer.launch();

    return EXIT_SUCCESS;
}
