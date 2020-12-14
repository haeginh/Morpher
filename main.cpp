// Because of Mosek complications, we don't use static library if Mosek is used.
#ifdef LIBIGL_WITH_MOSEK
#ifdef IGL_STATIC_LIBRARY
#undef IGL_STATIC_LIBRARY
#endif
#endif

#include <igl/boundary_conditions.h>
#include <igl/colon.h>
#include <igl/column_to_quats.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/jet.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/normalize_row_sums.h>
#include <igl/readDMAT.h>
#include <igl/readMESH.h>
#include <igl/readTGF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/bbw.h>
//#include <igl/embree/bone_heat.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <iostream>

//#include "tutorial_shared_path.h"
#include "functions.h"
#include "G4Tet.hh"
#include "G4ThreeVector.hh"

typedef
  std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> >
  RotationList;

const Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);
int selected = 0;
Eigen::MatrixXd V,W,U,C,M;
Eigen::MatrixXi T,F,BE;
Eigen::VectorXi P;
RotationList pose;
double anim_t = 1.0;
double anim_t_dir = -0.03;

bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
  if(viewer.core().is_animating)
  {
    // Interpolate pose and identity
    RotationList anim_pose(pose.size());
    for(int e = 0;e<pose.size();e++)
    {
      anim_pose[e] = pose[e].slerp(anim_t,Quaterniond::Identity());
    }
    // Propagate relative rotations via FK to retrieve absolute transformations
    RotationList vQ;
    vector<Vector3d> vT;
    igl::forward_kinematics(C,BE,P,anim_pose,vQ,vT);
    const int dim = C.cols();
    MatrixXd T(BE.rows()*(dim+1),dim);
    for(int e = 0;e<BE.rows();e++)
    {
      Affine3d a = Affine3d::Identity();
      a.translate(vT[e]);
      a.rotate(vQ[e]);
      T.block(e*(dim+1),0,dim+1,dim) =
        a.matrix().transpose().block(0,0,dim+1,dim);
    }
    // Compute deformation via LBS as matrix multiplication
    U = M*T;

    // Also deform skeleton edges
    MatrixXd CT;
    MatrixXi BET;
    igl::deform_skeleton(C,BE,T,CT,BET);

    viewer.data().set_vertices(U);
    viewer.data().set_edges(CT,BET,sea_green);
    viewer.data().compute_normals();
    anim_t += anim_t_dir;
    anim_t_dir *= (anim_t>=1.0 || anim_t<=0.0?-1.0:1.0);
  }
  return false;
}

#include "dual_quat_cu.hpp"
using namespace Tbx;
void dual_quat_deformer();
void dual_quat_deformer(vector<Dual_quat_cu>);
void LBS_deformer();
void PrintPLY();
void PrintPLY(int);
void ReadOBJFiles();
void PrintOBJ(string);
void PrintKINECTinput(string modelN);
vector<vector<Dual_quat_cu>> ReadBVH(string fileN);
vector<vector<Dual_quat_cu>> dual_quat_vec;
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case ' ':
      viewer.core().is_animating = !viewer.core().is_animating;
      break;
    case '.':
      selected++;
      selected = std::min(std::max(selected,0),(int)W.cols()-1);
      viewer.data().set_data(W.col(selected));
      break;
    case ',':
      selected--;
      selected = std::min(std::max(selected,0),(int)W.cols()-1);
      viewer.data().set_data(W.col(selected));
      break;
    case '/':
      dual_quat_deformer();
      viewer.data().set_mesh(U,F);
      break;
    case '5':
      LBS_deformer();
      viewer.data().set_mesh(U,F);
      break;
    case '1':
      PrintPLY();
      break;
    case '2':
      dual_quat_vec.clear();
      while(dual_quat_vec.size()==0){
          string fileN;
          cout<<"bvh file name? ";cin>>fileN;
          dual_quat_vec = ReadBVH(fileN);
      }
      break;
    case '3':
     {int frameNo;
      cout<<"frame No.? ";cin>>frameNo;
      dual_quat_deformer(dual_quat_vec[frameNo]);}
      viewer.data().set_vertices(U);
      viewer.data().compute_normals();
      break;
    case '6'://print all for bvh
      for(size_t i=0; i<dual_quat_vec.size(); i++){
          dual_quat_deformer(dual_quat_vec[i]);
          PrintPLY(i);
          viewer.data().set_vertices(U);
          viewer.data().compute_normals();
      }
      break;
    case '7'://print all obj for bvh
      ReadOBJFiles();
      for(size_t i=0; i<dual_quat_vec.size(); i++){
          dual_quat_deformer(dual_quat_vec[i]);
          PrintOBJ(to_string(i));
          viewer.data().set_vertices(U);
          viewer.data().compute_normals();
      }
      break;
    case '4':
     {string fileN;
      cout<<"Kinect input file name? ";cin>>fileN;
      PrintKINECTinput(fileN);}
  }
  return true;
}

map<int, int> m2p;
void dual_quat_deformer()
{
    vector<Dual_quat_cu> dual_quat; dual_quat.reserve(P.rows());
    for(int i=0;i<P.rows();i++)
        dual_quat[i] = Dual_quat_cu::identity();
    double x,y,z;
    int id; double angle;
    cout<<"id: "; cin>>id;
    cin.clear();cin.ignore(256, '\n');
    std::string axisStr;
    std::cout<<"axis: ";
    std::getline(std::cin, axisStr);
    stringstream ss(axisStr); ss>>x>>y>>z;
    cout<<"degree: "; cin>>angle; angle *= 3.14159/180;
    Transfo tf = Transfo::rotate(Vec3(C(BE(id,0),0),C(BE(id,0),1),C(BE(id,0),2)), Vec3(x,y,z), angle);
    dual_quat[id] = Dual_quat_cu(tf);
    for(int i=id+1;i<P.size();i++){
        if(P(i)<0) continue;
        dual_quat[i] = dual_quat[P[i]];
    }
    cout<<"Start deformation.."<<flush;
    U=V;
    for(auto p:m2p)
    {
        Dual_quat_cu dq_blend;
        bool first(true);
        Quat_cu q0;

        for(int i=0;i<W.cols();i++){
            if(first){
                dq_blend = dual_quat[i] * W(p.first,i);
                q0 = dual_quat[i].rotation();
                first = false;
                continue;
            }
            if( dual_quat[i].rotation().dot( q0 ) < 0.f )
                dq_blend = dq_blend + dual_quat[i] * (-W(p.first,i));
            else dq_blend = dq_blend + dual_quat[i] * W(p.first,i);
        }

        // Compute animated position
        Point3 vi = dq_blend.transform(Point3(V(p.first,0), V(p.first,1), V(p.first,2)));
        U(p.first,0) = vi.x;
        U(p.first,1) = vi.y;
        U(p.first,2) = vi.z;
    }
    cout<<"done"<<endl;
}

void LBS_deformer()
{
    vector<Transfo> tf_data; tf_data.reserve(P.rows());
    for(int i=0;i<P.rows();i++)
        tf_data[i] = Transfo::identity();
    double x,y,z;
    int id; double angle;
    cout<<"id: "; cin>>id;
    cin.clear();cin.ignore(256, '\n');
    std::string axisStr;
    std::cout<<"axis: ";
    std::getline(std::cin, axisStr);
    stringstream ss(axisStr); ss>>x>>y>>z;
    cout<<"degree: "; cin>>angle; angle *= 3.14159/180;
    Transfo tf = Transfo::rotate(Vec3(C(BE(id,0),0),C(BE(id,0),1),C(BE(id,0),2)), Vec3(x,y,z), angle);
    tf_data[id] = tf;
    for(int i=id+1;i<P.size();i++){
        if(P(i)<0) continue;
        tf_data[i] = tf_data[P[i]];
    }
    cout<<"Start deformation.."<<flush;
    U=V;
    for(auto p:m2p)
    {
        Transfo tf_blend;
        bool first(true);

        for(int i=0;i<W.cols();i++){
            if(first){
                tf_blend = tf_data[i] * W(p.first,i);
                first = false;
                continue;
            }
            tf_blend = tf_blend + tf_data[i] * W(p.first,i);
        }

        // Compute animated position
        Point3 vi = tf_blend*Point3(V(p.first,0), V(p.first,1), V(p.first,2));
        U(p.first,0) = vi.x;
        U(p.first,1) = vi.y;
        U(p.first,2) = vi.z;
    }
    cout<<"done"<<endl;
}

void PrintPLY(){
    string fileN;
    cout<<"*.PLY File name? "<<flush; cin>>fileN;
    ofstream ofs(fileN+".ply");
    ofs<<"ply"<<endl;
    ofs<<"format ascii 1.0"<<endl;
    ofs<<"comment Exported in BBW"<<endl;
    ofs<<"element vertex "<<m2p.size()<<endl;
    ofs<<"property float x"<<endl;
    ofs<<"property float y"<<endl;
    ofs<<"property float z"<<endl;
    ofs<<"element face "<<F.rows()<<endl;
    ofs<<"property list uchar int vertex_index"<<endl;
    ofs<<"end_header"<<endl;
    for(auto p:m2p){
        ofs<<U(p.first,0)<<" "<<U(p.first,1)<<" "<<U(p.first,2)<<endl;
    }
    for(int i=0;i<F.rows();i++)
        ofs<<"3 "<<m2p[F(i,0)]<<" "<<m2p[F(i,1)]<<" "<<m2p[F(i,2)]<<endl;
    ofs.close();
    cout<<fileN+".ply was exported"<<endl;
}

void PrintPLY(int i){
    ofstream ofs(to_string(i)+".ply");
    ofs<<"ply"<<endl;
    ofs<<"format ascii 1.0"<<endl;
    ofs<<"comment Exported in BBW"<<endl;
    ofs<<"element vertex "<<m2p.size()<<endl;
    ofs<<"property float x"<<endl;
    ofs<<"property float y"<<endl;
    ofs<<"property float z"<<endl;
    ofs<<"element face "<<F.rows()<<endl;
    ofs<<"property list uchar int vertex_index"<<endl;
    ofs<<"end_header"<<endl;
    for(auto p:m2p){
        ofs<<U(p.first,0)<<" "<<U(p.first,1)<<" "<<U(p.first,2)<<endl;
    }
    for(int i=0;i<F.rows();i++)
        ofs<<"3 "<<m2p[F(i,0)]<<" "<<m2p[F(i,1)]<<" "<<m2p[F(i,2)]<<endl;
    ofs.close();
    cout<<to_string(i)+".ply was exported"<<endl;
}

stringstream ss_obj;
vector<int> objID;

void ReadOBJFiles(){
    string dump;
    ifstream ifs("obj_face.txt");
    while(getline(ifs,dump)){
        ss_obj<<dump<<endl;
    }ifs.close();

    ifstream ifs2("obj_conv.txt");
    int id;
    while(ifs2>>id) objID.push_back(id);
    ifs2.close();
}

void PrintOBJ(string file){
    ofstream ofs(file+".obj");
    ofs<<"mtllib IR2_o2.mtl"<<endl<<endl;
    for(auto id:objID)
        ofs<<"v "<<U(id,0)<<" "<<U(id,1)<<" "<<U(id,2)<<endl;

    ofs<<endl;
    ofs<<ss_obj.str();
    ofs.close();
}

int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;
    //argument
    string modelN = string(argv[1]);
//   EraseInnerTri(modelN+".mesh");
    //GenerateTGF("KinectJoint", "mrcp.joint");

    igl::readMESH(modelN+".2.mesh",V,T,F);
    U=V;

    igl::readTGF("MRCP_AM.tgf",C,BE);
//    igl::readTGF("mrcp.tgf",C,BE);
//      igl::readTGF("skinUp.tgf",C,BE);
    // retrieve parents for forward kinematics
    igl::directed_edge_parents(BE,P);

    int opt(0);
    cout<<"1) Use BBW "<<endl;
    cout<<"2) Read m.weight file "<<endl;
    cout<<"3) Broadcast from other mesh "<<endl;
    cout<<"Choose weight calculation method ->"; cin>>opt;

    if(opt==1){
        // List of boundary indices (aka fixed value indices into VV)
        VectorXi b;
        // List of boundary conditions of each weight function
        MatrixXd bc;
        igl::boundary_conditions(V,T,C,VectorXi(),BE,MatrixXi(),b,bc);

        // compute BBW weights matrix
        igl::BBWData bbw_data;
        // only a few iterations for sake of demo
        bbw_data.active_set_params.max_iter = 8;
        bbw_data.verbosity = 2;
        if(!igl::bbw(V,T,b,bc,bbw_data,W))
        {
            return EXIT_FAILURE;
        }

        // Normalize weights to sum to one
        igl::normalize_row_sums(W,W);
        ofstream ofsWGT_m(modelN+".m.weight");
        ofsWGT_m<<W; ofsWGT_m.close();
        cout<<"Printed weight for mesh phantom"<<endl;
    }
    else if(opt==2){
        string fileN;
        ifstream ifsW;
        while(!ifsW.is_open()){
            cout<<"*.m.weight file name? "; cin>>fileN;
            ifsW.open(fileN);
        }
        double tmp;
        W.resize(V.rows(), BE.rows());
        std::vector<int> notAssigned;
        ofstream ofss("tmp.pts");
        for(int r=0;r<V.rows();r++){
            for(int c=0;c<BE.rows();c++){
                ifsW>>tmp;
                W(r,c) = tmp;
            }
            if(W.row(r).sum()<0.9){
                ofss<<V(r,0)<<" "<<V(r,1)<<" "<<V(r,2)<<endl;
                notAssigned.push_back(r);
            }
        }ifsW.close(); ofss.close();
        cout<<"imported "<<fileN<<" ("<<notAssigned.size()<<" X)"<<endl;
//        cout<<"finding neighbors for not assigned points.."<<flush;

//        std::map<int, vector<int>> neighbors;
//        for(auto v:notAssigned) neighbors[v] = {};
//        for(int i=0;i<T.rows();i++){
//            if(neighbors.find(T(i, 0))!=neighbors.end()){
//                neighbors[T(i, 0)].push_back(T(i,1));
//                neighbors[T(i, 0)].push_back(T(i,2));
//                neighbors[T(i, 0)].push_back(T(i,3));
//            }
//            else if(neighbors.find(T(i, 1))!=neighbors.end()){
//                neighbors[T(i, 1)].push_back(T(i,0));
//                neighbors[T(i, 1)].push_back(T(i,2));
//                neighbors[T(i, 1)].push_back(T(i,3));
//            }
//            else if(neighbors.find(T(i, 2))!=neighbors.end()){
//                neighbors[T(i, 2)].push_back(T(i,0));
//                neighbors[T(i, 2)].push_back(T(i,1));
//                neighbors[T(i, 2)].push_back(T(i,3));
//            }
//            else if(neighbors.find(T(i, 3))!=neighbors.end()){
//                neighbors[T(i, 3)].push_back(T(i,0));
//                neighbors[T(i, 3)].push_back(T(i,1));
//                neighbors[T(i, 3)].push_back(T(i,2));
//            }
//        }
//        cout<<"done"<<endl;
//        for(auto vec:neighbors){
//            sort(vec.second.begin(),vec.second.end());
//            vec.second.erase(unique(vec.second.begin(), vec.second.end()), vec.second.end());
//            vector<double> distInv;
//            double sum(0);
//            for(auto n:vec.second){
//                double dInv = 1./sqrt((V(vec.first,0)-V(n,0))*(V(vec.first,0)-V(n,0))+
//                                      (V(vec.first,1)-V(n,1))*(V(vec.first,1)-V(n,1))+
//                                      (V(vec.first,2)-V(n,2))*(V(vec.first,2)-V(n,2)));
//                distInv.push_back(dInv);
//                sum += dInv;
//            }
//            for(auto &d:distInv) d /= sum;
//            for(int i=0;i<W.cols();i++) W(vec.first,i) = 0;
//            for(int i=0;i<vec.second.size();i++) {
//                W.row(vec.first) = W.row(vec.first) + distInv[i]*W.row(vec.second[i]);
//                cout<<distInv[i]<<W.row(vec.second[i]).sum()<<endl;
//            }
//            cout<<W.row(vec.first).sum()<<endl;
////            cout<<"\rAssigning..."<<++count<<"/"<<notAssigned.size()<<flush;
//        }
//        ofstream ofsWGT_m(modelN+".m.weight2");
//        ofsWGT_m<<W; ofsWGT_m.close();
//        cout<<endl<<"Printed "+modelN+".m.weight2"<<endl;
    }
    else if(opt==3){
        string fileN;
        ifstream ifsW;
        while(!ifsW.is_open()){
            cout<<"model name? "; cin>>fileN;
            ifsW.open(fileN+".m.weight");
        }
        Eigen::MatrixXd V1,W1;
        Eigen::MatrixXi T1,F1;
        igl::readMESH(fileN+".mesh",V1,T1,F1);
        cout<<"imported "+fileN+".mesh"<<endl;
        double tmp;
        W1.resize(V1.rows(), BE.rows());
        for(int r=0;r<V1.rows();r++){
            for(int c=0;c<BE.rows();c++){
                ifsW>>tmp;
                W1(r,c) = tmp;
            }
        }ifsW.close();
        cout<<"imported "+fileN+".m.weight"<<endl;
        vector<G4ThreeVector> vVec;
        for(int i=0;i<V.rows();i++)
            vVec.push_back(G4ThreeVector(V(i,0),V(i,1),V(i,2)));
        int count(0);
        W.resize(V.rows(),BE.rows());
        for(int r=0;r<T1.rows();r++){
            int a=T1(r,0);
            int b=T1(r,1);
            int c=T1(r,2);
            int d=T1(r,3);
            G4ThreeVector p0(V1(a,0),V1(a,1),V1(a,2));
            G4ThreeVector p1(V1(b,0),V1(b,1),V1(b,2));
            G4ThreeVector p2(V1(c,0),V1(c,1),V1(c,2));
            G4ThreeVector p3(V1(d,0),V1(d,1),V1(d,2));
            G4ThreeVector n0 = -(p2-p1).cross(p3-p1);
            G4ThreeVector n1 = (p3-p2).cross(p0-p2);
            G4ThreeVector n2 = -(p0-p3).cross(p1-p3);
            G4ThreeVector n3 = (p1-p0).cross(p2-p0);
            double d0 = (p1-p0).dot(n0);
            double d1 = (p2-p1).dot(n1);
            double d2 = (p3-p2).dot(n2);
            double d3 = (p0-p3).dot(n3);
            G4Tet tet(".",p0,p1,p2,p3);
            int id(-1);
            for(G4ThreeVector v:vVec){
                id++;
                if(tet.Inside(v)==kOutside) continue;
                if(W.row(id).sum()>0.9) continue;
                if((v-p0).mag2()<1e-10||(v-p1).mag2()<1e-10||(v-p2).mag2()<1e-10||(v-p3).mag2()<1e-10){
                    vector<int> nVec = {a,b,c,d};
                    for(auto n:nVec){
                        if(fabs(V(id,0)-V1(n,0))>0.0001) continue;
                        if(fabs(V(id,1)-V1(n,1))>0.0001) continue;
                        if(fabs(V(id,2)-V1(n,2))>0.0001) continue;
                        W.row(id) = W1.row(n);
                        count++;
                        break;
                    }
                    continue;
                }
                double dd0 = d0/(v-p0).dot(n0);
                double dd1 = d1/(v-p1).dot(n1);
                double dd2 = d2/(v-p2).dot(n2);
                double dd3 = d3/(v-p3).dot(n3);
                double dist0 = (v-p0).mag()*dd0;
                double dist1 = (v-p1).mag()*dd1;
                double dist2 = (v-p2).mag()*dd2;
                double dist3 = (v-p3).mag()*dd3;
                double w0 = (dist0-(v-p0).mag())/dist0;
                double w1 = (dist1-(v-p1).mag())/dist1;
                double w2 = (dist2-(v-p2).mag())/dist2;
                double w3 = (dist3-(v-p3).mag())/dist3;
                count++;
                //grad[id] = {w0/wSum, w1/wSum, w2/wSum, w3/wSum};
                W.row(id)=w0*W1.row(a)+w1*W1.row(b)+w2*W1.row(c)+w3*W1.row(d);
            }
            cout<<"\rAssigned points.."+to_string(count)+"/"+to_string(vVec.size())<<flush;
        }
        cout<<endl;
        ofstream ofsWGT_m(modelN+".m.weight");
        ofsWGT_m<<W; ofsWGT_m.close();
        cout<<"Printed weight for mesh phantom"<<endl;
    }

    //print input files for DQS
    //map<int, int> m2p;
    for(int i=0;i<F.rows();i++){
        m2p[F(i,0)]=0;
        m2p[F(i,1)]=0;
        m2p[F(i,2)]=0;
    }
//    cout<<"Extract "<<m2p.size()<<"/"<<V.rows()<<endl;
//    cout<<"Print "+modelN+".obj/weight.."<<flush;
//    map<int, int> conv2Kinect;
//    conv2Kinect[0] = 0;
//    conv2Kinect[1] = 1;
//    conv2Kinect[2] = 2;
//    conv2Kinect[3] = 2;
//    conv2Kinect[4] = 4;
//    conv2Kinect[5] = 5;
//    conv2Kinect[6] = 6;
//    conv2Kinect[7] = 7;
//    conv2Kinect[8] = 2;
//    conv2Kinect[9] = 11;
//    conv2Kinect[10] = 12;
//    conv2Kinect[11] = 13;
//    conv2Kinect[12] = 14;
//    conv2Kinect[13] = 0;
//    conv2Kinect[14] = 18;
//    conv2Kinect[15] = 19;
//    conv2Kinect[16] = 20;
//    conv2Kinect[17] = 0;
//    conv2Kinect[18] = 22;
//    conv2Kinect[19] = 23;
//    conv2Kinect[20] = 24;
//    conv2Kinect[21] = 3;
//    conv2Kinect[22] = 26;

//    ofstream ofsOBJ(modelN+".obj");
//    ofstream ofsWGT(modelN+".weight");
//    double eps(1e-7);
    int count(0);
    for(auto &p:m2p){
//        ofsWGT<<count;
//        map<int, double> weights;
//        for(auto i=0;i<W.cols();i++){
//            if(W(p.first,i)>eps) weights[conv2Kinect[i]]+=W(p.first,i);
//        }
//        for(auto w:weights) ofsWGT<<" "<<w.first<<" "<<w.second;
//        ofsWGT<<endl;
        p.second = count++;
//        ofsOBJ<<"v "<<V(p.first,0)<<" "<<V(p.first,1)<<" "<<V(p.first,2)<<endl;
    }//ofsWGT.close();
//    cout<<"weight done"<<flush;
//    for(int i=0;i<F.rows();i++){
//        ofsOBJ<<"f "<<m2p[F(i,0)]<<" "<<m2p[F(i,1)]<<" "<<m2p[F(i,2)]<<endl;
//    }ofsOBJ.close();
//    cout<<"\rPrint "+modelN+".obj/weight..all done     "<<endl;


    // Plot the mesh with pseudocolors
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(U, F);
    //viewer.data().set_data(W.col(selected));
    viewer.data().set_edges(C,BE,sea_green);
    viewer.data().show_lines = false;
    viewer.data().show_overlay_depth = false;
    viewer.data().line_width = 1;
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &key_down;
    viewer.core().is_animating = false;
    viewer.core().animation_max_fps = 30.;
    cout<<
      "Press '.' to show next weight function."<<endl<<
      "Press ',' to show previous weight function."<<endl<<
      "Press [space] to toggle animation."<<endl;
    viewer.launch();

/*
  igl::readMESH("hand.mesh",V,T,F);
  U=V;
  igl::readTGF("hand.tgf",C,BE);
  // retrieve parents for forward kinematics
  igl::directed_edge_parents(BE,P);

  // Read pose as matrix of quaternions per row
  MatrixXd Q;
  igl::readDMAT("hand-pose.dmat",Q);
  igl::column_to_quats(Q,pose);
  assert(pose.size() == BE.rows());

  // List of boundary indices (aka fixed value indices into VV)
  VectorXi b;
  // List of boundary conditions of each weight function
  MatrixXd bc;
  igl::boundary_conditions(V,T,C,VectorXi(),BE,MatrixXi(),b,bc);

  // compute BBW weights matrix
  igl::BBWData bbw_data;
  // only a few iterations for sake of demo
  bbw_data.active_set_params.max_iter = 8;
  bbw_data.verbosity = 2;
  if(!igl::bbw(V,T,b,bc,bbw_data,W))
  {
    return EXIT_FAILURE;
  }

  //MatrixXd Vsurf = V.topLeftCorner(F.maxCoeff()+1,V.cols());
  //MatrixXd Wsurf;
  //if(!igl::bone_heat(Vsurf,F,C,VectorXi(),BE,MatrixXi(),Wsurf))
  //{
  //  return false;
  //}
  //W.setConstant(V.rows(),Wsurf.cols(),1);
  //W.topLeftCorner(Wsurf.rows(),Wsurf.cols()) = Wsurf = Wsurf = Wsurf = Wsurf;

  // Normalize weights to sum to one
  igl::normalize_row_sums(W,W);
  // precompute linear blend skinning matrix
  igl::lbs_matrix(V,W,M);
*/
  return EXIT_SUCCESS;
}

void dual_quat_deformer(vector<Dual_quat_cu> dual_quat)
{
    cout<<"Start deformation.."<<flush;
    U=V;
 //   for(auto p:m2p)
    for(auto id:objID)
    {
//        int id = p.first;
        Dual_quat_cu dq_blend;
        bool first(true);
        Quat_cu q0;

  //      cout<<p.first<<endl;
        for(int i=0;i<W.cols();i++){
  //          cout<<W(p.first,i)<<" "<<dual_quat[i].get_non_dual_part()<<" "<<dual_quat[i].get_dual_part()<<" "<<endl;
            if(first){
                dq_blend = dual_quat[i] * W(id,i);
                q0 = dual_quat[i].rotation();
                first = false;
                continue;
            }
            if( dual_quat[i].rotation().dot( q0 ) < 0.f )
                dq_blend = dq_blend + dual_quat[i] * (-W(id,i));
            else dq_blend = dq_blend + dual_quat[i] * W(id,i);
        }

        // Compute animated position
        Point3 vi = dq_blend.transform(Point3(V(id,0), V(id,1), V(id,2)));
        U(id,0) = vi.x;
        U(id,1) = vi.y;
        U(id,2) = vi.z;
    }
    cout<<"done"<<endl;
}

vector<vector<Dual_quat_cu>> ReadBVH(string fileN){
    vector<vector<Transfo>> data;
    vector<vector<Dual_quat_cu>> conv_data;

    ifstream ifs(fileN);
    if(!ifs.is_open()) return conv_data;

    string dump;
    while(getline(ifs,dump)){
        stringstream ss(dump);
        ss>>dump;
        if(dump=="MOTION") break;
    }
    int frameNo; double frameT;
    ifs>>dump>>frameNo;
    ifs>>dump>>dump>>frameT;
    cout<<"frame #: "<<frameNo<<" / frame time: "<<frameT<<endl;
    data.resize(frameNo);
    double x,y,z;

    for(int i=0;i<frameNo;i++){
        data[i].resize(23);
        vector<Transfo> bvhData;
        ifs>>x>>y>>z;
        Transfo trans = Transfo::translate(x,y,z);

        //collect a row
        for(int j=0;j<19;j++){
            ifs>>y>>x>>z;
            if(j==5) bvhData.push_back(transfo_from_eulerYXZ(y,z,x));
            else if(j==11) bvhData.push_back(transfo_from_eulerYXZ(y,-z,-x));
 //           else if(j==6) bvhData.push_back(transfo_from_eulerYXZ(y,z,x));
 //           else if(j==12) bvhData.push_back(transfo_from_eulerYXZ(y,z,x));
            else bvhData.push_back(transfo_from_eulerYXZ(y,x,z));

        }

//        bvhData[1] = Transfo::identity();
//        bvhData[2] = Transfo::identity();
        //bvhData[3] = Transfo::identity();
        //bvhData[9] = Transfo::identity();

        Transfo tf4 = Transfo::rotate(Vec3(-0.2,-1,0),85*deg);
        Transfo tf10 = Transfo::rotate(Vec3(-0.2,1,0),85*deg);
        Transfo tf6 = Transfo::rotate(Vec3(0,1,0),80*deg);
        Transfo tf12 = Transfo::rotate(Vec3(0,1,0),-80*deg);
        bvhData[4]  = bvhData[4]*tf4;
        bvhData[6]  = bvhData[6]*tf6;
        bvhData[10]  = bvhData[10]*tf10;
        bvhData[12]  = bvhData[12]*tf12;

        Transfo tf0 = trans*bvhData[0];
//        Transfo tf0 = bvhData[0];
        for(int j=0;j<BE.rows();j++){
            int p = BE(j,0);
            if(P(j)<0) {
                data[i][j] = Transfo::identity();
                continue;
            }
            Transfo tf = Transfo::translate(Vec3(C(p,0),C(p,1),C(p,2)))*bvhData[p]*Transfo::translate(-Vec3(C(p,0),C(p,1),C(p,2)));
            data[i][j]=data[i][P(j)]*tf;
//            if(p==4) data[i][j] = data[i][j];
        }
        for(Transfo &tf:data[i]){
            tf = tf0*tf;
        }
    }

    for(auto d:data){
        vector<Dual_quat_cu> du_quat;
        for(auto tf:d){
            du_quat.push_back(Dual_quat_cu(tf));
        }
        conv_data.push_back(du_quat);
    }
    return conv_data;
}

void PrintKINECTinput(string fileN){
    cout<<"Print "+fileN+".ply file.."<<flush;
    ofstream ofs(fileN+".ply");
    ofs<<"ply"<<endl;
    ofs<<"format ascii 1.0"<<endl;
    ofs<<"comment Exported in BBW"<<endl;
    ofs<<"element vertex "<<m2p.size()<<endl;
    ofs<<"property float x"<<endl;
    ofs<<"property float y"<<endl;
    ofs<<"property float z"<<endl;
    ofs<<"element face "<<F.rows()<<endl;
    ofs<<"property list uchar int vertex_index"<<endl;
    ofs<<"end_header"<<endl;
    for(auto v:m2p)
        ofs<<V(v.first,0)<<" "<<V(v.first,1)<<" "<<V(v.first,2)<<endl;
    for(int i=0;i<F.rows();i++)
        ofs<<"3 "<<m2p[F(i,0)]<<" "<<m2p[F(i,1)]<<" "<<m2p[F(i,2)]<<endl;
    ofs.close();
    cout<<"done"<<endl;
    cout<<"Read KinectJoint file.."<<flush;
    ifstream ifs("KinectJoint");
    string dump; int d,p;
    vector<int> bone2kinect;
    while(ifs>>d>>dump>>p){
        if(p<0) continue;
        bone2kinect.push_back(p);
    }ifs.close();
    cout<<"done"<<endl;
    cout<<"Print "+fileN+".weight file.."<<flush;
    ofstream ofsW(fileN+".weight");
    int id(0);
    bone2kinect[21]=26;
    for(auto v:m2p){
        map<int, double> w;
        double epsilon(1e-7),sum(0);
        for(int i=0;i<BE.rows();i++){
            double wVal = W(v.first,i);
            if(wVal<epsilon) continue;
            w[bone2kinect[i]] = w[bone2kinect[i]]+wVal;
            sum+=wVal;
        }
        ofsW<<id++;
        for(auto ww:w)
           ofsW<<" "<<ww.first<<" "<<ww.second/sum;
        ofsW<<endl;
    }ofs.close();
    cout<<"done"<<endl;
}
