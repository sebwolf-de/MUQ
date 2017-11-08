#include "MUQ/Utilities/LinearAlgebra/EigenMatrixAlgebra.h"

using namespace muq::Utilities;

EigenMatrixAlgebra::EigenMatrixAlgebra() {}

EigenMatrixAlgebra::~EigenMatrixAlgebra() {}

bool EigenMatrixAlgebra::IsEigenMatrix(std::type_info const& obj_type) {
  // is this an Eigen::Matrix type?
  return typeid(Eigen::Matrix2d)==obj_type
    || typeid(Eigen::Matrix2f)==obj_type
    || typeid(Eigen::Matrix2i)==obj_type
    || typeid(Eigen::Matrix3d)==obj_type
    || typeid(Eigen::Matrix3f)==obj_type
    || typeid(Eigen::Matrix3i)==obj_type
    || typeid(Eigen::Matrix4d)==obj_type
    || typeid(Eigen::Matrix4f)==obj_type
    || typeid(Eigen::Matrix4i)==obj_type
    || typeid(Eigen::MatrixXd)==obj_type
    || typeid(Eigen::MatrixXf)==obj_type
    || typeid(Eigen::MatrixXi)==obj_type;
}

bool EigenMatrixAlgebra::IsZero(boost::any const& obj) {
  if( typeid(Eigen::Matrix2d)==obj.type() ) { return IsZero<Eigen::Matrix2d>(obj); }
  if( typeid(Eigen::Matrix2f)==obj.type() ) { return IsZero<Eigen::Matrix2f>(obj); }
  if( typeid(Eigen::Matrix2i)==obj.type() ) { return IsZero<Eigen::Matrix2i>(obj); }

  if( typeid(Eigen::Matrix3d)==obj.type() ) { return IsZero<Eigen::Matrix3d>(obj); }
  if( typeid(Eigen::Matrix3f)==obj.type() ) { return IsZero<Eigen::Matrix3f>(obj); }
  if( typeid(Eigen::Matrix3i)==obj.type() ) { return IsZero<Eigen::Matrix3i>(obj); }

  if( typeid(Eigen::Matrix4d)==obj.type() ) { return IsZero<Eigen::Matrix4d>(obj); }
  if( typeid(Eigen::Matrix4f)==obj.type() ) { return IsZero<Eigen::Matrix4f>(obj); }
  if( typeid(Eigen::Matrix4i)==obj.type() ) { return IsZero<Eigen::Matrix4i>(obj); }

  if( typeid(Eigen::MatrixXd)==obj.type() ) { return IsZero<Eigen::MatrixXd>(obj); }
  if( typeid(Eigen::MatrixXf)==obj.type() ) { return IsZero<Eigen::MatrixXf>(obj); }
  if( typeid(Eigen::MatrixXi)==obj.type() ) { return IsZero<Eigen::MatrixXi>(obj); }

  // something went wrong
  assert(false);
  return false;
}

double EigenMatrixAlgebra::Norm(boost::any const& mat) {
  if( typeid(Eigen::Matrix2d)==mat.type() ) { return boost::any_cast<Eigen::Matrix2d const&>(mat).norm(); }
  if( typeid(Eigen::Matrix3d)==mat.type() ) { return boost::any_cast<Eigen::Matrix3d const&>(mat).norm(); }
  if( typeid(Eigen::Matrix4d)==mat.type() ) { return boost::any_cast<Eigen::Matrix4d const&>(mat).norm(); }
  if( typeid(Eigen::MatrixXd)==mat.type() ) { return boost::any_cast<Eigen::MatrixXd const&>(mat).norm(); }

  if( typeid(Eigen::Matrix2f)==mat.type() ) { return boost::any_cast<Eigen::Matrix2f const&>(mat).norm(); }
  if( typeid(Eigen::Matrix3f)==mat.type() ) { return boost::any_cast<Eigen::Matrix3f const&>(mat).norm(); }
  if( typeid(Eigen::Matrix4f)==mat.type() ) { return boost::any_cast<Eigen::Matrix4f const&>(mat).norm(); }
  if( typeid(Eigen::MatrixXf)==mat.type() ) { return boost::any_cast<Eigen::MatrixXf const&>(mat).norm(); }

  if( typeid(Eigen::Matrix2i)==mat.type() ) { return boost::any_cast<Eigen::Matrix2i const&>(mat).norm(); }
  if( typeid(Eigen::Matrix3i)==mat.type() ) { return boost::any_cast<Eigen::Matrix3i const&>(mat).norm(); }
  if( typeid(Eigen::Matrix4i)==mat.type() ) { return boost::any_cast<Eigen::Matrix4i const&>(mat).norm(); }
  if( typeid(Eigen::MatrixXi)==mat.type() ) { return boost::any_cast<Eigen::MatrixXi const&>(mat).norm(); }

  // something went wrong
  assert(false);
  return -1.0;
}

unsigned int EigenMatrixAlgebra::Size(boost::any const& mat, int const dim) {
  // fixed size vectors
  if( typeid(Eigen::Matrix2d)==mat.type() ) { return dim==-1? 4 : 2;}
  if( typeid(Eigen::Matrix3d)==mat.type() ) { return dim==-1? 9 : 3; }
  if( typeid(Eigen::Matrix4d)==mat.type() ) { return dim==-1? 16 : 4; }

  if( typeid(Eigen::Matrix2f)==mat.type() ) { return dim==-1? 4 : 2;}
  if( typeid(Eigen::Matrix3f)==mat.type() ) { return dim==-1? 9 : 3; }
  if( typeid(Eigen::Matrix4f)==mat.type() ) { return dim==-1? 16 : 4; }

  if( typeid(Eigen::Matrix2i)==mat.type() ) { return dim==-1? 4 : 2;}
  if( typeid(Eigen::Matrix3i)==mat.type() ) { return dim==-1? 9 : 3; }
  if( typeid(Eigen::Matrix4i)==mat.type() ) { return dim==-1? 16 : 4; }

  // generic Eigen::MatrixXd,f,i
  if( typeid(Eigen::MatrixXd)==mat.type() ) { return Size<Eigen::MatrixXd>(mat, dim); }
  if( typeid(Eigen::MatrixXf)==mat.type() ) { return Size<Eigen::MatrixXf>(mat, dim); }
  if( typeid(Eigen::MatrixXi)==mat.type() ) { return Size<Eigen::MatrixXi>(mat, dim); }

  // something went wrong
  assert(false);
  return 0;
}

boost::any EigenMatrixAlgebra::AccessElement(boost::any const& mat, unsigned int const i, unsigned int const j) {
  if( typeid(Eigen::Matrix2d)==mat.type() ) { return AccessElement<Eigen::Matrix2d>(mat, i, j); }
  if( typeid(Eigen::Matrix2f)==mat.type() ) { return AccessElement<Eigen::Matrix2f>(mat, i, j); }
  if( typeid(Eigen::Matrix2i)==mat.type() ) { return AccessElement<Eigen::Matrix2i>(mat, i, j); }

  if( typeid(Eigen::Matrix3d)==mat.type() ) { return AccessElement<Eigen::Matrix3d>(mat, i, j); }
  if( typeid(Eigen::Matrix3f)==mat.type() ) { return AccessElement<Eigen::Matrix3f>(mat, i, j); }
  if( typeid(Eigen::Matrix3i)==mat.type() ) { return AccessElement<Eigen::Matrix3i>(mat, i, j); }

  if( typeid(Eigen::Matrix4d)==mat.type() ) { return AccessElement<Eigen::Matrix4d>(mat, i, j); }
  if( typeid(Eigen::Matrix4f)==mat.type() ) { return AccessElement<Eigen::Matrix4f>(mat, i, j); }
  if( typeid(Eigen::Matrix4i)==mat.type() ) { return AccessElement<Eigen::Matrix4i>(mat, i, j); }

  if( typeid(Eigen::MatrixXd)==mat.type() ) { return AccessElement<Eigen::MatrixXd>(mat, i, j); }
  if( typeid(Eigen::MatrixXf)==mat.type() ) { return AccessElement<Eigen::MatrixXf>(mat, i, j); }
  if( typeid(Eigen::MatrixXi)==mat.type() ) { return AccessElement<Eigen::MatrixXi>(mat, i, j); }

  // something went wront
  assert(false);
  return boost::none;
}

boost::any EigenMatrixAlgebra::Identity(std::type_info const& type, unsigned int const rows, unsigned int const cols) {
  if( type==typeid(Eigen::Matrix2d) ) { return (Eigen::Matrix2d)Eigen::Matrix2d::Identity(); }
  if( type==typeid(Eigen::Matrix2f) ) { return (Eigen::Matrix2f)Eigen::Matrix2f::Identity(); }
  if( type==typeid(Eigen::Matrix2i) ) { return (Eigen::Matrix2i)Eigen::Matrix2i::Identity(); }
  
  if( type==typeid(Eigen::Matrix3d) ) { return (Eigen::Matrix3d)Eigen::Matrix3d::Identity(); }
  if( type==typeid(Eigen::Matrix3f) ) { return (Eigen::Matrix3f)Eigen::Matrix3f::Identity(); }
  if( type==typeid(Eigen::Matrix3i) ) { return (Eigen::Matrix3i)Eigen::Matrix3i::Identity(); }
  
  if( type==typeid(Eigen::Matrix4d) ) { return (Eigen::Matrix4d)Eigen::Matrix4d::Identity(); }
  if( type==typeid(Eigen::Matrix4f) ) { return (Eigen::Matrix4f)Eigen::Matrix4f::Identity(); }
  if( type==typeid(Eigen::Matrix4i) ) { return (Eigen::Matrix4i)Eigen::Matrix4i::Identity(); }

  if( type==typeid(Eigen::MatrixXd) ) { return (Eigen::MatrixXd)Eigen::MatrixXd::Identity(rows, cols); }
  if( type==typeid(Eigen::MatrixXf) ) { return (Eigen::MatrixXf)Eigen::MatrixXf::Identity(rows, cols); }
  if( type==typeid(Eigen::MatrixXi) ) { return (Eigen::MatrixXi)Eigen::MatrixXi::Identity(rows, cols); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenMatrixAlgebra::Add(boost::any const& in0, boost::any const& in1) {
  // 2D matrices
  if( in0.type()==typeid(Eigen::Matrix2d) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return Add<Eigen::Matrix2d, Eigen::Matrix2d>(in0, in1); }
    return Add<Eigen::Matrix2d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2f) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return Add<Eigen::Matrix2f, Eigen::Matrix2f>(in0, in1); }
    return Add<Eigen::Matrix2f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2i) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return Add<Eigen::Matrix2i, Eigen::Matrix2i>(in0, in1); }
    return Add<Eigen::Matrix2i, Eigen::MatrixXi>(in0, in1); 
  }

  // 3D matrices
  if( in0.type()==typeid(Eigen::Matrix3d) ) {
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return Add<Eigen::Matrix3d, Eigen::Matrix3d>(in0, in1); }
    return Add<Eigen::Matrix3d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3f) ) {
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return Add<Eigen::Matrix3f, Eigen::Matrix3f>(in0, in1); }
    return Add<Eigen::Matrix3f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3i) ) {
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return Add<Eigen::Matrix3i, Eigen::Matrix3i>(in0, in1); }
    return Add<Eigen::Matrix3i, Eigen::MatrixXi>(in0, in1); 
  }

  // 4D matrices
  if( in0.type()==typeid(Eigen::Matrix4d) ) {
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return Add<Eigen::Matrix4d, Eigen::Matrix4d>(in0, in1); }
    return Add<Eigen::Matrix4d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4f) ) {
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return Add<Eigen::Matrix4f, Eigen::Matrix4f>(in0, in1); }
    return Add<Eigen::Matrix4f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4i) ) {
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return Add<Eigen::Matrix4i, Eigen::Matrix4i>(in0, in1); }
    return Add<Eigen::Matrix4i, Eigen::MatrixXi>(in0, in1); 
  }

  // XD matrices
  if( in0.type()==typeid(Eigen::MatrixXd) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return Add<Eigen::MatrixXd, Eigen::Matrix2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return Add<Eigen::MatrixXd, Eigen::Matrix3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return Add<Eigen::MatrixXd, Eigen::Matrix4d>(in0, in1); }
    return Add<Eigen::MatrixXd, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXf) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return Add<Eigen::MatrixXf, Eigen::Matrix2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return Add<Eigen::MatrixXf, Eigen::Matrix3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return Add<Eigen::MatrixXf, Eigen::Matrix4f>(in0, in1); }
    return Add<Eigen::MatrixXf, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXi) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return Add<Eigen::MatrixXi, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return Add<Eigen::MatrixXi, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return Add<Eigen::MatrixXi, Eigen::Matrix4i>(in0, in1); }
    return Add<Eigen::MatrixXi, Eigen::MatrixXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenMatrixAlgebra::Subtract(boost::any const& in0, boost::any const& in1) {
  // 2D matrices
  if( in0.type()==typeid(Eigen::Matrix2d) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return Subtract<Eigen::Matrix2d, Eigen::Matrix2d>(in0, in1); }
    return Subtract<Eigen::Matrix2d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2f) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return Subtract<Eigen::Matrix2f, Eigen::Matrix2f>(in0, in1); }
    return Subtract<Eigen::Matrix2f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2i) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return Subtract<Eigen::Matrix2i, Eigen::Matrix2i>(in0, in1); }
    return Subtract<Eigen::Matrix2i, Eigen::MatrixXi>(in0, in1); 
  }

  // 3D matrices
  if( in0.type()==typeid(Eigen::Matrix3d) ) {
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return Subtract<Eigen::Matrix3d, Eigen::Matrix3d>(in0, in1); }
    return Subtract<Eigen::Matrix3d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3f) ) {
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return Subtract<Eigen::Matrix3f, Eigen::Matrix3f>(in0, in1); }
    return Subtract<Eigen::Matrix3f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3i) ) {
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return Subtract<Eigen::Matrix3i, Eigen::Matrix3i>(in0, in1); }
    return Subtract<Eigen::Matrix3i, Eigen::MatrixXi>(in0, in1); 
  }

  // 4D matrices
  if( in0.type()==typeid(Eigen::Matrix4d) ) {
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return Subtract<Eigen::Matrix4d, Eigen::Matrix4d>(in0, in1); }
    return Subtract<Eigen::Matrix4d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4f) ) {
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return Subtract<Eigen::Matrix4f, Eigen::Matrix4f>(in0, in1); }
    return Subtract<Eigen::Matrix4f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4i) ) {
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return Subtract<Eigen::Matrix4i, Eigen::Matrix4i>(in0, in1); }
    return Subtract<Eigen::Matrix4i, Eigen::MatrixXi>(in0, in1); 
  }

  // XD matrices
  if( in0.type()==typeid(Eigen::MatrixXd) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return Subtract<Eigen::MatrixXd, Eigen::Matrix2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return Subtract<Eigen::MatrixXd, Eigen::Matrix3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return Subtract<Eigen::MatrixXd, Eigen::Matrix4d>(in0, in1); }
    return Subtract<Eigen::MatrixXd, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXf) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return Subtract<Eigen::MatrixXf, Eigen::Matrix2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return Subtract<Eigen::MatrixXf, Eigen::Matrix3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return Subtract<Eigen::MatrixXf, Eigen::Matrix4f>(in0, in1); }
    return Subtract<Eigen::MatrixXf, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXi) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return Subtract<Eigen::MatrixXi, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return Subtract<Eigen::MatrixXi, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return Subtract<Eigen::MatrixXi, Eigen::Matrix4i>(in0, in1); }
    return Subtract<Eigen::MatrixXi, Eigen::MatrixXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenMatrixAlgebra::Multiply(boost::any const& in0, boost::any const& in1) {
  // 2D matrices
  if( in0.type()==typeid(Eigen::Matrix2d) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return Multiply<Eigen::Matrix2d, Eigen::Matrix2d>(in0, in1); }
    return Multiply<Eigen::Matrix2d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2f) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return Multiply<Eigen::Matrix2f, Eigen::Matrix2f>(in0, in1); }
    return Multiply<Eigen::Matrix2f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2i) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return Multiply<Eigen::Matrix2i, Eigen::Matrix2i>(in0, in1); }
    return Multiply<Eigen::Matrix2i, Eigen::MatrixXi>(in0, in1); 
  }

  // 3D matrices
  if( in0.type()==typeid(Eigen::Matrix3d) ) {
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return Multiply<Eigen::Matrix3d, Eigen::Matrix3d>(in0, in1); }
    return Multiply<Eigen::Matrix3d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3f) ) {
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return Multiply<Eigen::Matrix3f, Eigen::Matrix3f>(in0, in1); }
    return Multiply<Eigen::Matrix3f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3i) ) {
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return Multiply<Eigen::Matrix3i, Eigen::Matrix3i>(in0, in1); }
    return Multiply<Eigen::Matrix3i, Eigen::MatrixXi>(in0, in1); 
  }

  // 4D matrices
  if( in0.type()==typeid(Eigen::Matrix4d) ) {
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return Multiply<Eigen::Matrix4d, Eigen::Matrix4d>(in0, in1); }
    return Multiply<Eigen::Matrix4d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4f) ) {
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return Multiply<Eigen::Matrix4f, Eigen::Matrix4f>(in0, in1); }
    return Multiply<Eigen::Matrix4f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4i) ) {
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return Multiply<Eigen::Matrix4i, Eigen::Matrix4i>(in0, in1); }
    return Multiply<Eigen::Matrix4i, Eigen::MatrixXi>(in0, in1); 
  }

  // XD matrices
  if( in0.type()==typeid(Eigen::MatrixXd) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return Multiply<Eigen::MatrixXd, Eigen::Matrix2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return Multiply<Eigen::MatrixXd, Eigen::Matrix3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return Multiply<Eigen::MatrixXd, Eigen::Matrix4d>(in0, in1); }
    return Multiply<Eigen::MatrixXd, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXf) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return Multiply<Eigen::MatrixXf, Eigen::Matrix2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return Multiply<Eigen::MatrixXf, Eigen::Matrix3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return Multiply<Eigen::MatrixXf, Eigen::Matrix4f>(in0, in1); }
    return Multiply<Eigen::MatrixXf, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXi) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return Multiply<Eigen::MatrixXi, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return Multiply<Eigen::MatrixXi, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return Multiply<Eigen::MatrixXi, Eigen::Matrix4i>(in0, in1); }
    return Multiply<Eigen::MatrixXi, Eigen::MatrixXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenMatrixAlgebra::ScalarMultiply(boost::any const& in0, boost::any const& in1) {
  if( in0.type()==typeid(double) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return ScalarMultiply<double, Eigen::Matrix2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return ScalarMultiply<double, Eigen::Matrix3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return ScalarMultiply<double, Eigen::Matrix4d>(in0, in1); }

    return ScalarMultiply<double, Eigen::MatrixXd>(in0, in1); 
  }

  if( in0.type()==typeid(float) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return ScalarMultiply<float, Eigen::Matrix2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return ScalarMultiply<float, Eigen::Matrix3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return ScalarMultiply<float, Eigen::Matrix4f>(in0, in1); }

    return ScalarMultiply<float, Eigen::MatrixXf>(in0, in1); 
  }

  if( in0.type()==typeid(int) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return ScalarMultiply<int, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return ScalarMultiply<int, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return ScalarMultiply<int, Eigen::Matrix4i>(in0, in1); }

    return ScalarMultiply<int, Eigen::MatrixXi>(in0, in1); 
  }

  if( in0.type()==typeid(unsigned int) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return ScalarMultiply<unsigned int, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return ScalarMultiply<unsigned int, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return ScalarMultiply<unsigned int, Eigen::Matrix4i>(in0, in1); }

    return ScalarMultiply<unsigned int, Eigen::MatrixXi>(in0, in1); 
  }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenMatrixAlgebra::Zero(std::type_info const& type, unsigned int const rows, unsigned int const cols) {
  if( typeid(Eigen::Matrix2d)==type ) { return (Eigen::Matrix2d)Eigen::Matrix2d::Zero(); }
  if( typeid(Eigen::Matrix2f)==type ) { return (Eigen::Matrix2f)Eigen::Matrix2f::Zero(); }
  if( typeid(Eigen::Matrix2i)==type ) { return (Eigen::Matrix2i)Eigen::Matrix2i::Zero(); }
  
  if( typeid(Eigen::Matrix3d)==type ) { return (Eigen::Matrix3d)Eigen::Matrix3d::Zero(); }
  if( typeid(Eigen::Matrix3f)==type ) { return (Eigen::Matrix3f)Eigen::Matrix3f::Zero(); }
  if( typeid(Eigen::Matrix3i)==type ) { return (Eigen::Matrix3i)Eigen::Matrix3i::Zero(); }
  
  if( typeid(Eigen::Matrix4d)==type ) { return (Eigen::Matrix4d)Eigen::Matrix4d::Zero(); }
  if( typeid(Eigen::Matrix4f)==type ) { return (Eigen::Matrix4f)Eigen::Matrix4f::Zero(); }
  if( typeid(Eigen::Matrix4i)==type ) { return (Eigen::Matrix4i)Eigen::Matrix4i::Zero(); }

  if( typeid(Eigen::MatrixXd)==type ) { return (Eigen::MatrixXd)Eigen::MatrixXd::Zero(rows, cols); }
  if( typeid(Eigen::MatrixXf)==type ) { return (Eigen::MatrixXf)Eigen::MatrixXf::Zero(rows, cols); }
  if( typeid(Eigen::MatrixXi)==type ) { return (Eigen::MatrixXi)Eigen::MatrixXi::Zero(rows, cols); }

  // something went wrong
  assert(false);
  return boost::none;
}
