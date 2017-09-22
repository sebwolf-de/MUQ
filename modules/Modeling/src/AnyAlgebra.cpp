#include "MUQ/Modeling/AnyAlgebra.h"

using namespace muq::Modeling;

AnyAlgebra::AnyAlgebra() {}

unsigned int AnyAlgebra::VectorDimensionBase(boost::any const& vec) const {
  if( doubleType.compare(vec.type().name())==0 ) { // the type is a double
    return 1;
  }

  if( eigenVec2Type.compare(vec.type().name())==0 ) { // the type is an Eigen::Vector2d
    // return the size
    return 2;
  }

  if( eigenVec3Type.compare(vec.type().name())==0 ) { // the type is an Eigen::Vector3d
    // return the size
    return 3;
  }

  if( eigenVec4Type.compare(vec.type().name())==0 ) { // the type is an Eigen::Vector4d
    // return the size
    return 4;
  }

  if( eigenVecType.compare(vec.type().name())==0 ) { // the type is an Eigen::VectorXd
    // get a constant reference to the Eigen::VectorXd
    const Eigen::VectorXd& vecref = boost::any_cast<const Eigen::VectorXd&>(vec);
    
    // return the size
    return vecref.size();
  }

#if MUQ_HAS_SUNDIALS==1
  if( N_VectorType.compare(vec.type().name())==0 ) { // the type is a N_Vector
    // get a constant reference to the N_Vector
    const N_Vector& vecref = boost::any_cast<const N_Vector&>(vec);
    
    // return the size
    return NV_LENGTH_S(vecref);
  }
#endif

  return VectorDimension(vec);
}

boost::any AnyAlgebra::ZeroVectorBase(std::string const& type, unsigned int size) const {
  // note: we need to implicitly convert Eigen types or they a returned as: Eigen::CwiseNullaryOp<.,.>
  if( eigenVec2Type.compare(type)==0 ) { // 2D Eigen vector
    return (Eigen::Vector2d)Eigen::Vector2d::Zero();
  } else if( eigenVec3Type.compare(type)==0 ) { // 3D Eigen vector
    return (Eigen::Vector3d)Eigen::Vector3d::Zero();
  } else if( eigenVec4Type.compare(type)==0 ) { // 4D Eigen vector
    return (Eigen::Vector4d)Eigen::Vector4d::Zero();
  } else if( eigenVecType.compare(type)==0 ) { // XD Eigen vector
    return (Eigen::VectorXd)Eigen::VectorXd::Zero(size);
  } else if( doubleType.compare(type)==0 ) { // double
    return 0.0;
  } 
  
  return ZeroVector(type, size);
}

boost::any AnyAlgebra::ZeroVector(std::string const& type, unsigned int const size) const {
  std::cerr << std::endl << "ERROR: cannot compute zero of a vector with type " << boost::core::demangle(type.c_str()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::ZeroVector()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::ZeroVector()" << std::endl << std::endl;
  assert(false);

  return boost::none;
}

double AnyAlgebra::NormBase(boost::any const& obj) const {
  const std::string type = obj.type().name();
    
  if( eigenVec2Type.compare(type)==0 ) { // 2D Eigen vector
    const Eigen::Vector2d& v = boost::any_cast<Eigen::Vector2d const&>(obj);
    return v.norm();
  } else if( eigenVec3Type.compare(type)==0 ) { // 3D Eigen vector
    const Eigen::Vector3d& v = boost::any_cast<Eigen::Vector3d const&>(obj);
    return v.norm();
  } else if( eigenVec4Type.compare(type)==0 ) { // 4D Eigen vector
    const Eigen::Vector4d& v = boost::any_cast<Eigen::Vector4d const&>(obj);
    return v.norm();
  } else if( eigenVecType.compare(type)==0 ) { // XD Eigen vector
    const Eigen::VectorXd& v = boost::any_cast<Eigen::VectorXd const&>(obj);
    return v.norm();
  } else if( eigenMatType.compare(type)==0 ) { // XD Eigen matrix
    const Eigen::MatrixXd& v = boost::any_cast<Eigen::MatrixXd const&>(obj);
    return v.norm();
  } else if( doubleType.compare(type)==0 ) { // double
    return std::abs(boost::any_cast<double const>(obj));
  }
  
  return 0.0;
}

double AnyAlgebra::Norm(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: Cannot compute the norm of an object with type " << boost::core::demangle(obj.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::Norm()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::Norm()" << std::endl << std::endl;
  assert(false);

  return -1.0;
}

bool AnyAlgebra::IsZeroBase(boost::any const& obj) const {
  const std::string type = obj.type().name();
  
  if( eigenVec2Type.compare(type)==0 ) { // 2D Eigen vector
    const Eigen::Vector2d& v = boost::any_cast<Eigen::Vector2d const&>(obj);
    return (v.array()==Eigen::Array2d::Zero()).all();
  } else if( eigenVec3Type.compare(type)==0 ) { // 3D Eigen vector
    const Eigen::Vector3d& v = boost::any_cast<Eigen::Vector3d const&>(obj);
    return (v.array()==Eigen::Array3d::Zero()).all();
  } else if( eigenVec4Type.compare(type)==0 ) { // 4D Eigen vector
    const Eigen::Vector4d& v = boost::any_cast<Eigen::Vector4d const&>(obj);
    return (v.array()==Eigen::Array4d::Zero()).all();
  } else if( eigenVecType.compare(type)==0 ) { // XD Eigen vector
    const Eigen::VectorXd& v = boost::any_cast<Eigen::VectorXd const&>(obj);
    return (v.array()==Eigen::ArrayXd::Zero(v.size())).all();
  } else if( eigenMatType.compare(type)==0 ) { // XD Eigen matrix
    const Eigen::MatrixXd& v = boost::any_cast<Eigen::MatrixXd const&>(obj);
    return (v.array()==Eigen::ArrayXd::Zero(v.rows(), v.cols())).all();
  } else if( doubleType.compare(type)==0 ) { // double
    return boost::any_cast<double const>(obj)==0.0;
  }
  
  return IsZero(obj);
}

bool AnyAlgebra::IsZero(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: No way to determine if an object with type " << boost::core::demangle(obj.type().name()) << " is the zero vector." << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::IsZero()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::IsZero()" << std::endl << std::endl;
  assert(false);

  return false;
}

unsigned int AnyAlgebra::VectorDimension(boost::any const& vec) const {
  std::cerr << std::endl << "ERROR: No way to compute the dimension of a vector with type " << boost::core::demangle(vec.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::VectorDimension()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::VectorDimension()" << std::endl << std::endl;
  assert(false);
  
  return 0;
}

boost::any AnyAlgebra::AccessElementBase(unsigned int i, boost::any const& vec) const {
  if( doubleType.compare(vec.type().name())==0 ) { // the type is a double
    // check the size
    assert(i==0);

    // return the value (it is actually a scalar)
    return vec;
  }

  if( eigenVec2Type.compare(vec.type().name())==0 ) { // the type is an Eigen::Vector2d
    // check the size
    assert(i<2);

    // get a constant reference to the Eigen::Vector2d
    const Eigen::Vector2d& vecref = boost::any_cast<const Eigen::Vector2d&>(vec);

    // return ith element
    return vecref(i);
  }

  if( eigenVec3Type.compare(vec.type().name())==0 ) { // the type is an Eigen::Vector3d
    // check the size
    assert(i<3);

    // get a constant reference to the Eigen::Vector3d
    const Eigen::Vector3d& vecref = boost::any_cast<const Eigen::Vector3d&>(vec);

    // return ith element
    return vecref(i);
  }

  if( eigenVec4Type.compare(vec.type().name())==0 ) { // the type is an Eigen::Vector4d
    // check the size
    assert(i<4);

    // get a constant reference to the Eigen::Vector4d
    const Eigen::Vector4d& vecref = boost::any_cast<const Eigen::Vector4d&>(vec);

    // return ith element
    return vecref(i);
  }

  if( eigenVecType.compare(vec.type().name())==0 ) { // the type is an Eigen::VectorXd
    // get a constant reference to the Eigen::VectorXd
    const Eigen::VectorXd& vecref = boost::any_cast<const Eigen::VectorXd&>(vec);

    // check the size
    assert(i<vecref.size());
    
    // return the ith element
    return vecref(i);
  }

#if MUQ_HAS_SUNDIALS==1
  if( N_VectorType.compare(vec.type().name())==0 ) { // the type is a N_Vector
    // get a constant reference to the N_Vector
    const N_Vector& vecref = boost::any_cast<const N_Vector&>(vec);

    // check the size
    assert(i<NV_LENGTH_S(vecref));

    // return the ith element
    return NV_Ith_S(vecref, i);
  }
#endif
  
  return AccessElement(i, vec);
}

boost::any AnyAlgebra::AccessElement(unsigned int i, boost::any const& vec) const {
  std::cerr << std::endl << "ERROR: No way to access element " << i << " of a vector with type " << boost::core::demangle(vec.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::AccessElement()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::AccessElement()" << std::endl << std::endl;
  assert(false);

  return boost::none;
}

boost::any AnyAlgebra::IdentityBase(std::reference_wrapper<const boost::any> const& in) const {
  // Eigen::VectorXd type
  if( eigenVecType.compare(in.get().type().name())==0 ) {
    const Eigen::VectorXd& inVec = boost::any_cast<const Eigen::VectorXd&>(in);
    
    return (Eigen::MatrixXd)Eigen::MatrixXd::Identity(inVec.size(), inVec.size());
  }
  
  // double type
  if( doubleType.compare(in.get().type().name())==0 ) {
    return 1.0;
  }
  
  return Identity(in);
}

boost::any AnyAlgebra::Identity(std::reference_wrapper<const boost::any> const& in) const {
  std::cerr << std::endl << "ERROR: No way to compute identy object with type " << boost::core::demangle(in.get().type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::Identity()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::Identity()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::AddBase(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  // the first type is boost::none
  if( noneType.compare(in0.get().type().name())==0 ) {
    // return the second
    return in1.get();
  }

  // both in/out are double
  if( doubleType.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const double in0d = boost::any_cast<const double>(in0);
    const double in1d = boost::any_cast<const double>(in1);

    return in0d+in1d;
  }

  // both in/out are Eigen::MatrixXd
  if( eigenMatType.compare(in0.get().type().name())==0 && eigenMatType.compare(in1.get().type().name())==0 ) {
    const Eigen::MatrixXd& in0Mat = boost::any_cast<const Eigen::MatrixXd>(in0);
    const Eigen::MatrixXd& in1Mat = boost::any_cast<const Eigen::MatrixXd>(in1);

    // make sure the sizes match
    assert(in0Mat.rows()==in1Mat.rows());
    assert(in0Mat.cols()==in1Mat.cols());
	
    return (Eigen::MatrixXd)(in0Mat+in1Mat);
  }

  // both in/out are Eigen::VectorXd
  if( eigenVecType.compare(in0.get().type().name())==0 && eigenVecType.compare(in1.get().type().name())==0 ) {
    const Eigen::VectorXd& in0Vec = boost::any_cast<const Eigen::VectorXd>(in0);
    const Eigen::VectorXd& in1Vec = boost::any_cast<const Eigen::VectorXd>(in1);

    // make sure the sizes match
    assert(in0Vec.size()==in1Vec.size());
	
    return (Eigen::VectorXd)(in0Vec+in1Vec);
  }

  // both in/out are Eigen::Vector2d
  if( eigenVec2Type.compare(in0.get().type().name())==0 && eigenVec2Type.compare(in1.get().type().name())==0 ) {
    const Eigen::Vector2d& in0Vec = boost::any_cast<const Eigen::Vector2d>(in0);
    const Eigen::Vector2d& in1Vec = boost::any_cast<const Eigen::Vector2d>(in1);

    return (Eigen::Vector2d)(in0Vec+in1Vec);
  }

  // both in/out are Eigen::Vector3d
  if( eigenVec3Type.compare(in0.get().type().name())==0 && eigenVec3Type.compare(in1.get().type().name())==0 ) {
    const Eigen::Vector3d& in0Vec = boost::any_cast<const Eigen::Vector3d>(in0);
    const Eigen::Vector3d& in1Vec = boost::any_cast<const Eigen::Vector3d>(in1);

    return (Eigen::Vector3d)(in0Vec+in1Vec);
  }

  // both in/out are Eigen::Vector4d
  if( eigenVec4Type.compare(in0.get().type().name())==0 && eigenVec4Type.compare(in1.get().type().name())==0 ) {
    const Eigen::Vector4d& in0Vec = boost::any_cast<const Eigen::Vector4d>(in0);
    const Eigen::Vector4d& in1Vec = boost::any_cast<const Eigen::Vector4d>(in1);

    return (Eigen::Vector4d)(in0Vec+in1Vec);
  }

  return Add(in0, in1);
}

boost::any AnyAlgebra::Add(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  std::cerr << std::endl << "ERROR: No way to add type " << boost::core::demangle(in0.get().type().name()) << " and type " << boost::core::demangle(in1.get().type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::Add()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::Add()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::SubtractBase(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  // the first type is boost::none
  if( noneType.compare(in0.get().type().name())==0 ) {
    // return the second
    return in1.get();
  }

  // both in/out are double
  if( doubleType.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const double in0d = boost::any_cast<const double>(in0);
    const double in1d = boost::any_cast<const double>(in1);

    return in0d-in1d;
  }

  // both in/out are Eigen::MatrixXd
  if( eigenMatType.compare(in0.get().type().name())==0 && eigenMatType.compare(in1.get().type().name())==0 ) {
    const Eigen::MatrixXd& in0Mat = boost::any_cast<const Eigen::MatrixXd>(in0);
    const Eigen::MatrixXd& in1Mat = boost::any_cast<const Eigen::MatrixXd>(in1);

    // make sure the sizes match
    assert(in0Mat.rows()==in1Mat.rows());
    assert(in0Mat.cols()==in1Mat.cols());
	
    return (Eigen::MatrixXd)(in0Mat-in1Mat);
  }

  // both in/out are Eigen::VectorXd
  if( eigenVecType.compare(in0.get().type().name())==0 && eigenVecType.compare(in1.get().type().name())==0 ) {
    const Eigen::VectorXd& in0Vec = boost::any_cast<const Eigen::VectorXd>(in0);
    const Eigen::VectorXd& in1Vec = boost::any_cast<const Eigen::VectorXd>(in1);

    // make sure the sizes match
    assert(in0Vec.size()==in1Vec.size());
	
    return (Eigen::VectorXd)(in0Vec-in1Vec);
  }

  // both in/out are Eigen::Vector2d
  if( eigenVec2Type.compare(in0.get().type().name())==0 && eigenVec2Type.compare(in1.get().type().name())==0 ) {
    const Eigen::Vector2d& in0Vec = boost::any_cast<const Eigen::Vector2d>(in0);
    const Eigen::Vector2d& in1Vec = boost::any_cast<const Eigen::Vector2d>(in1);

    return (Eigen::Vector2d)(in0Vec-in1Vec);
  }

  // both in/out are Eigen::Vector3d
  if( eigenVec3Type.compare(in0.get().type().name())==0 && eigenVec3Type.compare(in1.get().type().name())==0 ) {
    const Eigen::Vector3d& in0Vec = boost::any_cast<const Eigen::Vector3d>(in0);
    const Eigen::Vector3d& in1Vec = boost::any_cast<const Eigen::Vector3d>(in1);

    return (Eigen::Vector3d)(in0Vec-in1Vec);
  }

  // both in/out are Eigen::Vector4d
  if( eigenVec4Type.compare(in0.get().type().name())==0 && eigenVec4Type.compare(in1.get().type().name())==0 ) {
    const Eigen::Vector4d& in0Vec = boost::any_cast<const Eigen::Vector4d>(in0);
    const Eigen::Vector4d& in1Vec = boost::any_cast<const Eigen::Vector4d>(in1);

    return (Eigen::Vector4d)(in0Vec-in1Vec);
  }

  return Subtract(in0, in1);
}

boost::any AnyAlgebra::Subtract(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  std::cerr << std::endl << "ERROR: No way to subtract type " << boost::core::demangle(in0.get().type().name()) << " and type " << boost::core::demangle(in1.get().type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::Subtract()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::Subtract()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::MultiplyBase(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  // both in/out are Eigen::MatrixXd
  if( eigenMatType.compare(in0.get().type().name())==0 && eigenMatType.compare(in1.get().type().name())==0 ) {
    const Eigen::MatrixXd& in0Mat = boost::any_cast<const Eigen::MatrixXd>(in0);
    const Eigen::MatrixXd& in1Mat = boost::any_cast<const Eigen::MatrixXd>(in1);

    // make sure the sizes match
    assert(in0Mat.cols()==in1Mat.rows());
	
    return (Eigen::MatrixXd)(in0Mat*in1Mat);
  }

  // double times Eigen::MatrixXd
  if( doubleType.compare(in0.get().type().name())==0 && eigenMatType.compare(in1.get().type().name())==0 ) {
    const double in0Mat = boost::any_cast<double>(in0);
    const Eigen::MatrixXd& in1Mat = boost::any_cast<const Eigen::MatrixXd>(in1);

    return (Eigen::MatrixXd)(in0Mat*in1Mat);
  }

  // Eigen::MatrixXd times double 
  if( eigenMatType.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const Eigen::MatrixXd& in0Mat = boost::any_cast<const Eigen::MatrixXd>(in1);
    const double in1Mat = boost::any_cast<double>(in0);
    
    return (Eigen::MatrixXd)(in0Mat*in1Mat);
  }

  // Eigen::Vector2d times double 
  if( eigenVec2Type.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const Eigen::Vector2d& in0Vec = boost::any_cast<const Eigen::Vector2d>(in0);
    const double in1 = boost::any_cast<double>(in1);
	
    return (Eigen::Vector2d)(in0Vec*in1);
  }

  // double times Eigen::Vector2d
  if( eigenVec2Type.compare(in1.get().type().name())==0 && doubleType.compare(in0.get().type().name())==0 ) {
    const Eigen::Vector2d& in0Vec = boost::any_cast<const Eigen::Vector2d>(in1);
    const double in1 = boost::any_cast<double>(in0);
	
    return (Eigen::Vector2d)(in0Vec*in1);
  }

  // Eigen::Vector3d times double 
  if( eigenVec3Type.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const Eigen::Vector3d& in0Vec = boost::any_cast<const Eigen::Vector3d>(in0);
    const double in1 = boost::any_cast<double>(in1);
	
    return (Eigen::Vector3d)(in0Vec*in1);
  }

  // double times Eigen::Vector3d
  if( eigenVec3Type.compare(in1.get().type().name())==0 && doubleType.compare(in0.get().type().name())==0 ) {
    const Eigen::Vector3d& in0Vec = boost::any_cast<const Eigen::Vector3d>(in1);
    const double in1 = boost::any_cast<double>(in0);
	
    return (Eigen::Vector3d)(in0Vec*in1);
  }

  // Eigen::Vector4d times double 
  if( eigenVec4Type.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const Eigen::Vector4d& in0Vec = boost::any_cast<const Eigen::Vector4d>(in0);
    const double in1 = boost::any_cast<double>(in1);
	
    return (Eigen::Vector4d)(in0Vec*in1);
  }

  // double times Eigen::Vector4d
  if( eigenVec4Type.compare(in1.get().type().name())==0 && doubleType.compare(in0.get().type().name())==0 ) {
    const Eigen::Vector4d& in0Vec = boost::any_cast<const Eigen::Vector4d>(in1);
    const double in1 = boost::any_cast<double>(in0);
	
    return (Eigen::Vector4d)(in0Vec*in1);
  }

  // Eigen::VectorXd times double 
  if( eigenVecType.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const Eigen::VectorXd& in0Vec = boost::any_cast<const Eigen::VectorXd>(in0);
    const double in1 = boost::any_cast<double>(in1);
	
    return (Eigen::Vector4d)(in0Vec*in1);
  }

  // double times Eigen::Vector4d
  if( eigenVecType.compare(in1.get().type().name())==0 && doubleType.compare(in0.get().type().name())==0 ) {
    const Eigen::VectorXd& in0Vec = boost::any_cast<const Eigen::VectorXd>(in1);
    const double in1 = boost::any_cast<double>(in0);
	
    return (Eigen::VectorXd)(in0Vec*in1);
  }

  return Multiply(in0, in1);
}

boost::any AnyAlgebra::Multiply(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  std::cerr << std::endl << "ERROR: No way to multiply type " << boost::core::demangle(in0.get().type().name()) << " and type " << boost::core::demangle(in1.get().type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::Multiply()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::Multiply()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}
