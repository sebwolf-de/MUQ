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

  return VectorDimension(vec);
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

  return Add(in0, in1);
}

boost::any AnyAlgebra::Add(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  std::cerr << std::endl << "ERROR: No way to add type " << boost::core::demangle(in0.get().type().name()) << " and type " << boost::core::demangle(in1.get().type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::Add()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::Add()" << std::endl << std::endl;
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

    // make sure the sizes match
    assert(in1Mat.rows()==1);
	
    return (Eigen::MatrixXd)(in0Mat*in1Mat);
  }

  // Eigen::MatrixXd times double 
  if( eigenMatType.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const Eigen::MatrixXd& in0Mat = boost::any_cast<const Eigen::MatrixXd>(in1);
    const double in1Mat = boost::any_cast<double>(in0);
    
    // make sure the sizes match
    assert(in0Mat.cols()==1);
	
    return (Eigen::MatrixXd)(in0Mat*in1Mat);
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
