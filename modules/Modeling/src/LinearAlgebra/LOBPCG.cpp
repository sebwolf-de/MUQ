#include "MUQ/Modeling/LinearAlgebra/LOBPCG.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include <boost/property_tree/ptree.hpp>

using namespace muq::Modeling;
using namespace muq::Utilities;

LOBPCG::LOBPCG(int    numEigsIn,
               double tolIn,
               int    maxItsIn,
               bool   largestIn,
               int    verbosityIn) : numEigs(numEigsIn),
                                     tol(tolIn),
                                     maxIts(maxItsIn),
                                     largest(largestIn),
                                     verbosity(verbosityIn)
{
  assert(numEigs>0);
}


LOBPCG::LOBPCG(boost::property_tree::ptree const& opts) : LOBPCG(opts.get<int>("NumEigs"),
                                                                 opts.get("Tolerance",-1.0),
                                                                 opts.get("MaxIts", -1),
                                                                 opts.get("Largest",true),
                                                                 opts.get("Verbosity",0))
{}

Eigen::MatrixXd LOBPCG::Orthonormalizer::Compute(Eigen::MatrixXd const& V)
{
  Eigen::MatrixXd output(V);
  if(B!=nullptr){
    ComputeInPlace(output,B->Apply(V));
  }else{
    ComputeInPlace(output,V);
  }
  return output;
}

Eigen::MatrixXd LOBPCG::Orthonormalizer::Compute(Eigen::MatrixXd const& V, Eigen::MatrixXd const& BVin)
{
  Eigen::MatrixXd output(V);
  ComputeInPlace(output, BVin);
  return output;
}

void LOBPCG::Orthonormalizer::ComputeInPlace(Eigen::MatrixXd& V)
{
  if(B!=nullptr){
    ComputeInPlace(V, B->Apply(V));
  }else{
    ComputeInPlace(V,V);
  }
}

void LOBPCG::Orthonormalizer::ComputeInPlace(Eigen::MatrixXd& V, Eigen::MatrixXd const& BVin)
{
  vDim = V.cols();

  VBV_Chol = (V.transpose()*BVin).eval().llt().matrixL();
  V = VBV_Chol.triangularView<Eigen::Lower>().solve(V.transpose()).transpose();

  if(B!=nullptr){
    BV = VBV_Chol.triangularView<Eigen::Lower>().solve(BV.transpose()).transpose();
  }
}

Eigen::MatrixXd LOBPCG::Orthonormalizer::InverseVBV() const
{
  return VBV_Chol.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(vDim,vDim));
}


LOBPCG::Constraints::Constraints(std::shared_ptr<LinearOperator> const& B,
                                 Eigen::MatrixXd                 const& constMat) : Y(constMat)
{
  if(B!=nullptr){
    BY = B->Apply(constMat);
  }else{
    BY = constMat;
  }

  YBY_llt = (constMat.transpose() * BY).eval().llt();
}

void LOBPCG::Constraints::ApplyInPlace(Eigen::MatrixXd& x)
{
    Eigen::MatrixXd YBX = BY.transpose()*x;
    x -= Y*YBY_llt.solve(YBX);
}


LOBPCG& LOBPCG::compute(std::shared_ptr<LinearOperator> const& A,
                        Eigen::MatrixXd                 const& constMat,
                        std::shared_ptr<LinearOperator>        B,
                        std::shared_ptr<LinearOperator>        M)
{

  const int dim = A->rows();
  if(numEigs>int(std::ceil(0.2*dim))){
    std::cerr << "ERROR: LOBPCG can fail when more than 20\% of the eigenvalues are computed." << std::endl;
    assert(numEigs>int(std::ceil(0.2*dim)));
  }

  // Initial message
  if(verbosity>0){
    std::cout << "Solving ";
    if(B==nullptr){
      std::cout << "standard ";
    }else{
      std::cout << "generalized ";
    }
    std::cout << "eigenvalue problem." << std::endl;

    if(verbosity>1)
      std::cout << "  matrix size = " << A->rows() << " x " << A->cols() << std::endl;
  }

  // If this is the first call to compute, initialize everything
  if(eigVals.size()==0)
    reset(dim); // This will initialize the eigenvectors to random values

  Eigen::MatrixXd X = eigVecs;

  // To start with, all of the indices are active
  unsigned int numActive = numEigs;

  // Get ready to apply constraints (e.g., build Y, BY, )
  Constraints consts(B, constMat);
  if(consts.size()>0)
    consts.ApplyInPlace(X); // this will remove components of the eigenvectors that lie in the constraint subspace

  // Make the current vectors orthogonal wrt B
  Orthonormalizer bOrtho(B);
  bOrtho.ComputeInPlace(X);

  Eigen::MatrixXd BX;
  if(B != nullptr){
    BX = bOrtho.GetBV();
  }else{
    BX = X;
  }

  // Compute the initial Ritz vectors by solving the reduced eigenproblem.
  Eigen::MatrixXd AX = A->Apply(X);

  Eigen::MatrixXd XAX = X.transpose() * AX;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ritzSolver(XAX);
  Eigen::VectorXd subEigVals = ritzSolver.eigenvalues();
  Eigen::MatrixXd const& subEigVecs = ritzSolver.eigenvectors();

  X = (X*subEigVecs).eval();
  AX = (AX*subEigVecs).eval();
  if(B!=nullptr)
    BX = (BX*subEigVecs).eval();

  //
  Eigen::MatrixXd resids;
  Eigen::VectorXd residNorms;
  Eigen::MatrixXd AR, P, BP, AP;

  std::vector<int> sortInds(numEigs);
  for(int i=0; i<numEigs; ++i)
    sortInds[i] = i;

  for(unsigned it=0; it<maxIts; ++it){
    std::cout << "Iteration " << it << std::endl;
    std::cout << "  starting eigVals = " << subEigVals.transpose() << std::endl;

    // Compute the residuals
    resids = AX - BX*subEigVals.asDiagonal();

    // Compute the norm of the residuals
    residNorms = resids.colwise().norm();
    std::cout << "  residNorms = " << residNorms.transpose() << std::endl;
    std::cout << "  X=\n" << X << std::endl;

    // Figure out what vectors we're still working on by sorting
    std::sort(sortInds.begin(), sortInds.end(), [&residNorms](int i1, int i2) {return residNorms(i1) < residNorms(i2);});
    std::cout << "Sorted residual norm = \n";
    for(unsigned int k=0; k<numEigs; ++k){
      std::cout << residNorms(sortInds[k]) << "  ";
    }
    std::cout << std::endl;

    // Make residuals B-orthogonal to current values
    resids -= (X*(BX.transpose()*resids)).eval();
    std::cout << "Residual ortho=\n" << X.transpose()*resids << std::endl;

    // Apply constraints

    // Orthonormalize residuals wrt B-inner product
    bOrtho.ComputeInPlace(resids);
    Eigen::MatrixXd BR;
    if(B != nullptr){
      BR = bOrtho.GetBV();
    } else {
      BR = resids;
    }



    AR = A->Apply(resids);

    if(it!=0){
      // Orthonormalize P
      bOrtho.ComputeInPlace(P,BP);
      AP = bOrtho.VBV_Chol.triangularView<Eigen::Lower>().solve(AP.transpose()).transpose();
    }

    Eigen::MatrixXd XAR = X.transpose()*AR;
    Eigen::MatrixXd RAR = resids.transpose()*AR;
    Eigen::MatrixXd XAX = X.transpose()*AX;
    Eigen::MatrixXd XBX, RBR, XBR;

    XBX = BX.transpose() * X;
    RBR = BR.transpose() * resids;
    XBR = X.transpose() * BR;


    //gramXAR=full(blockVectorAX'*blockVectorR(:,activeMask));
    //gramRAR=full(blockVectorAR(:,activeMask)'*blockVectorR(:,activeMask));
    //gramRAR=(gramRAR'+gramRAR)*0.5;

    Eigen::MatrixXd gramA(3*numEigs,3*numEigs);
    Eigen::MatrixXd gramB(3*numEigs,3*numEigs);

    if(it!=0){
      gramA.resize(3*numEigs,3*numEigs);
      gramB.resize(3*numEigs,3*numEigs);

      Eigen::MatrixXd XAP = X.transpose() * AP;
      Eigen::MatrixXd RAP = resids.transpose()*AP;
      Eigen::MatrixXd PAP = P.transpose()*AP;
      Eigen::MatrixXd XBP = X.transpose()*BP;
      Eigen::MatrixXd RBP = resids.transpose()*BP;

      gramA.block(0,0,numEigs,numEigs) = subEigVals.asDiagonal();
      gramA.block(0,numEigs,numEigs,numEigs) = XAR;
      gramA.block(0,2*numEigs,numEigs,numEigs) = XAP;
      gramA.block(numEigs,0,numEigs,numEigs) = XAR.transpose();
      gramA.block(numEigs,numEigs,numEigs,numEigs) = RAR;
      gramA.block(numEigs,2*numEigs,numEigs,numEigs) = RAP;
      gramA.block(2*numEigs,0,numEigs,numEigs) = XAP.transpose();
      gramA.block(2*numEigs,numEigs,numEigs,numEigs) = RAP.transpose();
      gramA.block(2*numEigs,2*numEigs,numEigs,numEigs) = PAP;

      gramB.block(0,0,numEigs,numEigs) = Eigen::MatrixXd::Identity(numEigs,numEigs);
      gramB.block(0,numEigs,numEigs,numEigs) = XBR;
      gramB.block(0,2*numEigs,numEigs,numEigs) = XBP;
      gramB.block(numEigs,0,numEigs,numEigs) = XBR.transpose();
      gramB.block(numEigs,numEigs,numEigs,numEigs) = Eigen::MatrixXd::Identity(numEigs,numEigs);
      gramB.block(numEigs,2*numEigs,numEigs,numEigs) = RBP;
      gramB.block(2*numEigs,0,numEigs,numEigs) = XBP.transpose();
      gramB.block(2*numEigs,numEigs,numEigs,numEigs) = RBP.transpose();
      gramB.block(2*numEigs,2*numEigs,numEigs,numEigs) = Eigen::MatrixXd::Identity(numEigs,numEigs);

    }else{

      gramA.resize(2*numEigs,2*numEigs);
      gramB.resize(2*numEigs,2*numEigs);

      gramA.block(0,0,numEigs,numEigs) = subEigVals.asDiagonal();
      gramA.block(0,numEigs,numEigs,numEigs) = XAR;
      gramA.block(numEigs,0,numEigs,numEigs) = XAR.transpose();
      gramA.block(numEigs,numEigs,numEigs,numEigs) = RAR;

      gramB.block(0,0,numEigs,numEigs) = Eigen::MatrixXd::Identity(numEigs,numEigs);
      gramB.block(0,numEigs,numEigs,numEigs) = XBR;
      gramB.block(numEigs,0,numEigs,numEigs) = XBR.transpose();
      gramB.block(numEigs,numEigs,numEigs,numEigs) = Eigen::MatrixXd::Identity(numEigs,numEigs);
    }

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ritzSolver(gramA, gramB);
    subEigVals = ritzSolver.eigenvalues().tail(numEigs);
    assert(subEigVals(0)<subEigVals(1));
    Eigen::Ref<const Eigen::MatrixXd> subEigVecs = ritzSolver.eigenvectors().rightCols(numEigs);

    std::cout << "  subEigVals = " << subEigVals.transpose() << std::endl;


    Eigen::Ref<const Eigen::MatrixXd> eigX = subEigVecs.topRows(numEigs);
    Eigen::Ref<const Eigen::MatrixXd> eigR = subEigVecs.block(numEigs,0,numEigs,subEigVecs.cols());

    Eigen::MatrixXd pp  = resids * eigR;
    Eigen::MatrixXd app = AR*eigR;
    Eigen::MatrixXd bpp = BR*eigR;

    if(it!=0){

      Eigen::Ref<const Eigen::MatrixXd> eigP = subEigVecs.bottomRows(numEigs);

      pp  += P*eigP;
      app += AP*eigP;
      bpp += BP*eigP;

    }


    X = (X*eigX + pp).eval();
    AX = (AX*eigX + app).eval();

    if(B!=nullptr){
      BX = (BX*eigX + bpp).eval();
    }else{
      BX = X;
    }

    P = pp;
    AP = app;
    BP = bpp;
  }

  return *this;
}


LOBPCG& LOBPCG::compute(std::shared_ptr<LinearOperator> const& A,
                        std::shared_ptr<LinearOperator>        B,
                        std::shared_ptr<LinearOperator>        M)
{
  return compute(A,Eigen::MatrixXd(), B,M);
};

LOBPCG& LOBPCG::reset(int dim)
{
  if(tol<0.0)
    tol = dim * std::sqrt(std::numeric_limits<double>::epsilon());

  if(maxIts<0)
    maxIts = std::min(dim,20);

  eigVals.resize(numEigs);
  eigVecs = RandomGenerator::GetNormal(dim,numEigs);

  return *this;
}
