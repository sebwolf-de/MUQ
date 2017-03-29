#include "MUQ/Modeling/WorkPiece.h"

/// An object used to test using user-created inputs to WorkPieces
struct AnObject {
  /// Construct the object
  /**
     @param[in] b The value of the object
   */
  AnObject(double const b) : value(b) {}
  
  /// The object's value
  const double value;
  
  /// A flag that changes the behavior of the object
  bool flag;
};

/// A WorkPiece with no fixed input/output number or type
class UnfixedMod : public muq::Modeling::WorkPiece {
public:

  /// Default constructor
  UnfixedMod() : WorkPiece() {}

  /// Default destructor
  virtual ~UnfixedMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the number of inputs
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
         
    switch( inputs.size() ) {
    case 0 : { // there are no inputs
      const std::string hi = "hello!";

      // the output size is 3 
      outputs.resize(3);

      outputs[0] = hi;
      outputs[1] = 3.0;
      outputs[2] = 6;

      return;
    } case 1: { // there is one input

      // the output size is 1
      outputs.resize(1);
      // first input and first ouptut must be a string but the parent does not check this
      outputs[0] = boost::any_cast<std::string>(inputs[0]);

      return;
      } default: // there is more than one input (there are no outputs)
      return;
    }
  }    
};

/// A WorkPiece with a fixed number of inputs but no fixed input types and no fixed output number or type
class FixedInsMod : public muq::Modeling::WorkPiece {
public:

  /// Constructor 
  /**
     @param[in] numIns The fixed number of inputs
   */
 FixedInsMod(unsigned int numIns) : WorkPiece(numIns) {}

  /// Default destructor
  virtual ~FixedInsMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the type of the third input
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    if( boost::any_cast<bool>(inputs[2]) ) { // if the third input is a bool
      // there are 2 outputs
      outputs.resize(2);

      // first input and first output must be a string but the parent does not check this
      outputs[0] = boost::any_cast<std::string>(inputs[0]);
      // second input and second output must be an unsigned int but the parent does not check this
      outputs[1] = boost::any_cast<unsigned int>(inputs[1]);
    } else { // if the third input is not a bool
      // there is 1 output
      outputs.resize(1);

      // first output and second input must be an unsigned int but the parent does not check this
      outputs[0] = boost::any_cast<unsigned int>(inputs[1]);
    }
  }    
};

/// A WorkPiece with no fixed input number or type and a fixed number of ouputs but no fixed output type
class FixedOutsMod : public muq::Modeling::WorkPiece {
public:

  /// Constructor
  /**
     @param[in] numOuts The fixed number of outputs
   */
  FixedOutsMod(unsigned int numOuts) : WorkPiece(numOuts, WorkPiece::Fix::Outputs) {}

  /// Default destructor
  virtual ~FixedOutsMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the number of inputs
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    outputs.resize(numOutputs);
    
    outputs[0] = 1.0;
    // if there is at least 1 input (must be an int but the parent does not check this)
    outputs[1] = inputs.size()>0? boost::any_cast<int>(inputs[0]) : 2;
  }    
};

/// A WorkPiece with a fixed number of inputs and outputs but no fixed types
class FixedInOutMod : public muq::Modeling::WorkPiece {
public:

  /// Constructor
  /**
     @param[in] numIns The fixed number of inputs
     @param[in] numOuts The fixed number of outputs
   */
  FixedInOutMod(unsigned int const numIns, unsigned int numOuts) : WorkPiece(numIns, numOuts) {}

  /// Default destructor
  virtual ~FixedInOutMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
   We don't actually do anything.
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    outputs = std::vector<boost::any>();
  }
};

/// A WorkPiece with a fixed input number and type but niether the number nor type of the output is fixed
class FixedInTypeMod : public muq::Modeling::WorkPiece {
public:

  /// Constructor
  /**
     @param[in] types A list of the input types (the fixed number of inputs is determined by the size of this vector)
   */
  FixedInTypeMod(std::vector<std::string> const& types) : WorkPiece(types) {}

  /// Default destructor
  virtual ~FixedInTypeMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the flag value of the input AnObject
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    // the first input must be a string and the parent checks this
    const std::string s = boost::any_cast<std::string>(inputs[0]);
    // the second input must be a shared pointer to AnObject and the parent checks this
    auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[2]);

    if( obj->flag ) { // of the flag in AnObject is true ...
      // ... the output size is 2
      outputs.resize(2);
      outputs[0] = s;
      // the second output must a double but the parent does not check this
      outputs[1] = obj->value*boost::any_cast<double>(inputs[1]);

      return;
    }

    // if the flag in AnObject is false, the ouput size is 1
    outputs.resize(1);
    outputs[0] = s;
  }
};

/// A WorkPiece with a fixed input number and type and a fixed number of outputs
class FixedInTypeOutNumMod : public muq::Modeling::WorkPiece {
public:

  /// Constructor
  /**
     @param[in] types A list of the input types (the fixed number of inputs is determined by the size of this vector)
     @param[in] numOuts The number of outputs
   */
 FixedInTypeOutNumMod(std::vector<std::string> const& types, unsigned int const numOuts) : WorkPiece(types, numOuts) {}

  /// Default destructor
  virtual ~FixedInTypeOutNumMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the flag value of the input AnObject
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    // the first input must be a string and the parent checks this
    const std::string s = boost::any_cast<std::string>(inputs[0]);
    // the second input must be a shared pointer to AnObject and the parent checks this
    auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[2]);

    if( obj->flag ) { // of the flag in AnObject is true ...
      // ... the output size is 2
      outputs.resize(2);
      outputs[0] = s;
      // the second output must a double but the parent does not check this
      outputs[1] = obj->value*boost::any_cast<double>(inputs[1]);

      return;
    }

    std::string anotherString = "second string";

    // if the flag in AnObject is false, the ouput size is still 2
    outputs.resize(2);
    outputs[0] = s;
    outputs[1] = anotherString;
  }
};

/// A WorkPiece with a fixed output number and type but niether the number nor type of the input is fixed
class FixedOutTypeMod : public muq::Modeling::WorkPiece {
public:

  /// Constructor
  /**
     @param[in] types A list of the output types (the fixed number of outputs is determined by the size of this vector)
   */
 FixedOutTypeMod(std::vector<std::string> const& types) : WorkPiece(types, WorkPiece::Fix::Outputs) {}

  /// Default destructor
  virtual ~FixedOutTypeMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the flag value of the input AnObject
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    // the first input must be a string but the parent does not check this
    const std::string s = boost::any_cast<std::string>(inputs[0]);
    // the second input must be a shared pointer to an object but the parent does not check this
    auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[1]);

    // there are 2 outputs
    outputs.resize(2);

    // the first output is a string and the parent checks this
    outputs[0] = s;
    
    // the second output is a double and the parent checks this --- but the value depends on the flag in AnObject
    if( obj->flag ) { // if the flag in AnObject is true
      outputs[1] = obj->value*boost::any_cast<double>(inputs[2]);
    } else { // if the flag in AnObject is false
      outputs[1] = (double)obj->value;
    }
  }
};

/// A WorkPiece with a fixed output number and type and a fixed number of inputs
class FixedOutTypeInNumMod : public muq::Modeling::WorkPiece {
public:

  /// Constructor
  /**
     @param[in] types A list of the output types (the fixed number of outputs is determined by the size of this vector)
     @param[in] numOuts The number of inputs
   */
 FixedOutTypeInNumMod(std::vector<std::string> const& types, unsigned int const numIns) : WorkPiece(types, numIns, WorkPiece::Fix::Outputs) {}

  /// Default destructor
  virtual ~FixedOutTypeInNumMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the flag value of the input AnObject
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    // the first input must be a string but the parent does not check this
    const std::string s = boost::any_cast<std::string>(inputs[0]);
    // the second input must be a shared pointer to AnObject but the parent does not check this
    auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[1]);

    if( obj->flag ) { // of the flag in AnObject is true ...
      // ... the output size is 2
      outputs.resize(2);
      outputs[0] = s;
      // the second output must a double and the parent checks this
      outputs[1] = obj->value*boost::any_cast<double>(inputs[2]);

      return;
    }

    // if the flag in AnObject is false, the ouput size is still 2
    outputs.resize(2);
    outputs[0] = s;
    outputs[1] = obj->value;
  }
};

/// A WorkPiece with a fixed number and type for both the inputs and the outputs
class FixedTypesMod : public muq::Modeling::WorkPiece {
public:

  /// Constructor
  /**
     @param[in] inTypes A list of the input types (the fixed number of inputs is determined by the size of this vector)
     @param[in] outTypes A list of the output types (the fixed number of outputs is determined by the size of this vector)
   */
 FixedTypesMod(std::vector<std::string> const& inTypes, std::vector<std::string> const& outTypes) : WorkPiece(inTypes, outTypes) {}

  /// Default destructor
  virtual ~FixedTypesMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the flag value of the input AnObject
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    // the first input must be a string but the parent does not check this
    const std::string s = boost::any_cast<std::string>(inputs[0]);
    // the second input must be a shared pointer to AnObject but the parent does not check this
    auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[1]);

    if( obj->flag ) { // of the flag in AnObject is true ...
      // ... the output size is 2
      outputs.resize(2);
      outputs[0] = s;
      // the second output must a double and the parent checks this
      outputs[1] = obj->value*boost::any_cast<double>(inputs[2]);

      return;
    }

    // if the flag in AnObject is false, the ouput size is still 2
    outputs.resize(2);
    outputs[0] = s;
    outputs[1] = obj->value;
  }
};
