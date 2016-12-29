#include "MUQ/Modelling/WorkPiece.h"

/// An opject used to test passing pointers to user-created things
struct AnObject {
  /// Construct an object
  /**
     @param[in] b The value of the object
   */
  AnObject(double const b) : value(b) {}

  /// A flag that changes the object's behavior
  bool flag; 

  /// The value of the object
  const double value;
};

