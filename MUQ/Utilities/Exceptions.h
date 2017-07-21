#ifndef MUQEXCEPTIONS_H
#define MUQEXCEPTIONS_H

#include <exception>

namespace muq
{

    /** @defgroup Exceptions
       
     */
    /** @ingroup Exceptions
        @class NotImplementedError
        @brief Class for virtual base functions that are not implemented.
        @details In general, it's best to implement abstract class interfaces with pure virtual functions.  However,
                 there are some situations where not all children of the base class will implement a function.  This
                 exception is meant to be used in such a case.  It should be raised in the base classes virtual 
                 function.  When children override this function, no exception will be thrown.
    */
    class NotImplementedError : public std::logic_error
    {
    public:
        NotImplementedError(std::string const& message) : std::logic_error(message){};

    };

    
};



#endif // #ifndef MUQEXCEPTIONS
