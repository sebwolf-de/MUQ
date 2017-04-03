# make sure that the boost graph library is available
set(CMAKE_REQUIRED_LIBRARIES ${BOOST_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
"
#include <boost/graph/adjacency_list.hpp>
typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::bidirectionalS, int, int> Graph;
int main(){
Graph temp;
return 0; 
}
"
BOOST_GRAPH_COMPILES)

# make sure boost serialization works
set(CMAKE_REQUIRED_LIBRARIES ${BOOST_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIR})
CHECK_CXX_SOURCE_COMPILES(
"
#include <iostream>
#include <fstream>
#include <string>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/shared_ptr_helper.hpp>
int main(){
int temp = 10;
std::string filename;
{
std::ofstream ofs(filename.c_str());
boost::archive::text_oarchive oa(ofs);
oa << temp;
}
{
std::ifstream ifs(filename.c_str());
boost::archive::text_iarchive ia(ifs);
ia >> temp;
}
return 0; 
 }
"
BOOST_SERIALIZATION_OLD_COMPILES)

CHECK_CXX_SOURCE_COMPILES(
"
#include <iostream>
#include <fstream>
#include <string>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr_helper.hpp>
int main(){
int temp = 10;
std::string filename;
{
std::ofstream ofs(filename.c_str());
boost::archive::text_oarchive oa(ofs);
oa << temp;
}
{
std::ifstream ifs(filename.c_str());
boost::archive::text_iarchive ia(ifs);
ia >> temp;
}
return 0; 
 }
"
BOOST_SERIALIZATION_NEW_COMPILES)




CHECK_CXX_SOURCE_COMPILES(
"
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/filesystem.hpp>
using namespace std;
using namespace boost::filesystem;
int main(int argc, char* argv[])
{
path p (argv[1]);
directory_iterator temp(p);
directory_iterator temp2(temp);
return 0;
}
"
BOOST_DIRECTORY_ITERATOR_COMPILES)


if(MUQ_USE_PYTHON)
	set(CMAKE_REQUIRED_LIBRARIES ${BOOST_LIBRARIES})
	set(CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIR})
	CHECK_CXX_SOURCE_COMPILES(
	"
	#include <memory>
	#include <boost/get_pointer.hpp>
	int main(){
	std::shared_ptr<double> tempDbl = std::make_shared<double>(1.0);
	boost::get_pointer(tempDbl);
	return 0; 
	 }
	"
	BOOST_GETPOINTER_COMPILES)
else()
	set(BOOST_GETPOINTER_COMPILES 1)
endif()



if(NOT (BOOST_SERIALIZATION_OLD_COMPILES OR BOOST_SERIALIZATION_NEW_COMPILES) OR NOT BOOST_GRAPH_COMPILES OR NOT BOOST_GETPOINTER_COMPILES OR NOT BOOST_DIRECTORY_ITERATOR_COMPILES)
	set(BOOST_TEST_FAIL 1)
else()
	set(BOOST_TEST_FAIL 0)
endif()


if(NOT BOOST_TEST_FAIL)
    if(BOOST_SERIALIZATION_NEW_COMPILES)
        set(MUQ_USE_NEW_BOOST ON)
    else()
        set(MUQ_USE_NEW_BOOST OFF)
    endif()
endif()