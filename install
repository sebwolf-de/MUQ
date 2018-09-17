echo -n "Do you want to build docmentaion (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
    cd build/
    make doc
    cd ..
fi

cd build/
make -j10 install
cd ..

echo -n "Do you want to run tests (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
    ./build/RunAllTests
fi
