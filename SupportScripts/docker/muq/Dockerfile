ARG tag
FROM mparno/muq-build:$tag
MAINTAINER M. Parno

# Run with "docker build --build-arg tag=v0.3.1 -t mparno/muq:v0.3.1 ."
ARG tag

# Download, compile, and install MUQ
RUN git clone --depth=1 --branch $tag https://bitbucket.org/mituq/muq2.git && \
    cd muq2/build && \
    cmake -DMUQ_USE_GTEST=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DMUQ_USE_PYTHON=ON -DCMAKE_INSTALL_PREFIX=/home/muq-user/muq-install .. && \
    make install && \
    cd ../../ && cp -r muq2/examples /home/muq-user/muq-examples && \
    rm -r muq2

ENV PYTHONPATH "${PYTHONPATH}:/home/muq-user/muq-install/python"
ENV LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/home/muq-user/muq-install/muq_external/lib:/home/muq-user/muq-install/lib"

CMD ["bash"]
