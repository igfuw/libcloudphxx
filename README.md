libcloudph++ - a cloud (micro)physics library

First, that's all work in progress!


...libcloudphxx$ mkdir build
...libcloudphxx$ cd build
...libcloudphxx/build$ cmake    \
  -DCMAKE_BUILD_TYPE=Release    \
  -DCMAKE_CXX_COMPILER=clang++  \
  -DCMAKE_INSTALL_PREFIX=/usr   \
  ..
...libcloudphxx/build$ make 
...libcloudphxx/build$ make test
...libcloudphxx/build$ sudo make install

