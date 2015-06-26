libcloudph++ - a cloud (micro)physics library

Quick install guide:

...libcloudphxx$ mkdir build
...libcloudphxx$ cd build
...libcloudphxx/build$ cmake .. -DCMAKE_CXX_COMPILER=clang++  
...libcloudphxx/build$ cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_INSTALL_PREFIX=/usr 
...libcloudphxx/build$ make 
...libcloudphxx/build$ make test
...libcloudphxx/build$ sudo make install
