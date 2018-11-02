#!/bin/bash
cd /tmp
wget https://dl.bintray.com/boostorg/release/1.68.0/source/boost_1_68_0.tar.bz2 -O boost.tar.bz2
mkdir -p boost && tar -xvf boost.tar.bz2 -C boost --strip-components 1
ls
cd boost/
./bootstrap --with-libraries=program_options,test
./b2 install