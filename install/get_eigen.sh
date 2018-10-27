#!/bin/bash
wget https://bitbucket.org/eigen/eigen/get/3.3.2.tar.gz -O /tmp/eigen.tar.gz
tar -xvzf /tmp/eigen.tar.gz
mkdir -p eigen/build && tar -xvzf eigen.tar.gz -C eigen --strip-components 1
cd eigen/build && cmake .. && make install