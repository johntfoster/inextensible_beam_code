FROM andrewosh/binder-base

MAINTAINER John Foster <johntfosterjr@gmail.com>

USER root

RUN apt-get update && \
    apt-get install gcc

RUN conda install --quiet --yes \
          'cffi' \
          'numpy' \
          'scipy' \
          'matplotlib'

ADD residual.cpp residual.cpp
ADD bspline.hpp bspline.hpp

RUN c++ -g -o $HOME/notebooks/_residual.so residual.cpp -std=c++14 -fPIC -shared 
