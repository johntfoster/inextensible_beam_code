FROM andrewosh/binder-base

MAINTAINER John Foster <johntfosterjr@gmail.com>

USER root

ADD residual.cpp residual.cpp
ADD bspline.hpp bspline.hpp

RUN c++ -g -o $HOME/notebooks/_residual.so residual.cpp -std=c++14 -fPIC -shared 
