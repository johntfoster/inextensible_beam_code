FROM andrewosh/binder-base

MAINTAINER John Foster <johntfosterjr@gmail.com>

USER root

ADD residual.cpp residual.cpp
ADD bsplines.hpp bsplines.hpp

RUN c++ -g -o /usr/local/lib/_residual.so residual.cpp -std=c++14 -fPIC -shared
