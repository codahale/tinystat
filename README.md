tinystat
========

[![Build Status](https://travis-ci.org/codahale/tinystat.png?branch=master)](https://travis-ci.org/codahale/tinystat)

A Go library and command for evaluating whether two or more sets of measurements
are statistically different. It does this by performing a *Student's t-test* at
a particular confidence level, making it suitable for small sets of measurements
(e.g., multiple runs of a benchmark).

For documentation, check [godoc](http://godoc.org/github.com/codahale/tinystat).
