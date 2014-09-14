tinystat
========

[![Build Status](https://travis-ci.org/codahale/tinystat.png?branch=master)](https://travis-ci.org/codahale/tinystat)

A Go library and command for evaluating whether two or more sets of measurements
are statistically different. It does this by performing a *Student's t-test* at
a particular confidence level, making it suitable for small sets of measurements
(e.g., multiple runs of a benchmark). It's inspired largely by FreeBSD's
`ministat` (written by Poul-Henning Kamp).

## Installation

    $ go get github.com/codahale/tinystat/...

## Documentation

For documentation, check the
[godoc for the Go package](http://godoc.org/github.com/codahale/tinystat) or the
[godoc for the `tinystat` command](http://godoc.org/github.com/codahale/tinystat/cmd/tinystat).
