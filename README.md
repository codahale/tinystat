tinystat
========

A Go library and command for evaluating whether two or more sets of measurements are statistically
different. It does this by performing a *Welch's t-test* at a particular confidence level, making
it suitable for small sets of measurements (e.g., multiple runs of a benchmark). It's inspired
largely by FreeBSD's `ministat` (written by Poul-Henning Kamp).

Imagine we have the results of different animals' SAT scores. Each animal took the SATs multiple
times, and we're assuming that differences between each animal's attempts are measurement error
(i.e., normally distributed). We can test for differences as follows:

```
$ tinystat iguana chameleon leopard

 1.5 k  +
        |
        |
        |
        |
        |
  1000  +                                             |
        |                              |              |
        |                        +-----------+  +-----------+
        |              |         |           |  |           |
        |              |         |           |  +-----*-----+
        |              |         |     *     |  |           |
   500  +              |         +-----------+  +-----------+
        |        +-----------+   |           |        |
        |        |     *     |   +-----------+
        |        +-----------+         |
        |        +-----------+         |
     0  +--------------|-----------------------------------------------
                    iguana         chameleon       leopard

Experiment                Results
../../examples/chameleon  No difference proven at 95% confidence, p = 0.0500
../../examples/leopard    Difference at 95% confidence, p = 0.0500
                            643.5000   >   300.0000 ± 293.9689
                            114.5000%  >    97.9896%
```

As you can see, despite the superficial differences between the iguana's scores and the chameleon's
scores, there is no statistically significant difference between the two at a 95% confidence level.
The leopard, on the other hand, has statistically significantly different scores.
