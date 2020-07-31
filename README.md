tinystat
========

A Go library and CLI tool for evaluating whether two or more sets of measurements are statistically
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

File       N  Mean    Stddev
iguana     7  300.00  238.05  (control)
chameleon  5  540.00  299.08  (no difference, p = 0.050)
leopard    6  643.50  240.09  (643.50 > 300.00 ± 293.97, p = 0.050)
```

As you can see, despite the superficial differences between the iguana's scores and the chameleon's
scores, there is no statistically significant difference between the two at a 95% confidence level.
The leopard, on the other hand, has statistically significantly different scores.
