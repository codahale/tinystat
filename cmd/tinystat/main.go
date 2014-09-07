// tinystat is used to compare two or more sets of measurements (e.g., runs of a
// multiple runs of benchmarks of two possible implementations) and determine if
// they are statistically different.
//
// Imagine we have the results of different animals' SAT scores. Each animal
// took the SATs multiple times, and we're assuming that differences between
// each animal's attempts are measurement error (i.e., normally distributed). We
// can test for differences as follows:
//
//     $ tinystat examples/iguana examples/chameleon examples/leopard
//
//     Filename            N  Min         Max         Median      Mean        StdDev
//     examples/iguana     7  50.000000   750.000000  200.000000  300.000000  238.047614
//     examples/chameleon  5  150.000000  930.000000  500.000000  540.000000  299.081929
//                         No difference proven at 95% confidence.
//     examples/leopard    6  700.000000  700.000000  700.000000  700.000000  0.000000
//                         Difference at 95% confidence!
//                            400.000000 +/- 215.283224
//                            57.142857% +/- 30.754746%
//                            (Student's t, pooled s = 175.809815)
//
// As you can see, despite the superficial differences between the iguana's
// scores and the chameleon's scores, there is no statistically significant
// difference between the two at a 95% confidence level. The leopard, on the
// other hand, has statistically significantly different scores.
package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"os"
	"strconv"
	"text/tabwriter"

	"github.com/codahale/tinystat"
)

func main() {
	var (
		conf  = flag.Float64("confidence", 95, "confidence level")
		col   = flag.Int("column", 0, "the column of data to analyze")
		delim = flag.String("delimiter", "tab", "the column delimiter (tab, space)")
	)
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "tinystat helps you conduct experiments.\n\n")
		fmt.Fprintln(os.Stderr, "Usage: tinystat [OPTIONS] <CONTROL> [EXP1 EXP2 EXP3...]")
		fmt.Fprintln(os.Stderr, "\nOptions:")
		flag.PrintDefaults()
	}
	flag.Parse()

	var delimiter rune
	switch *delim {
	case "tab":
		delimiter = '\t'
	case "space":
		delimiter = ' '
	default:
		panic(fmt.Sprintf("invalid delimiter: %v", *delim))
	}

	files := flag.Args()
	if len(files) < 1 {
		panic("needs more args")
	}

	w := new(tabwriter.Writer)
	w.Init(os.Stdout, 2, 0, 2, ' ', 0)
	fmt.Fprintln(w, "Filename\tN\tMin\tMax\tMedian\tMean\tStdDev\t")

	control, err := readFile(files[0], *col, delimiter)
	if err != nil {
		panic(err)
	}
	fmt.Fprintf(w,
		"%s\t%.0f\t%f\t%f\t%f\t%f\t%f\t\n",
		files[0],
		control.N,
		control.Min,
		control.Max,
		control.Median,
		control.Mean,
		control.StdDev,
	)

	for _, filename := range files[1:] {
		experimental, err := readFile(filename, *col, delimiter)
		if err != nil {
			panic(err)
		}

		fmt.Fprintf(w,
			"%s\t%.0f\t%f\t%f\t%f\t%f\t%f\t\n",
			filename,
			experimental.N,
			experimental.Min,
			experimental.Max,
			experimental.Median,
			experimental.Mean,
			experimental.StdDev,
		)

		d := tinystat.Compare(control, experimental, *conf)
		if d.Significant() {
			fmt.Fprintf(w, "\tDifference at %v%% confidence!\n", *conf)
			fmt.Fprintf(w, "\t\t%10f +/- %10f\n", d.Delta, d.Error)
			fmt.Fprintf(w, "\t\t%9f%% +/- %9f%%\n", d.PctDelta, d.PctError)
			fmt.Fprintf(w, "\t\t(Student's t, pooled s = %f)\n", d.PooledStdDev)
		} else {
			fmt.Fprintf(w, "\tNo difference proven at %v%% confidence.\n", *conf)
		}
	}
	_ = w.Flush()
}

func readFile(filename string, col int, del rune) (tinystat.Summary, error) {
	f, err := os.Open(filename)
	if err != nil {
		return tinystat.Summary{}, err
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.Comma = del
	records, err := r.ReadAll()
	if err != nil {
		return tinystat.Summary{}, err
	}

	data := make([]float64, len(records))
	for i, s := range records {
		n, err := strconv.ParseFloat(s[col], 64)
		if err != nil {
			return tinystat.Summary{}, err
		}
		data[i] = n
	}

	return tinystat.Summarize(data), nil
}
