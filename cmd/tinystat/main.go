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
	"fmt"
	"os"
	"strconv"
	"text/tabwriter"

	"github.com/codahale/tinystat"
	"github.com/docopt/docopt-go"
)

func main() {
	usage := `tinystat helps you conduct experiments.

Usage:
  tinystat [options] <control> [<experiment>...]

Arguments:
  <control>     Input file containing control measurements.
  <experiment>  Input file containing experimental measurements.

Options:
  -h --help               Display this usage information.
  -c --confidence=<pct>   Specify confidence level for analysis. [default: 95]
  -C --column=<num>       The column to analyze. [default: 0]
  -d --delimiter=(t|s|c)  Tab, space, or comma delimited data. [default: t]
`

	args, err := docopt.Parse(usage, nil, true, "", true)
	if err != nil {
		panic(err)
	}

	confidence, err := strconv.ParseFloat(args["--confidence"].(string), 64)
	if err != nil {
		panic(err)
	}

	column, err := strconv.Atoi(args["--column"].(string))
	if err != nil {
		panic(err)
	}

	var delimiter rune
	switch d := args["--delimiter"]; d {
	case "t":
		delimiter = '\t'
	case "c":
		delimiter = ','
	case "s":
		delimiter = ' '
	default:
		panic(fmt.Errorf("bad delimiter: %#v", d))
	}

	controlFilename := args["<control>"].(string)
	experimentFilenames := args["<experiment>"].([]string)

	w := new(tabwriter.Writer)
	w.Init(os.Stdout, 2, 0, 2, ' ', 0)
	fmt.Fprintln(w, "Filename\tN\tMin\tMax\tMedian\tMean\tStdDev\t")

	control, err := readFile(controlFilename, column, delimiter)
	if err != nil {
		panic(err)
	}
	fmt.Fprintf(w,
		"%s\t%.0f\t%f\t%f\t%f\t%f\t%f\t\n",
		controlFilename,
		control.N,
		control.Min,
		control.Max,
		control.Median,
		control.Mean,
		control.StdDev,
	)

	for _, filename := range experimentFilenames {
		experimental, err := readFile(filename, column, delimiter)
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

		d := tinystat.Compare(control, experimental, confidence)
		if d.Significant() {
			fmt.Fprintf(w, "\tDifference at %v%% confidence!\n", confidence)
			fmt.Fprintf(w, "\t\t%10f +/- %10f\n", d.Delta, d.Error)
			fmt.Fprintf(w, "\t\t%9f%% +/- %9f%%\n", d.PctDelta, d.PctError)
			fmt.Fprintf(w, "\t\t(Student's t, pooled s = %f)\n", d.PooledStdDev)
		} else {
			fmt.Fprintf(w, "\tNo difference proven at %v%% confidence.\n", confidence)
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
