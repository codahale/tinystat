// tinystat is used to compare two or more sets of measurements (e.g., runs of a
// multiple runs of benchmarks of two possible implementations) and determine if
// they are statistically different. It's inspired largely by FreeBSD's ministat
// (written by Poul-Henning Kamp).
//
// Imagine we have the results of different animals' SAT scores. Each animal
// took the SATs multiple times, and we're assuming that differences between
// each animal's attempts are measurement error (i.e., normally distributed). We
// can test for differences as follows:
//
//     $ tinystat iguana chameleon leopard
//
//      1.5 k  +
//             |
//             |
//             |
//             |
//             |
//       1000  +                                             |
//             |                              |              |
//             |                        +-----------+  +-----------+
//             |              |         |           |  |           |
//             |              |         |           |  +-----*-----+
//             |              |         |     *     |  |           |
//        500  +              |         +-----------+  +-----------+
//             |        +-----------+   |           |        |
//             |        |     *     |   +-----------+
//             |        +-----------+         |
//             |        +-----------+         |
//          0  +--------------|-----------------------------------------------
//                         iguana         chameleon       leopard
//
//     Filename   N  Min         Max         Median      Mean        StdDev
//     iguana     7  50.000000   750.000000  200.000000  300.000000  238.047614
//     chameleon  5  150.000000  930.000000  500.000000  540.000000  299.081929
//                No difference proven at 95% confidence.
//     leopard    6  353.000000  1057.000000  619.000000  643.500000  240.093940
//                Difference at 95% confidence!
//                  343.500000 +/- 292.634539
//                  53.379953% +/- 45.475453%
//                  (Student's t, pooled s = 238.979934)
//
// As you can see, despite the superficial differences between the iguana's
// scores and the chameleon's scores, there is no statistically significant
// difference between the two at a 95% confidence level. The leopard, on the
// other hand, has statistically significantly different scores.
package main

import (
	"bytes"
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
	"text/tabwriter"

	"github.com/codahale/tinystat"
	"github.com/docopt/docopt-go"
	"github.com/vdobler/chart"
	"github.com/vdobler/chart/txtg"
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
  -n --no-chart           Don't display the box chart.
  -W --width=<chars>      Width of box chart, in characters. [default: 74]
  -H --height=<chars>     Height of box chart, in characters. [default: 20]
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

	width, err := strconv.Atoi(args["--width"].(string))
	if err != nil {
		panic(err)
	}

	height, err := strconv.Atoi(args["--height"].(string))
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

	c := chart.BoxChart{}
	c.XRange.Fixed(-1, 3, 1)
	c.XRange.Category = append([]string{controlFilename}, experimentFilenames...)

	table := bytes.NewBuffer(nil)
	w := new(tabwriter.Writer)
	w.Init(table, 2, 0, 2, ' ', 0)
	fmt.Fprintln(w, "Filename\tN\tMin\tMax\tMedian\tMean\tStdDev\t")

	control, controlData, err := readFile(controlFilename, column, delimiter)
	if err != nil {
		panic(err)
	}
	c.AddSet(0, controlData, true)
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

	for i, filename := range experimentFilenames {
		experimental, expData, err := readFile(filename, column, delimiter)
		if err != nil {
			panic(err)
		}

		c.AddSet(float64(i+1), expData, true)

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

	if args["--no-chart"] != true {
		txt := txtg.New(width, height)
		c.Plot(txt)
		fmt.Println(txt)
	}

	fmt.Println(table.String())
}

func readFile(filename string, col int, del rune) (tinystat.Summary, []float64, error) {
	f, err := os.Open(filename)
	if err != nil {
		return tinystat.Summary{}, nil, err
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.Comma = del
	records, err := r.ReadAll()
	if err != nil {
		return tinystat.Summary{}, nil, err
	}

	data := make([]float64, len(records))
	for i, s := range records {
		n, err := strconv.ParseFloat(s[col], 64)
		if err != nil {
			return tinystat.Summary{}, nil, err
		}
		data[i] = n
	}

	return tinystat.Summarize(data), data, nil
}
