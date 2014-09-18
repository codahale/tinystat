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
//     Experiment  Results
//     chameleon   No difference proven at 95% confidence.
//     leopard     Difference at 95% confidence!
//                   343.5 +/- 292.63453863922877
//                   114.5% +/- 97.54484621307626%
//                   (Student's t, pooled s = 238.9799344943192)
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
	"path"
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

	// read the data
	controlData, err := readFile(controlFilename, column, delimiter)
	if err != nil {
		panic(err)
	}

	experimentData := make(map[string][]float64)
	for _, filename := range experimentFilenames {
		expData, err := readFile(filename, column, delimiter)
		if err != nil {
			panic(err)
		}
		experimentData[filename] = expData
	}

	// chart the data
	if args["--no-chart"] != true {
		c := chart.BoxChart{}
		c.XRange.Fixed(-1, float64(len(experimentFilenames))+1, 1)
		c.XRange.Category = make([]string, len(experimentFilenames)+1)
		c.XRange.Category[0] = path.Base(controlFilename)
		for i, file := range experimentFilenames {
			c.XRange.Category[i+1] = path.Base(file)
		}

		c.AddSet(0, controlData, true)
		for i, filename := range experimentFilenames {
			c.AddSet(float64(i+1), experimentData[filename], true)
		}

		txt := txtg.New(width, height)
		c.Plot(txt)
		fmt.Println(txt)
	}

	// compare the data
	if len(experimentFilenames) > 0 {
		control := tinystat.Summarize(controlData)
		table := new(tabwriter.Writer)
		table.Init(os.Stdout, 2, 0, 2, ' ', 0)
		fmt.Fprintln(table, "Experiment\tResults\t")
		for _, filename := range experimentFilenames {
			experimental := tinystat.Summarize(experimentData[filename])
			d := tinystat.Compare(control, experimental, confidence)

			if d.Significant() {
				fmt.Fprintf(table,
					"%s\tDifference at %v%% confidence!\t\n",
					filename, confidence)
				fmt.Fprintf(table,
					"\t  %v +/- %v\t\n",
					d.Delta, d.Error)
				fmt.Fprintf(table,
					"\t  %v%% +/- %v%%\t\n",
					d.RelDelta*100, d.RelError*100)
				fmt.Fprintf(table,
					"\t  (Student's t, pooled s = %v)\t\n",
					d.StdDev)
			} else {
				fmt.Fprintf(table,
					"%s\tNo difference proven at %v%% confidence.\t\n",
					filename, confidence)
			}
		}
		_ = table.Flush()
	}
}

func readFile(filename string, col int, del rune) ([]float64, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.Comma = del
	records, err := r.ReadAll()
	if err != nil {
		return nil, err
	}

	data := make([]float64, len(records))
	for i, s := range records {
		n, err := strconv.ParseFloat(s[col], 64)
		if err != nil {
			return nil, err
		}
		data[i] = n
	}

	return data, nil
}
