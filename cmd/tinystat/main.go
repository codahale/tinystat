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
	"flag"
	"fmt"
	"os"
	"path"
	"strconv"
	"text/tabwriter"

	"github.com/codahale/tinystat"
	"github.com/vdobler/chart"
	"github.com/vdobler/chart/txtg"
)

//nolint:gochecknoglobals // main package
var (
	confidence = flag.Float64("confidence", 95, "confidence level for analysis (0,100)")
	column     = flag.Int("column", 0, "the CSV column to analyze")
	delimiter  = flag.String("delimiter", ",", "the CSV delimiter to use")
	noChart    = flag.Bool("no_chart", false, "don't display the box chart")
	width      = flag.Int("width", 74, "width of box chart, in chars")
	height     = flag.Int("height", 20, "height of box chart, in chars")
)

func main() {
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
	}

	controlFilename := flag.Arg(0)
	experimentFilenames := flag.Args()[1:]

	// read the data
	controlData, experimentData := readData(controlFilename, experimentFilenames)

	// chart the data
	if !*noChart {
		printChart(experimentFilenames, controlFilename, controlData, experimentData)
	}

	// compare the data
	if len(experimentFilenames) > 0 {
		control := tinystat.Summarize(controlData)
		table := new(tabwriter.Writer)
		table.Init(os.Stdout, 2, 0, 2, ' ', 0)
		_, _ = fmt.Fprintln(table, "Experiment\tResults\t")

		for _, filename := range experimentFilenames {
			experimental := tinystat.Summarize(experimentData[filename])
			d := tinystat.Compare(control, experimental, *confidence)

			if d.Significant() {
				_, _ = fmt.Fprintf(table,
					"%s\tDifference at %v%% confidence!\t\n",
					filename, *confidence)
				_, _ = fmt.Fprintf(table,
					"\t  %v +/- %v\t\n",
					d.Delta, d.Error)
				_, _ = fmt.Fprintf(table,
					"\t  %v%% +/- %v%%\t\n",
					d.RelDelta*100, d.RelError*100)
				_, _ = fmt.Fprintf(table,
					"\t  (Student's t, pooled s = %v)\t\n",
					d.StdDev)
			} else {
				_, _ = fmt.Fprintf(table,
					"%s\tNo difference proven at %v%% confidence.\t\n",
					filename, *confidence)
			}
		}

		_ = table.Flush()
	}
}

func readData(controlFilename string, experimentFilenames []string) ([]float64, map[string][]float64) {
	controlData, err := readFile(controlFilename, *column, *delimiter)
	if err != nil {
		panic(err)
	}

	experimentData := make(map[string][]float64)

	for _, filename := range experimentFilenames {
		expData, err := readFile(filename, *column, *delimiter)
		if err != nil {
			panic(err)
		}

		experimentData[filename] = expData
	}

	return controlData, experimentData
}

func printChart(
	experimentFilenames []string, controlFilename string, controlData []float64,
	experimentData map[string][]float64,
) {
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

	txt := txtg.New(*width, *height)
	c.Plot(txt)
	fmt.Println(txt)
}

func readFile(filename string, col int, del string) ([]float64, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	defer func() { _ = f.Close() }()

	r := csv.NewReader(f)
	r.Comma = []rune(del)[0]

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
