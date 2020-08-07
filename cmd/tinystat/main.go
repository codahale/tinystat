// tinystat is used to compare two or more sets of measurements (e.g., runs of a
// multiple runs of benchmarks of two possible implementations) and determine if
// they are statistically different. It's inspired largely by FreeBSD's ministat
// (written by Poul-Henning Kamp).
package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"os"
	"path"
	"strconv"
	"strings"
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
	v          = flag.Bool("version", false, "display the version number")
	version    = "dev" // injected by goreleaser
)

func main() {
	flag.Usage = func() {
		_, _ = fmt.Fprintf(os.Stderr, "Usage: tinystat <control.csv> [<experiment.csv>...] [options] \n\n")

		flag.PrintDefaults()
	}
	flag.Parse()

	if *v {
		fmt.Println(version)

		return
	}

	args := flag.Args()
	if len(args) < 1 {
		flag.Usage()
		os.Exit(-1)
	}

	controlFilename := args[0]
	experimentFilenames := args[1:]

	// read the data
	controlData, experimentData, err := readData(controlFilename, experimentFilenames)
	if err != nil {
		_, _ = fmt.Fprintln(os.Stderr, err)
		os.Exit(-1)
	}

	// chart the data
	if !*noChart {
		printChart(experimentFilenames, controlFilename, controlData, experimentData)
	}

	// compare the data
	if len(experimentFilenames) > 0 {
		printComparison(controlFilename, controlData, experimentFilenames, experimentData)
	}
}

func printComparison(
	controlFilename string, controlData []float64,
	experimentFilenames []string, experimentData map[string][]float64,
) {
	t := tabwriter.NewWriter(os.Stdout, 2, 0, 2, ' ', 0)
	_, _ = fmt.Fprintf(t, "File\tN\tMean\tStddev\t\n")

	control := tinystat.Summarize(controlData)
	_, _ = fmt.Fprintf(t, "%s\t%.0f\t%.2f\t%0.2f\t%s\n", path.Base(controlFilename),
		control.N, control.Mean, control.StdDev(), "(control)")

	for _, filename := range experimentFilenames {
		experiment := tinystat.Summarize(experimentData[filename])
		d := tinystat.Compare(control, experiment, *confidence)
		p := strings.TrimLeft(fmt.Sprintf("%.3f", d.PValue), "0")

		var results string

		if d.Significant() {
			operator := ">"
			if experiment.Mean < control.Mean {
				operator = "<"
			}

			results = fmt.Sprintf("(%.2f %s %.2f ± %.2f, p = %s)",
				experiment.Mean, operator, control.Mean, d.CriticalValue, p)
		} else {
			results = fmt.Sprintf("(no difference, p = %s)", p)
		}

		_, _ = fmt.Fprintf(t, "%s\t%.0f\t%.2f\t%0.2f\t%s\n",
			path.Base(filename), experiment.N, experiment.Mean, experiment.StdDev(), results)
	}

	_ = t.Flush()
}

func readData(controlFilename string, experimentFilenames []string) ([]float64, map[string][]float64, error) {
	controlData, err := readFile(controlFilename, *column, *delimiter)
	if err != nil {
		return nil, nil, err
	}

	experimentData := make(map[string][]float64)

	for _, filename := range experimentFilenames {
		expData, err := readFile(filename, *column, *delimiter)
		if err != nil {
			return nil, nil, err
		}

		experimentData[filename] = expData
	}

	return controlData, experimentData, nil
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
