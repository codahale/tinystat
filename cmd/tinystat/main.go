// tinystat is used to compare two or more sets of measurements (e.g., runs of a multiple runs of
// benchmarks of two possible implementations) and determine if they are statistically different.
// It's inspired largely by FreeBSD's ministat (written by Poul-Henning Kamp).
package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"path"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/alecthomas/kong"
	"github.com/codahale/tinystat"
	"github.com/vdobler/chart"
	"github.com/vdobler/chart/txtg"
)

var version = "dev"

func main() {
	//nolint:maligned // ordering of fields matters
	var cli struct {
		//nolint:lll // can't format struct field tags
		Confidence      float64          `short:"C" default:"95" help:"Confidence level for statistical significance (0,100)."`
		Column          int              `short:"c" default:"0" help:"The CSV column to analyze."`
		Delimiter       string           `short:"d" default:"," help:"The CSV delimiter to use."`
		NoChart         bool             `default:"false" help:"Don't display the box chart.'"`
		Width           int              `default:"74" help:"The width of the box chart in chars."`
		Height          int              `default:"20" help:"The height of the box chart in chars."`
		Version         kong.VersionFlag `help:"Display the application version."`
		ControlPath     string           `arg:"" type:"existingfile" help:"The CSV file containing measurements of the control group."`            //nolint:lll // can't format struct field tags
		ExperimentPaths []string         `arg:"" optional:"" type:"existingfile" help:"CSV files containing measurements of experimental groups."` //nolint:lll // can't format struct field tags
	}

	ctx := kong.Parse(&cli, kong.Vars{"version": version})
	if ctx.Error != nil {
		_, _ = fmt.Fprintln(os.Stderr, ctx.Error)
		os.Exit(1)
	}

	// read the data
	controlData, experimentData, err := readData(cli.ControlPath, cli.ExperimentPaths, cli.Column, cli.Delimiter)
	if err != nil {
		_, _ = fmt.Fprintln(os.Stderr, err)
		os.Exit(-1)
	}

	// chart the data
	if !cli.NoChart {
		printChart(cli.ExperimentPaths, cli.ControlPath, controlData, experimentData, cli.Width, cli.Height)
	}

	// compare the data
	if len(cli.ExperimentPaths) > 0 {
		printComparison(cli.ControlPath, controlData, cli.ExperimentPaths, experimentData, cli.Confidence)
	}
}

func printComparison(
	controlFilename string, controlData []float64,
	experimentFilenames []string, experimentData map[string][]float64,
	confidence float64,
) {
	t := tabwriter.NewWriter(os.Stdout, 2, 0, 2, ' ', 0)
	_, _ = fmt.Fprintf(t, "File\tN\tMean\tStddev\t\n")

	control := tinystat.Summarize(controlData)
	_, _ = fmt.Fprintf(t, "%s\t%.0f\t%.2f\t%0.2f\t%s\n", path.Base(controlFilename),
		control.N, control.Mean, control.StdDev(), "(control)")

	for _, filename := range experimentFilenames {
		experiment := tinystat.Summarize(experimentData[filename])
		d := tinystat.Compare(control, experiment, confidence)
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

func readData(
	controlFilename string, experimentFilenames []string,
	column int, delimiter string,
) ([]float64, map[string][]float64, error) {
	controlData, err := readFile(controlFilename, column, delimiter)
	if err != nil {
		return nil, nil, err
	}

	experimentData := make(map[string][]float64)

	for _, filename := range experimentFilenames {
		expData, err := readFile(filename, column, delimiter)
		if err != nil {
			return nil, nil, err
		}

		experimentData[filename] = expData
	}

	return controlData, experimentData, nil
}

func printChart(
	experimentFilenames []string, controlFilename string, controlData []float64,
	experimentData map[string][]float64, width, height int,
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

	txt := txtg.New(width, height)
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
