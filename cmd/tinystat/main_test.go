package main

import (
	"io/ioutil"
	"os"
	"strings"
	"testing"
	"unicode"
)

func TestControlOnly(t *testing.T) {
	expected := `
   800  +
        |
        |                              |
        |                              |
        |                              |
   600  +                              |
        |                              |
        |                              |
        |                              |
   400  +-------------------------------------------------------------+
        |                                                             |
        |                              *                              |
        |                                                             |
   200  +-------------------------------------------------------------+
        +-------------------------------------------------------------+
        |                              |
        |                              |
     0  +--------------------------------------------------------------
                                    iguana

`
	actual := mainTest(t, "../../examples/iguana")
	if actual != expected {
		t.Errorf("Output was \n%q\n but expected \n%q", actual, expected)
	}
}

func TestOneExperiment(t *testing.T) {
	expected := `
  1000  +
        |
        |                                        |
        |                              +-------------------+
        |                              |                   |
        |                   |          |                   |
        |                   |          |                   |
        |                   |          |                   |
        |                   |          |         *         |
   500  +                   |          +-------------------+
        |                   |          |                   |
        |         +-------------------+|                   |
        |         |         *         ||                   |
        |         |                   |+-------------------+
        |         +-------------------+          |
        |         +-------------------+          |
        |                   |
     0  +-------------------|------------------------------------------
                         iguana              chameleon

Experiment                Results
../../examples/chameleon  No difference proven at 95% confidence.
`
	actual := mainTest(t, "../../examples/iguana", "../../examples/chameleon")
	if actual != expected {
		t.Errorf("Output was \n%q\n but expected \n%q", actual, expected)
	}
}

func TestTwoExperiments(t *testing.T) {
	expected := `
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
../../examples/chameleon  No difference proven at 95% confidence.
../../examples/leopard    Difference at 95% confidence!
                            343.5 +/- 292.63453863922877
                            114.5% +/- 97.54484621307626%
                            (Student's t, pooled s = 238.9799344943192)
`
	actual := mainTest(t,
		"../../examples/iguana",
		"../../examples/chameleon",
		"../../examples/leopard",
	)
	if actual != expected {
		t.Errorf("Output was \n%q\n but expected \n%q", actual, expected)
	}
}

func mainTest(t *testing.T, args ...string) string {
	os.Args = append([]string{"tinystat"}, args...)

	f, err := ioutil.TempFile(os.TempDir(), "tinystat")
	if err != nil {
		t.Fatal(err)
	}
	defer func() {
		f.Close()
	}()

	oldStdout := os.Stdout
	defer func() {
		os.Stdout = oldStdout
	}()
	os.Stdout = f
	main()

	stdout, err := ioutil.ReadFile(f.Name())
	if err != nil {
		t.Fatal(err)
	}

	// strip everything of trailing whitespace
	lines := strings.Split(string(stdout), "\n")
	for i, line := range lines {
		lines[i] = strings.TrimRightFunc(line, unicode.IsSpace)
	}
	return strings.Join(lines, "\n")
}
