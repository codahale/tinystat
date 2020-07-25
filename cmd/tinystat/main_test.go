package main

import (
	"io/ioutil"
	"os"
	"strings"
	"testing"
	"unicode"

	"github.com/codahale/tinystat/internal/assert"
)

func TestControlOnly(t *testing.T) {
	want := `
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
	assert.Equal(t, "Output", want, mainTest(t, "../../examples/iguana"))
}

func TestOneExperiment(t *testing.T) {
	want := `
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
	assert.Equal(t, "Output", want,
		mainTest(t, "../../examples/iguana", "../../examples/chameleon"))
}

func TestTwoExperiments(t *testing.T) {
	want := `
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
                            343.5 +/- 292.6345386382977
                            114.5% +/- 97.54484621276592%
                            (Student's t, pooled s = 238.9799344943192)
`
	assert.Equal(t, "Output", want,
		mainTest(t,
			"../../examples/iguana",
			"../../examples/chameleon",
			"../../examples/leopard",
		))
}

func mainTest(t *testing.T, args ...string) string {
	os.Args = append([]string{"tinystat"}, args...)

	f, err := ioutil.TempFile(os.TempDir(), "tinystat")
	if err != nil {
		t.Fatal(err)
	}

	defer func() {
		_ = f.Close()
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
