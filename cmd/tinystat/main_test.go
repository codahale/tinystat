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
	assert.Equal(t, "Output", want,
		mainTest(t, "../../examples/iguana"))
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

File       N  Mean    Stddev
iguana     7  300.00  238.05  (control)
chameleon  5  540.00  299.08  (no difference, p = 0.050)
`
	assert.Equal(t, "Output", want,
		mainTest(t,
			"../../examples/iguana",
			"../../examples/chameleon",
		))
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

File       N  Mean    Stddev
iguana     7  300.00  238.05  (control)
chameleon  5  540.00  299.08  (no difference, p = 0.050)
leopard    6  643.50  240.09  (643.50 > 300.00 ± 293.97, p = 0.050)
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
