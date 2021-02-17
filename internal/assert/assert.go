package assert

import (
	"testing"

	"github.com/google/go-cmp/cmp"
)

func Equal(tb testing.TB, name string, want, got interface{}, opts ...cmp.Option) {
	tb.Helper()

	if diff := cmp.Diff(want, got, opts...); diff != "" {
		tb.Errorf("%s mismatch (-want +got):\n%s", name, diff)
	}
}
