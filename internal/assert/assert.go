package assert

import (
	"testing"

	"github.com/google/go-cmp/cmp"
)

func Equal(t testing.TB, name string, want, got interface{}, opts ...cmp.Option) {
	t.Helper()

	if diff := cmp.Diff(want, got, opts...); diff != "" {
		t.Errorf("%s mismatch (-want +got):\n%s", name, diff)
	}
}
