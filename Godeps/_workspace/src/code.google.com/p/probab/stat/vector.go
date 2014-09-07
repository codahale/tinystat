package stat

type Vector struct {
	X []float64 // data
	L int       // length
}

func NewVector(length int) (v *Vector) {
	v = new(Vector)
	v.L = length
	v.X = make([]float64, length)
	return v
}

func (v Vector) Set(i int, x float64) {
	v.X[i] = x
}

func (v Vector) Get(i int) float64 {
	return v.X[i]
}

func (v Vector) Swap(i int, j int) {
	x := v.X[i]
	v.X[i] = v.X[j]
	v.X[j] = x
}

func (v Vector) Len() int {
	return v.L
}
