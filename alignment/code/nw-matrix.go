// Copyright ©2014 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"errors"
	"fmt"
	"os"
	"text/tabwriter"

	"github.com/biogo/biogo/align/matrix"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
)

// Setting debugNeedle to true gives verbose scoring table output for the dynamic programming.
const debugNeedle = true

type AlphabetSlicer interface {
	Alphabet() alphabet.Alphabet
	Slice() alphabet.Slice
}

// An Aligner aligns the sequence data of two type-matching Slicers, returning an ordered
// slice of features describing matching and mismatching segments. The sequences to be aligned
// must have a valid gap letter in the first position of their alphabet; the alphabets
// {DNA,RNA}{gapped,redundant} and Protein provided by the alphabet package satisfy this.
type Aligner interface {
	Align(reference, query AlphabetSlicer) ([]feat.Pair, error)
}

// A Linear is a basic linear gap penalty alignment description.
// It is a square scoring matrix with the first column and first row specifying gap penalties.
type Linear [][]int

// An Affine is a basic affine gap penalty alignment description.
type Affine struct {
	Matrix  Linear
	GapOpen int
}

const (
	diag = iota
	up
	left

	gap = 0

	minInt = -int(^uint(0)>>1) - 1
)

var (
	ErrMismatchedTypes     = errors.New("align: mismatched sequence types")
	ErrMismatchedAlphabets = errors.New("align: mismatched alphabets")
	ErrNoAlphabet          = errors.New("align: no alphabet")
	ErrNotGappedAlphabet   = errors.New("align: alphabet does not have gap at position 0")
	ErrTypeNotHandled      = errors.New("align: sequence type not handled")
	ErrMatrixNotSquare     = errors.New("align: scoring matrix is not square")
)

func max(a *[3]int) int {
	m := minInt
	for _, v := range a {
		if v > m {
			m = v
		}
	}
	return m
}

func max2(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func add(a, b int) int {
	if a == minInt || b == minInt {
		return minInt
	}
	return a + b
}

type feature struct {
	start, end int
	loc        feat.Feature
}

func (f feature) Name() string {
	if f.loc != nil {
		return f.loc.Name()
	}
	return ""
}
func (f feature) Description() string {
	if f.loc != nil {
		return f.loc.Description()
	}
	return ""
}
func (f feature) Location() feat.Feature { return f.loc }
func (f feature) Start() int             { return f.start }
func (f feature) End() int               { return f.end }
func (f feature) Len() int               { return f.end - f.start }

type featPair struct {
	a, b  feature
	score int
}

func (fp *featPair) Features() [2]feat.Feature { return [2]feat.Feature{fp.a, fp.b} }
func (fp *featPair) Score() int                { return fp.score }
func (fp *featPair) Invert()                   { fp.a, fp.b = fp.b, fp.a }
func (fp *featPair) String() string {
	switch {
	case fp.a.start == fp.a.end:
		return fmt.Sprintf("-/%s[%d,%d)=%d",
			fp.b.Name(), fp.b.start, fp.b.end,
			fp.score)
	case fp.b.start == fp.b.end:
		return fmt.Sprintf("%s[%d,%d)/-=%d",
			fp.a.Name(), fp.a.start, fp.a.end,
			fp.score)
	}
	return fmt.Sprintf("%s[%d,%d)/%s[%d,%d)=%d",
		fp.a.Name(), fp.a.start, fp.a.end,
		fp.b.Name(), fp.b.start, fp.b.end,
		fp.score)
}

// Format returns a [2]alphabet.Slice representing the formatted alignment of a and b described by the
// list of feature pairs in f, with gap used to fill gaps in the alignment.
func Format(a, b seq.Slicer, f []feat.Pair, gap alphabet.Letter) [2]alphabet.Slice {
	var as, aln [2]alphabet.Slice
	for i, s := range [2]seq.Slicer{a, b} {
		as[i] = s.Slice()
		aln[i] = as[i].Make(0, 0)
	}
	for _, fs := range f {
		fc := fs.Features()
		for i := range aln {
			if fc[i].Len() == 0 {
				switch aln[i].(type) {
				case alphabet.Letters:
					aln[i] = aln[i].Append(alphabet.Letters(gap.Repeat(fc[1-i].Len())))
				case alphabet.QLetters:
					aln[i] = aln[i].Append(alphabet.QLetters(alphabet.QLetter{L: gap}.Repeat(fc[1-i].Len())))
				}
			} else {
				aln[i] = aln[i].Append(as[i].Slice(fc[i].Start(), fc[i].End()))
			}
		}
	}
	return aln
}

// NW is the linear gap penalty Needleman-Wunsch aligner type.
type NW Linear

// Align aligns two sequences using the Needleman-Wunsch algorithm. It returns an alignment description
// or an error if the scoring matrix is not square, or the sequence data types or alphabets do not match.
func (a NW) Align(reference, query AlphabetSlicer) ([]feat.Pair, error) {
	alpha := reference.Alphabet()
	if alpha == nil {
		return nil, ErrNoAlphabet
	}
	if alpha != query.Alphabet() {
		return nil, ErrMismatchedAlphabets
	}
	if alpha.IndexOf(alpha.Gap()) != 0 {
		return nil, ErrNotGappedAlphabet
	}
	switch rSeq := reference.Slice().(type) {
	case alphabet.Letters:
		qSeq, ok := query.Slice().(alphabet.Letters)
		if !ok {
			return nil, ErrMismatchedTypes
		}
		return a.alignLetters(rSeq, qSeq, alpha)
	default:
		return nil, ErrTypeNotHandled
	}
}

func drawNWTableLetters(rSeq, qSeq alphabet.Letters, index alphabet.Index, table []int, a [][]int) {
	tw := tabwriter.NewWriter(os.Stdout, 0, 0, 0, ' ', tabwriter.AlignRight|tabwriter.Debug)
	fmt.Printf("rSeq: %s\n", rSeq)
	fmt.Printf("qSeq: %s\n", qSeq)
	fmt.Fprint(tw, "\tqSeq\t")
	for _, l := range qSeq {
		fmt.Fprintf(tw, "%c\t", l)
	}
	fmt.Fprintln(tw)

	r, c := rSeq.Len()+1, qSeq.Len()+1
	fmt.Fprint(tw, "rSeq\t")
	for i := 0; i < r; i++ {
		if i != 0 {
			fmt.Fprintf(tw, "%c\t", rSeq[i-1])
		}

		for j := 0; j < c; j++ {
			p := pointerNWLetters(rSeq, qSeq, i, j, table, index, a, c)
			if p != "" {
				fmt.Fprintf(tw, "%s % 3v\t", p, table[i*c+j])
			} else {
				fmt.Fprintf(tw, "%v\t", table[i*c+j])
			}
		}
		fmt.Fprintln(tw)
	}
	tw.Flush()
}

func pointerNWLetters(rSeq, qSeq alphabet.Letters, i, j int, table []int, index alphabet.Index, a [][]int, c int) string {
	switch {
	case i == 0 && j == 0:
		return " "
	case i == 0:
		return "⬅"
	case j == 0:
		return "⬆"
	}
	rVal := index[rSeq[i-1]]
	qVal := index[qSeq[j-1]]
	if rVal < 0 || qVal < 0 {
		return " "
	}
	switch p := i*c + j; table[p] {
	case table[p-c-1] + a[rVal][qVal]:
		return "⬉"
	case table[p-c] + a[rVal][gap]:
		return "⬆"
	case table[p-1] + a[gap][qVal]:
		return "⬅"
	default:
		return " "
	}
}

func (a NW) alignLetters(rSeq, qSeq alphabet.Letters, alpha alphabet.Alphabet) ([]feat.Pair, error) {
	let := len(a)
	la := make([]int, 0, let*let)
	for _, row := range a {
		if len(row) != let {
			return nil, ErrMatrixNotSquare
		}
		la = append(la, row...)
	}

	index := alpha.LetterIndex()
	r, c := rSeq.Len()+1, qSeq.Len()+1
	table := make([]int, r*c)
	for j := range table[1:c] {
		table[j+1] = table[j] + la[index[qSeq[j]]]
	}
	for i := 1; i < r; i++ {
		table[i*c] = table[(i-1)*c] + la[index[rSeq[i-1]]*let]
	}

	var scores [3]int
	for i := 1; i < r; i++ {
		for j := 1; j < c; j++ {
			var (
				rVal = index[rSeq[i-1]]
				qVal = index[qSeq[j-1]]
			)
			if rVal < 0 || qVal < 0 {
				continue
			}
			p := i*c + j
			scores = [3]int{
				diag: table[p-c-1] + la[rVal*let+qVal],
				up:   table[p-c] + la[rVal*let],
				left: table[p-1] + la[qVal],
			}
			table[p] = max(&scores)
		}
	}
	if debugNeedle {
		drawNWTableLetters(rSeq, qSeq, index, table, a)
	}

	var aln []feat.Pair
	score, last := 0, diag
	i, j := r-1, c-1
	maxI, maxJ := i, j
	for i > 0 && j > 0 {
		var (
			rVal = index[rSeq[i-1]]
			qVal = index[qSeq[j-1]]
		)
		if rVal < 0 || qVal < 0 {
			continue
		}
		switch p := i*c + j; table[p] {
		case table[p-c-1] + la[rVal*let+qVal]:
			if last != diag {
				aln = append(aln, &featPair{
					a:     feature{start: i, end: maxI},
					b:     feature{start: j, end: maxJ},
					score: score,
				})
				maxI, maxJ = i, j
				score = 0
			}
			score += table[p] - table[p-c-1]
			i--
			j--
			last = diag
		case table[p-c] + la[rVal*let]:
			if last != up && p != len(table)-1 {
				aln = append(aln, &featPair{
					a:     feature{start: i, end: maxI},
					b:     feature{start: j, end: maxJ},
					score: score,
				})
				maxI, maxJ = i, j
				score = 0
			}
			score += table[p] - table[p-c]
			i--
			last = up
		case table[p-1] + la[qVal]:
			if last != left && p != len(table)-1 {
				aln = append(aln, &featPair{
					a:     feature{start: i, end: maxI},
					b:     feature{start: j, end: maxJ},
					score: score,
				})
				maxI, maxJ = i, j
				score = 0
			}
			score += table[p] - table[p-1]
			j--
			last = left
		default:
			panic(fmt.Sprintf("align: nw internal error: no path at row: %d col:%d\n", i, j))
		}
	}

	aln = append(aln, &featPair{
		a:     feature{start: i, end: maxI},
		b:     feature{start: j, end: maxJ},
		score: score,
	})
	if i != j {
		aln = append(aln, &featPair{
			a:     feature{start: 0, end: i},
			b:     feature{start: 0, end: j},
			score: table[i*c+j],
		})
	}

	for i, j := 0, len(aln)-1; i < j; i, j = i+1, j-1 {
		aln[i], aln[j] = aln[j], aln[i]
	}

	return aln, nil
}

type Scorer interface {
	Score() int
}

func init() {
	for i := 1; i < len(matrix.IDENTITY); i++ {
		matrix.IDENTITY[i] = make([]int, len(matrix.IDENTITY))
		matrix.IDENTITY[i][i] = 1
	}
}

func setGap(m [][]int, gap int) NW {
	for i := 1; i < len(m); i++ {
		m[0][i] = gap
		m[i][0] = gap
	}
	return m
}

func main() {
	nwsa := &linear.Seq{Seq: alphabet.BytesToLetters([]byte("VLVVFVL"))}
	nwsa.Alpha = alphabet.Protein
	nwsb := &linear.Seq{Seq: alphabet.BytesToLetters([]byte("ILIFIPMLQ"))}
	nwsb.Alpha = alphabet.Protein

	needle := setGap(matrix.IDENTITY, -10) // (Match: 1, Mismatch: 0, Gap: -10)

	aln, err := needle.Align(nwsa, nwsb)
	if err == nil {
		var score int
		for _, fp := range aln {
			score += fp.(Scorer).Score()
		}
		fmt.Printf("\n%s Total: %d\n\n", aln, score)
		fa := Format(nwsa, nwsb, aln, '-')
		fmt.Printf("%s\n%s\n", fa[0], fa[1])
	}
}
