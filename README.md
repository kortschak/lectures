# Bioinformatics course lecture slides 2014

These slide are essentially those used in my 2014 Bioinformatics course.

Changes have been made to avoid copyright issues; no images are included
except (mainly) for those I have made myself or that are freely available,
or with permission*. Others are linked directly from their source.

In some cases this means that the images will not be available unless viewed
from an institutional address - and even not then in some cases. *Wiley, I'm
looking at you.*

\* Exceptions to this are Swofford's, Altschul's, and Vingron-Stoye-and-Luz's
lecture notes.


# Using the slides

Slides are writting in present markup, so you will need to install the
[present tool](https://godoc.org/golang.org/x/tools/cmd/present). This
requires a [Go installation](https://golang.org/doc/install).

Then invoke present in the lectures directory with `present -base base`. If
runnable Go code examples are to be used you will need to install all the
dependencies listed in `imports.go` - `go get github.com/kortschak/lectures`
is the easiest way to do this.

Alternatively you can link to the [Introduction](http://talks.godoc.org/github.com/kortschak/lectures/01-bioinformatics.slide) from here (other slide decks depend on HTML which is not
supported by talks.godoc.org).
