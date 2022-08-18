// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//            National Center for Biotechnology Information (NCBI)
//
//  This software/database is a "United States Government Work" under the
//  terms of the United States Copyright Act. It was written as part of
//  the author's official duties as a United States Government employee and
//  thus cannot be copyrighted. This software/database is freely available
//  to the public for use. The National Library of Medicine and the U.S.
//  Government do not place any restriction on its use or reproduction.
//  We would, however, appreciate having the NCBI and the author cited in
//  any work or product based on this material.
//
//  Although all reasonable efforts have been taken to ensure the accuracy
//  and reliability of the software and data, the NLM and the U.S.
//  Government do not and cannot warrant the performance or results that
//  may be obtained by using this software or data. The NLM and the U.S.
//  Government disclaim all warranties, express or implied, including
//  warranties of performance, merchantability or fitness for any particular
//  purpose.
//
// ===========================================================================
//
// File Name:  split.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"strconv"
	"strings"
)

// PartitionPattern splits XML input from <pattern> to </pattern> and sends
// individual records to a callback. Requiring the input to be an XMLBlock
// channel of trimmed strings, generated by CreateXMLStreamer, simplifies the
// code by eliminating the need to check for an incomplete object tag at the end.
func PartitionPattern(pat, star string, turbo bool, inp <-chan XMLBlock, proc func(string)) {

	if pat == "" || inp == nil || proc == nil {
		return
	}

	// Scanner stores the precomputed Boyer-Moore-Horspool pattern matching table.
	// By experiment, this was slightly (but reproducibly) faster than the Boyer-Moore-Sunday variant.
	type Scanner struct {
		Pattern   string
		PatLength int
		CharSkip  [256]int
	}

	// newScanner initializes <pattern> to </pattern> scanner
	newScanner := func(pattern string) *Scanner {

		if pattern == "" {
			return nil
		}

		scr := &Scanner{Pattern: pattern}

		patlen := len(pattern)
		scr.PatLength = patlen

		// position of last character in pattern
		last := patlen - 1

		// initialize bad character displacement table
		for i := range scr.CharSkip {
			scr.CharSkip[i] = patlen
		}
		for i := 0; i < last; i++ {
			ch := pattern[i]
			scr.CharSkip[ch] = last - i
		}

		return scr
	}

	// isAnElement checks surroundings of match candidate
	isAnElement := func(text string, lf, rt, mx int) bool {

		if (lf >= 0 && text[lf] == '<') || (lf > 0 && text[lf] == '/' && text[lf-1] == '<') {
			if (rt < mx && (text[rt] == '>' || text[rt] == ' ' || text[rt] == '\n')) || (rt+1 < mx && text[rt] == '/' && text[rt+1] == '>') {
				return true
			}
		}

		return false
	}

	// findNextMatch is a modified Boyer-Moore-Horspool search function for maximum partitioning speed
	findNextMatch := func(scr *Scanner, text string, offset int) (int, int, int) {

		if scr == nil || text == "" {
			return -1, -1, -1
		}

		// copy values into local variables for speed
		txtlen := len(text)
		pattern := scr.Pattern[:]
		patlen := scr.PatLength
		max := txtlen - patlen
		last := patlen - 1
		skip := scr.CharSkip[:]

		i := offset

		for i <= max {
			j := last
			k := i + last
			for j >= 0 && text[k] == pattern[j] {
				j--
				k--
			}
			// require match candidate to be element name, i.e., <pattern ... >, </pattern ... >, or <pattern ... />
			if j < 0 && isAnElement(text, i-1, i+patlen, txtlen) {
				// find positions of flanking brackets
				lf := i - 1
				for lf > 0 && text[lf] != '<' {
					lf--
				}
				rt := i + patlen
				for rt < txtlen && text[rt] != '>' {
					rt++
				}
				return i + 1, lf, rt + 1
			}
			// find character in text above last character in pattern
			ch := text[i+last]
			// displacement table can shift pattern by one or more positions
			i += skip[ch]
		}

		return -1, -1, -1
	}

	// PatternType is the integer type for XML tag classification
	type PatternType int

	// PatternType keys for XML parsing
	const (
		NOPATTERN PatternType = iota
		STARTPATTERN
		SELFPATTERN
		STOPPATTERN
	)

	// nextPattern finds next element with pattern name
	nextPattern := func(scr *Scanner, text string, pos int) (PatternType, int, int, int) {

		if scr == nil || text == "" {
			return NOPATTERN, 0, 0, 0
		}

		prev := pos

		for {
			next, start, stop := findNextMatch(scr, text, prev)
			if next < 0 {
				return NOPATTERN, 0, 0, 0
			}

			prev = next + 1

			if text[start+1] == '/' {
				return STOPPATTERN, start, stop, prev
			} else if text[stop-2] == '/' {
				return SELFPATTERN, start, stop, prev
			} else {
				return STARTPATTERN, start, stop, prev
			}
		}
	}

	// doTurbo reads an XML file that has NEXT_RECORD_SIZE objects with
	// the number of bytes to read to get the next indexed pattern
	doTurbo := func() {

		line := ""
		var accumulator strings.Builder

		for {

			// read next XMLBlock ending with '>' character
			line = string(<-inp)
			if line == "" {
				return
			}

			// should be at next NEXT_RECORD_SIZE object
			for {

				// find start tag of next record size object
				idx := strings.Index(line, "<NEXT_RECORD_SIZE>")
				if idx < 0 {
					break
				}
				line = line[idx+18:]

				// if end of buffer, read next XMLBlock
				if line == "" {
					line = string(<-inp)
					if line == "" {
						return
					}
				}

				// find stop tag of next record size object
				idx = strings.Index(line, "</NEXT_RECORD_SIZE>")
				if idx < 0 {
					break
				}
				str := line[:idx]
				line = line[idx+19:]
				if strings.HasPrefix(line, "\n") {
					line = line[1:]
				}

				// convert object value to integer
				size, err := strconv.Atoi(str)
				if err != nil {
					break
				}

				accumulator.Reset()
				for {
					// size of remaining text in block
					max := len(line)

					// is record completely contained in current block
					if size < max {
						rec := line[:size]
						line = line[size:]
						// prepend any text collected from previous blocks
						prev := accumulator.String()
						res := prev + rec
						res = strings.TrimPrefix(res, "\n")
						res = strings.TrimSuffix(res, "\n")
						proc(res[:])
						break
					}

					// otherwise record remainder of block
					accumulator.WriteString(line)
					// decrement remaining size
					size -= len(line)
					// read next block
					line = string(<-inp)
					if line == "" {
						// last record on final block
						res := accumulator.String()
						res = strings.TrimPrefix(res, "\n")
						res = strings.TrimSuffix(res, "\n")
						proc(res[:])
						return
					}
					// and keep going until desired size is collected
				}
			}
		}
	}

	// doNormal handles -pattern Object construct, keeping track of nesting level
	doNormal := func() {

		// current depth of -pattern objects
		level := 0

		begin := 0
		inPattern := false

		line := ""
		var accumulator strings.Builder

		match := NOPATTERN
		start := 0
		stop := 0
		next := 0

		scr := newScanner(pat)
		if scr == nil {
			return
		}

		for {

			begin = 0
			next = 0

			line = string(<-inp)
			if line == "" {
				return
			}

			for {
				match, start, stop, next = nextPattern(scr, line, next)
				if match == STARTPATTERN {
					if level == 0 {
						inPattern = true
						begin = start
					}
					level++
				} else if match == STOPPATTERN {
					level--
					if level == 0 {
						inPattern = false
						accumulator.WriteString(line[begin:stop])
						// read and process one -pattern object at a time
						str := accumulator.String()
						if str != "" {
							proc(str[:])
						}
						// reset accumulator
						accumulator.Reset()
					}
				} else if match == SELFPATTERN {
					if level == 0 {
						str := line[start:stop]
						if str != "" {
							proc(str[:])
						}
					}
				} else {
					if inPattern {
						accumulator.WriteString(line[begin:])
					}
					break
				}
			}
		}
	}

	// doStar handles -pattern Parent/* construct for heterogeneous objects. It now works
	// with concatenated files, but not if components are recursive or self-closing objects.
	// Process the latter through transmute -format -self first.
	doStar := func() {

		// current depth of -pattern objects
		level := 0

		begin := 0
		inPattern := false

		line := ""
		var accumulator strings.Builder

		match := NOPATTERN
		start := 0
		stop := 0
		next := 0

		scr := newScanner(pat)
		if scr == nil {
			return
		}

		last := pat

		// read to first <pattern> element
		for {

			next = 0

			line = string(<-inp)
			if line == "" {
				break
			}

			match, start, stop, next = nextPattern(scr, line, next)
			if match == STARTPATTERN {
				break
			}
		}

		if match != STARTPATTERN {
			return
		}

		// find next element in XML
		nextElement := func(text string, pos int) string {

			txtlen := len(text)

			tag := ""
			for i := pos; i < txtlen; i++ {
				if text[i] == '<' {
					tag = text[i+1:]
					break
				}
			}
			if tag == "" {
				return ""
			}
			if tag[0] == '/' {
				if strings.HasPrefix(tag[1:], pat) {
					//should be </pattern> at end, want to continue if concatenated files
					return "/"
				}
				return ""
			}
			for i, ch := range tag {
				if ch == '>' || ch == ' ' || ch == '/' {
					return tag[0:i]
				}
			}

			return ""
		}

		// read and process heterogeneous objects immediately below <pattern> parent
		for {
			tag := nextElement(line, next)
			if tag == "" {

				begin = 0
				next = 0

				line = string(<-inp)
				if line == "" {
					break
				}

				tag = nextElement(line, next)
			}
			if tag == "" {
				return
			}

			// check for concatenated parent set files
			if tag[0] == '/' {
				scr = newScanner(pat)
				if scr == nil {
					return
				}
				last = pat
				// confirm end </pattern> just found
				match, start, stop, next = nextPattern(scr, line, next)
				if match != STOPPATTERN {
					return
				}
				// now look for a new start <pattern> tag
				for {
					match, start, stop, next = nextPattern(scr, line, next)
					if match == STARTPATTERN {
						break
					}
					next = 0
					line = string(<-inp)
					if line == "" {
						break
					}
				}
				if match != STARTPATTERN {
					return
				}
				// continue with processing loop
				continue
			}

			if tag != last {
				scr = newScanner(tag)
				if scr == nil {
					return
				}
				last = tag
			}

			for {
				match, start, stop, next = nextPattern(scr, line, next)
				if match == STARTPATTERN {
					if level == 0 {
						inPattern = true
						begin = start
					}
					level++
				} else if match == STOPPATTERN {
					level--
					if level == 0 {
						inPattern = false
						accumulator.WriteString(line[begin:stop])
						// read and process one -pattern/* object at a time
						str := accumulator.String()
						if str != "" {
							proc(str[:])
						}
						// reset accumulator
						accumulator.Reset()
						break
					}
				} else if match == SELFPATTERN {
					if level == 0 {
						str := line[start:stop]
						if str != "" {
							proc(str[:])
						}
					}
				} else {
					if inPattern {
						accumulator.WriteString(line[begin:])
					}

					begin = 0
					next = 0

					line = string(<-inp)
					if line == "" {
						break
					}
				}
			}
		}
	}

	// call appropriate handler
	if turbo {
		doTurbo()
	} else if star == "" {
		doNormal()
	} else if star == "*" {
		doStar()
	}
}
