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
// File Name:  gbff.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"bufio"
	"bytes"
	"fmt"
	"github.com/klauspost/pgzip"
	"github.com/surgebase/porter2"
	"html"
	"io"
	"os"
	"path"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"sync"
	"unicode"
)

// GenBankConverter reads flatfiles and sends INSDSeq XML records down a channel
func GenBankConverter(inp io.Reader) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "Unable to create GenBank converter channel\n")
		os.Exit(1)
	}

	const twelvespaces = "            "
	const twentyonespaces = "                     "

	var rec strings.Builder
	var alt strings.Builder
	var con strings.Builder
	var seq strings.Builder

	convertGenBank := func(inp io.Reader, out chan<- string) {

		// close channel when all records have been sent
		defer close(out)

		scanr := bufio.NewScanner(inp)

		row := 0

		nextLine := func() string {

			for scanr.Scan() {
				line := scanr.Text()
				if line == "" {
					continue
				}
				return line
			}
			return ""

		}

		for {

			rec.Reset()

			// read first line of next record
			line := nextLine()
			if line == "" {
				break
			}

			row++

			for {
				if line == "" {
					break
				}
				if !strings.HasPrefix(line, "LOCUS") {
					// skip release file header information
					line = nextLine()
					row++
					continue
				}
				break
			}

			readContinuationLines := func(str string) string {

				for {
					// read next line
					line = nextLine()
					row++
					if !strings.HasPrefix(line, twelvespaces) {
						// if not continuation line, break out of loop
						break
					}
					// append subsequent line and continue with loop
					txt := strings.TrimPrefix(line, twelvespaces)
					str += " " + txt
				}

				str = CompressRunsOfSpaces(str)
				str = strings.TrimSpace(str)

				return str
			}

			writeOneElement := func(spaces, tag, value string) {

				rec.WriteString(spaces)
				rec.WriteString("<")
				rec.WriteString(tag)
				rec.WriteString(">")
				value = html.EscapeString(value)
				rec.WriteString(value)
				rec.WriteString("</")
				rec.WriteString(tag)
				rec.WriteString(">\n")
			}

			// each section will exit with the next line ready to process

			if strings.HasPrefix(line, "LOCUS") {

				// start of record
				rec.WriteString("  <INSDSeq>\n")

				// do not break if given artificial multi-line LOCUS
				str := readContinuationLines(line)

				cols := strings.Fields(str)
				ln := len(cols)
				if ln == 8 {
					moleculetype := cols[4]
					strandedness := ""
					if strings.HasPrefix(moleculetype, "ds-") {
						moleculetype = strings.TrimPrefix(moleculetype, "ds-")
						strandedness = "double"
					} else if strings.HasPrefix(moleculetype, "ss-") {
						moleculetype = strings.TrimPrefix(moleculetype, "ss-")
						strandedness = "single"
					} else if strings.HasPrefix(moleculetype, "ms-") {
						moleculetype = strings.TrimPrefix(moleculetype, "ms-")
						strandedness = "mixed"
					} else if strings.HasSuffix(moleculetype, "DNA") {
						strandedness = "double"
					} else if strings.HasSuffix(moleculetype, "RNA") {
						strandedness = "single"
					}

					writeOneElement("    ", "INSDSeq_locus", cols[1])
					writeOneElement("    ", "INSDSeq_length", cols[2])

					if strandedness != "" {
						writeOneElement("    ", "INSDSeq_strandedness", strandedness)
					}

					writeOneElement("    ", "INSDSeq_moltype", moleculetype)
					writeOneElement("    ", "INSDSeq_topology", cols[5])
					writeOneElement("    ", "INSDSeq_division", cols[6])
					writeOneElement("    ", "INSDSeq_update-date", cols[7])

				} else if ln == 7 {

					writeOneElement("    ", "INSDSeq_locus", cols[1])
					writeOneElement("    ", "INSDSeq_length", cols[2])
					writeOneElement("    ", "INSDSeq_moltype", "AA")
					writeOneElement("    ", "INSDSeq_topology", cols[4])
					writeOneElement("    ", "INSDSeq_division", cols[5])
					writeOneElement("    ", "INSDSeq_update-date", cols[6])

				} else {
					fmt.Fprintf(os.Stderr, "ERROR: "+str+"\n")
				}

				// read next line and continue - handled by readContinuationLines above
				// line = nextLine()
				// row++
			}

			if strings.HasPrefix(line, "DEFINITION") {

				txt := strings.TrimPrefix(line, "DEFINITION")
				def := readContinuationLines(txt)
				def = strings.TrimSuffix(def, ".")

				writeOneElement("    ", "INSDSeq_definition", def)
			}

			var secondaries []string

			if strings.HasPrefix(line, "ACCESSION") {

				txt := strings.TrimPrefix(line, "ACCESSION")
				str := readContinuationLines(txt)
				accessions := strings.Fields(str)
				ln := len(accessions)
				if ln > 1 {

					writeOneElement("    ", "INSDSeq_primary-accession", accessions[0])

					// skip past primary accession, collect secondaries
					secondaries = accessions[1:]

				} else if ln == 1 {

					writeOneElement("    ", "INSDSeq_primary-accession", accessions[0])

				} else {
					fmt.Fprintf(os.Stderr, "ERROR: ACCESSION "+str+"\n")
				}
			}

			accnver := ""
			gi := ""

			if strings.HasPrefix(line, "VERSION") {

				cols := strings.Fields(line)
				if len(cols) == 2 {

					accnver = cols[1]
					writeOneElement("    ", "INSDSeq_accession-version", accnver)

				} else if len(cols) == 3 {

					accnver = cols[1]
					writeOneElement("    ", "INSDSeq_accession-version", accnver)

					// collect gi for other-seqids
					if strings.HasPrefix(cols[2], "GI:") {
						gi = strings.TrimPrefix(cols[2], "GI:")
					}

				} else {
					fmt.Fprintf(os.Stderr, "ERROR: "+line+"\n")
				}

				// read next line and continue
				line = nextLine()
				row++
			}

			if gi != "" {

				rec.WriteString("    <INSDSeq_other-seqids>\n")
				writeOneElement("      ", "INSDSeqid", "gi|"+gi)
				rec.WriteString("    </INSDSeq_other-seqids>\n")
			}

			if len(secondaries) > 0 {

				rec.WriteString("    <INSDSeq_secondary-accessions>\n")

				for _, secndry := range secondaries {

					if strings.HasPrefix(secndry, "REGION") {
						break
					}
					writeOneElement("      ", "INSDSecondary-accn", secndry)
				}

				rec.WriteString("    </INSDSeq_secondary-accessions>\n")

				// reset secondaries slice
				secondaries = nil
			}

			dblink := ""

			if strings.HasPrefix(line, "DBLINK") {

				txt := strings.TrimPrefix(line, "DBLINK")
				dblink = readContinuationLines(txt)
			}

			srcdb := ""

			if strings.HasPrefix(line, "DBSOURCE") {

				txt := strings.TrimPrefix(line, "DBSOURCE")
				srcdb = readContinuationLines(txt)
			}

			if strings.HasPrefix(line, "KEYWORDS") {

				txt := strings.TrimPrefix(line, "KEYWORDS")
				key := readContinuationLines(txt)
				key = strings.TrimSuffix(key, ".")

				if key != "" {
					rec.WriteString("    <INSDSeq_keywords>\n")
					kywds := strings.Split(key, ";")
					for _, kw := range kywds {
						kw = strings.TrimSpace(kw)
						if kw == "" || kw == "." {
							continue
						}

						writeOneElement("      ", "INSDKeyword", kw)
					}
					rec.WriteString("    </INSDSeq_keywords>\n")
				}
			}

			if strings.HasPrefix(line, "SOURCE") {

				txt := strings.TrimPrefix(line, "SOURCE")
				src := readContinuationLines(txt)

				writeOneElement("    ", "INSDSeq_source", src)
			}

			if strings.HasPrefix(line, "  ORGANISM") {

				org := strings.TrimPrefix(line, "  ORGANISM")
				org = CompressRunsOfSpaces(org)
				org = strings.TrimSpace(org)

				writeOneElement("    ", "INSDSeq_organism", org)

				line = nextLine()
				row++
				if strings.HasPrefix(line, twelvespaces) {
					txt := strings.TrimPrefix(line, twelvespaces)
					tax := readContinuationLines(txt)
					tax = strings.TrimSuffix(tax, ".")

					writeOneElement("    ", "INSDSeq_taxonomy", tax)
				}
			}

			rec.WriteString("    <INSDSeq_references>\n")
			for {
				if !strings.HasPrefix(line, "REFERENCE") {
					// exit out of reference section
					break
				}

				ref := "0"

				rec.WriteString("      <INSDReference>\n")

				txt := strings.TrimPrefix(line, "REFERENCE")
				str := readContinuationLines(txt)
				str = CompressRunsOfSpaces(str)
				str = strings.TrimSpace(str)
				idx := strings.Index(str, "(")
				if idx > 0 {
					ref = strings.TrimSpace(str[:idx])

					writeOneElement("        ", "INSDReference_reference", ref)

					posn := str[idx+1:]
					posn = strings.TrimSuffix(posn, ")")
					posn = strings.TrimSpace(posn)
					if posn == "sites" {

						writeOneElement("        ", "INSDReference_position", posn)

					} else {
						var arry []string
						cls := strings.Split(posn, ";")
						for _, item := range cls {
							item = strings.TrimPrefix(item, "bases ")
							item = strings.TrimPrefix(item, "residues ")
							item = strings.TrimSpace(item)
							cols := strings.Fields(item)
							if len(cols) == 3 && cols[1] == "to" {
								arry = append(arry, cols[0]+".."+cols[2])
							}
						}
						if len(arry) > 0 {
							posit := strings.Join(arry, ",")
							writeOneElement("        ", "INSDReference_position", posit)
						} else {
							fmt.Fprintf(os.Stderr, "ERROR: "+posn+"\n")
						}
					}
				} else {
					ref = strings.TrimSpace(str)

					writeOneElement("        ", "INSDReference_reference", ref)
				}
				row++

				if strings.HasPrefix(line, "  AUTHORS") {

					txt := strings.TrimPrefix(line, "  AUTHORS")
					auths := readContinuationLines(txt)

					rec.WriteString("        <INSDReference_authors>\n")
					authors := strings.Split(auths, ", ")
					for _, auth := range authors {
						auth = strings.TrimSpace(auth)
						if auth == "" {
							continue
						}
						pair := strings.Split(auth, " and ")
						for _, name := range pair {

							writeOneElement("          ", "INSDAuthor", name)
						}
					}
					rec.WriteString("        </INSDReference_authors>\n")
				}

				if strings.HasPrefix(line, "  CONSRTM") {

					txt := strings.TrimPrefix(line, "  CONSRTM")
					cons := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_consortium", cons)
				}

				if strings.HasPrefix(line, "  TITLE") {

					txt := strings.TrimPrefix(line, "  TITLE")
					titl := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_title", titl)
				}

				if strings.HasPrefix(line, "  JOURNAL") {

					txt := strings.TrimPrefix(line, "  JOURNAL")
					jour := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_journal", jour)
				}

				if strings.HasPrefix(line, "   PUBMED") {

					txt := strings.TrimPrefix(line, "   PUBMED")
					pmid := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_pubmed", pmid)
				}

				if strings.HasPrefix(line, "  MEDLINE") {

					txt := strings.TrimPrefix(line, "  MEDLINE")
					// old MEDLINE uid not supported
					readContinuationLines(txt)
				}

				if strings.HasPrefix(line, "  REMARK") {

					txt := strings.TrimPrefix(line, "  REMARK")
					rem := readContinuationLines(txt)

					writeOneElement("        ", "INSDReference_remark", rem)
				}

				// end of this reference
				rec.WriteString("      </INSDReference>\n")
				// continue to next reference
			}
			rec.WriteString("    </INSDSeq_references>\n")

			if strings.HasPrefix(line, "COMMENT") {

				txt := strings.TrimPrefix(line, "COMMENT")
				com := readContinuationLines(txt)

				writeOneElement("    ", "INSDSeq_comment", com)
			}

			if strings.HasPrefix(line, "PRIMARY") {

				txt := strings.TrimPrefix(line, "PRIMARY")
				pmy := readContinuationLines(txt)

				writeOneElement("    ", "INSDSeq_primary", pmy)
			}

			if srcdb != "" {
				writeOneElement("    ", "INSDSeq_source-db", srcdb)
			}

			rec.WriteString("    <INSDSeq_feature-table>\n")
			if strings.HasPrefix(line, "FEATURES") {

				line = nextLine()
				row++

				for {
					if !strings.HasPrefix(line, "     ") {
						// exit out of features section
						break
					}
					if len(line) < 22 {
						fmt.Fprintf(os.Stderr, "ERROR: "+line+"\n")
						line = nextLine()
						row++
						continue
					}

					rec.WriteString("      <INSDFeature>\n")

					// read feature key and start of location
					fkey := line[5:21]
					fkey = strings.TrimSpace(fkey)

					writeOneElement("        ", "INSDFeature_key", fkey)

					loc := line[21:]
					loc = strings.TrimSpace(loc)
					for {
						line = nextLine()
						row++
						if !strings.HasPrefix(line, twentyonespaces) {
							break
						}
						txt := strings.TrimPrefix(line, twentyonespaces)
						if strings.HasPrefix(txt, "/") {
							// if not continuation of location, break out of loop
							break
						}
						// append subsequent line and continue with loop
						loc += strings.TrimSpace(txt)
					}

					writeOneElement("        ", "INSDFeature_location", loc)

					locationOperator := ""
					isComp := false
					prime5 := false
					prime3 := false

					// parseloc recursive definition
					var parseloc func(string) []string

					parseloc = func(str string) []string {

						var acc []string

						if strings.HasPrefix(str, "join(") && strings.HasSuffix(str, ")") {

							locationOperator = "join"

							str = strings.TrimPrefix(str, "join(")
							str = strings.TrimSuffix(str, ")")
							items := strings.Split(str, ",")

							for _, thisloc := range items {
								inner := parseloc(thisloc)
								for _, sub := range inner {
									acc = append(acc, sub)
								}
							}

						} else if strings.HasPrefix(str, "order(") && strings.HasSuffix(str, ")") {

							locationOperator = "order"

							str = strings.TrimPrefix(str, "order(")
							str = strings.TrimSuffix(str, ")")
							items := strings.Split(str, ",")

							for _, thisloc := range items {
								inner := parseloc(thisloc)
								for _, sub := range inner {
									acc = append(acc, sub)
								}
							}

						} else if strings.HasPrefix(str, "complement(") && strings.HasSuffix(str, ")") {

							isComp = true

							str = strings.TrimPrefix(str, "complement(")
							str = strings.TrimSuffix(str, ")")
							items := parseloc(str)

							// reverse items
							for i, j := 0, len(items)-1; i < j; i, j = i+1, j-1 {
								items[i], items[j] = items[j], items[i]
							}

							// reverse from and to positions, flip direction of angle brackets (partial flags)
							for _, thisloc := range items {
								pts := strings.Split(thisloc, "..")
								ln := len(pts)
								if ln == 2 {
									fst := pts[0]
									scd := pts[1]
									lf := ""
									rt := ""
									if strings.HasPrefix(fst, "<") {
										fst = strings.TrimPrefix(fst, "<")
										rt = ">"
									}
									if strings.HasPrefix(scd, ">") {
										scd = strings.TrimPrefix(scd, ">")
										lf = "<"
									}
									acc = append(acc, lf+scd+".."+rt+fst)
								} else if ln > 0 {
									acc = append(acc, pts[0])
								}
							}

						} else {

							// save individual interval or point if no leading accession
							if strings.Index(str, ":") < 0 {
								acc = append(acc, str)
							}
						}

						return acc
					}

					items := parseloc(loc)

					rec.WriteString("        <INSDFeature_intervals>\n")

					numIvals := 0

					// report individual intervals
					for _, thisloc := range items {
						if thisloc == "" {
							continue
						}

						numIvals++

						rec.WriteString("          <INSDInterval>\n")
						pts := strings.Split(thisloc, "..")
						if len(pts) == 2 {

							// fr..to
							fr := pts[0]
							to := pts[1]
							if strings.HasPrefix(fr, "<") {
								fr = strings.TrimPrefix(fr, "<")
								prime5 = true
							}
							if strings.HasPrefix(to, ">") {
								to = strings.TrimPrefix(to, ">")
								prime3 = true
							}
							writeOneElement("            ", "INSDInterval_from", fr)
							writeOneElement("            ", "INSDInterval_to", to)
							if isComp {
								rec.WriteString("            <INSDInterval_iscomp value=\"true\"/>\n")
							}
							writeOneElement("            ", "INSDInterval_accession", accnver)

						} else {

							crt := strings.Split(thisloc, "^")
							if len(crt) == 2 {

								// fr^to
								fr := crt[0]
								to := crt[1]
								writeOneElement("            ", "INSDInterval_from", fr)
								writeOneElement("            ", "INSDInterval_to", to)
								if isComp {
									rec.WriteString("            <INSDInterval_iscomp value=\"true\"/>\n")
								}
								rec.WriteString("            <INSDInterval_interbp value=\"true\"/>\n")
								writeOneElement("            ", "INSDInterval_accession", accnver)

							} else {

								// pt
								pt := pts[0]
								if strings.HasPrefix(pt, "<") {
									pt = strings.TrimPrefix(pt, "<")
									prime5 = true
								}
								if strings.HasPrefix(pt, ">") {
									pt = strings.TrimPrefix(pt, ">")
									prime3 = true
								}
								writeOneElement("            ", "INSDInterval_point", pt)
								writeOneElement("            ", "INSDInterval_accession", accnver)
							}
						}
						rec.WriteString("          </INSDInterval>\n")
					}

					rec.WriteString("        </INSDFeature_intervals>\n")

					if numIvals > 1 {
						writeOneElement("        ", "INSDFeature_operator", locationOperator)
					}
					if prime5 {
						rec.WriteString("        <INSDFeature_partial5 value=\"true\"/>\n")
					}
					if prime3 {
						rec.WriteString("        <INSDFeature_partial3 value=\"true\"/>\n")
					}

					hasQual := false
					for {
						if !strings.HasPrefix(line, twentyonespaces) {
							// if not qualifier line, break out of loop
							break
						}
						txt := strings.TrimPrefix(line, twentyonespaces)
						qual := ""
						val := ""
						if strings.HasPrefix(txt, "/") {
							if !hasQual {
								hasQual = true
								rec.WriteString("        <INSDFeature_quals>\n")
							}
							// read new qualifier and start of value
							qual = strings.TrimPrefix(txt, "/")
							qual = strings.TrimSpace(qual)
							idx := strings.Index(qual, "=")
							if idx > 0 {
								val = qual[idx+1:]
								qual = qual[:idx]
							}

							for {
								line = nextLine()
								row++
								if !strings.HasPrefix(line, twentyonespaces) {
									break
								}
								txt := strings.TrimPrefix(line, twentyonespaces)
								if strings.HasPrefix(txt, "/") {
									// if not continuation of qualifier, break out of loop
									break
								}
								// append subsequent line to value and continue with loop
								if qual == "transcription" || qual == "translation" || qual == "peptide" || qual == "anticodon" {
									val += strings.TrimSpace(txt)
								} else {
									val += " " + strings.TrimSpace(txt)
								}
							}

							rec.WriteString("          <INSDQualifier>\n")

							writeOneElement("            ", "INSDQualifier_name", qual)

							val = strings.TrimPrefix(val, "\"")
							val = strings.TrimSuffix(val, "\"")
							val = strings.TrimSpace(val)
							if val != "" {

								writeOneElement("            ", "INSDQualifier_value", val)
							}

							rec.WriteString("          </INSDQualifier>\n")
						}
					}
					if hasQual {
						rec.WriteString("        </INSDFeature_quals>\n")
					}

					// end of this feature
					rec.WriteString("      </INSDFeature>\n")
					// continue to next feature
				}
			}
			rec.WriteString("    </INSDSeq_feature-table>\n")

			// TSA, TLS, WGS, or CONTIG lines may be next

			altName := ""

			if strings.HasPrefix(line, "TSA") ||
				strings.HasPrefix(line, "TLS") ||
				strings.HasPrefix(line, "WGS") {

				alt.Reset()

				altName = line[:3]
				line = line[3:]
			}

			if strings.HasPrefix(line, "WGS_CONTIG") ||
				strings.HasPrefix(line, "WGS_SCAFLD") {

				alt.Reset()

				altName = line[:3]
				line = line[10:]
			}

			if altName != "" {

				altName = strings.ToLower(altName)
				txt := strings.TrimSpace(line)
				alt.WriteString(txt)
				for {
					// read next line
					line = nextLine()
					row++
					if !strings.HasPrefix(line, twelvespaces) {
						// if not continuation of contig, break out of loop
						break
					}
					// append subsequent line and continue with loop
					txt = strings.TrimPrefix(line, twelvespaces)
					txt = strings.TrimSpace(txt)
					alt.WriteString(txt)
				}
			}

			if strings.HasPrefix(line, "CONTIG") {

				// pathological records can have over 90,000 components, use strings.Builder
				con.Reset()

				txt := strings.TrimPrefix(line, "CONTIG")
				txt = strings.TrimSpace(txt)
				con.WriteString(txt)
				for {
					// read next line
					line = nextLine()
					row++
					if !strings.HasPrefix(line, twelvespaces) {
						// if not continuation of contig, break out of loop
						break
					}
					// append subsequent line and continue with loop
					txt = strings.TrimPrefix(line, twelvespaces)
					txt = strings.TrimSpace(txt)
					con.WriteString(txt)
				}
			}

			if strings.HasPrefix(line, "BASE COUNT") {

				txt := strings.TrimPrefix(line, "BASE COUNT")
				readContinuationLines(txt)
				// not supported
			}

			if strings.HasPrefix(line, "ORIGIN") {

				line = nextLine()
				row++
			}

			// remainder should be sequence

			// sequence can be millions of bases, use strings.Builder
			seq.Reset()

			for line != "" {

				if strings.HasPrefix(line, "//") {

					// end of record, print collected sequence
					str := seq.String()
					if str != "" {

						writeOneElement("    ", "INSDSeq_sequence", str)
					}
					seq.Reset()

					// print contig section
					str = con.String()
					str = strings.TrimSpace(str)
					if str != "" {
						writeOneElement("    ", "INSDSeq_contig", str)
					}
					con.Reset()

					if altName != "" {
						rec.WriteString("    <INSDSeq_alt-seq>\n")
						rec.WriteString("      <INSDAltSeqData>\n")
						str = alt.String()
						str = strings.TrimSpace(str)
						if str != "" {
							writeOneElement("        ", "INSDAltSeqData_name", altName)
							rec.WriteString("        <INSDAltSeqData_items>\n")
							fr, to := SplitInTwoLeft(str, "-")
							if fr != "" && to != "" {
								rec.WriteString("          <INSDAltSeqItem>\n")
								writeOneElement("            ", "INSDAltSeqItem_first-accn", fr)
								writeOneElement("            ", "INSDAltSeqItem_last-accn", to)
								rec.WriteString("          </INSDAltSeqItem>\n")
							} else {
								writeOneElement("          ", "INSDAltSeqItem_value", str)
							}
							rec.WriteString("        </INSDAltSeqData_items>\n")
						}
						alt.Reset()
						rec.WriteString("      </INSDAltSeqData>\n")
						rec.WriteString("    </INSDSeq_alt-seq>\n")
					}

					if dblink != "" {
						rec.WriteString("    <INSDSeq_xrefs>\n")
						// collect for database-reference
						flds := strings.Fields(dblink)
						for len(flds) > 1 {
							tag := flds[0]
							val := flds[1]
							flds = flds[2:]
							tag = strings.TrimSuffix(tag, ":")
							rec.WriteString("      <INSDXref>\n")
							writeOneElement("        ", "INSDXref_dbname", tag)
							writeOneElement("        ", "INSDXref_id", val)
							rec.WriteString("      </INSDXref>\n")
						}
						rec.WriteString("    </INSDSeq_xrefs>\n")
					}

					// end of record
					rec.WriteString("  </INSDSeq>\n")

					// send formatted record down channel
					txt := rec.String()
					out <- txt
					rec.Reset()
					// go to top of loop for next record
					break
				}

				// read next sequence line

				cols := strings.Fields(line)
				if len(cols) > 0 && !IsAllDigits(cols[0]) {
					fmt.Fprintf(os.Stderr, "ERROR: Unrecognized section "+cols[0]+"\n")
				}

				for _, str := range cols {

					if IsAllDigits(str) {
						continue
					}

					// append letters to sequence
					seq.WriteString(str)
				}

				// read next line and continue
				line = nextLine()
				row++

			}

			// continue to next record
		}
	}

	// launch single converter goroutine
	go convertGenBank(inp, out)

	return out
}

// GenBankRefIndex reads flatfiles and sends reference index XML records down a channel
func GenBankRefIndex(inp io.Reader, deStop, doStem bool) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, chanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "Unable to create GenBank converter channel\n")
		os.Exit(1)
	}

	const twelvespaces = "            "
	const twentyonespaces = "                     "

	var rec strings.Builder

	var arry []string

	indexGenBank := func(inp io.Reader, out chan<- string) {

		// close channel when all records have been sent
		defer close(out)

		scanr := bufio.NewScanner(inp)

		row := 0

		prevTitles := make(map[string]bool)

		nextLine := func() string {

			for scanr.Scan() {
				line := scanr.Text()
				if line == "" {
					continue
				}
				return line
			}
			return ""

		}

		for {
			arry = nil

			currTitle := ""
			firstAuth := ""
			division := ""
			dirSub := false

			// read first line of next record
			line := nextLine()
			if line == "" {
				break
			}

			row++

			for {
				if line == "" {
					break
				}
				if !strings.HasPrefix(line, "LOCUS") {
					// skip release file header information
					line = nextLine()
					row++
					continue
				}
				break
			}

			readContinuationLines := func(str string) string {

				for {
					// read next line
					line = nextLine()
					row++
					if !strings.HasPrefix(line, twelvespaces) {
						// if not continuation line, break out of loop
						break
					}
					// append subsequent line and continue with loop
					txt := strings.TrimPrefix(line, twelvespaces)
					str += " " + txt
				}

				str = CompressRunsOfSpaces(str)
				str = strings.TrimSpace(str)

				return str
			}

			writeOneElement := func(spaces, tag, value string) {

				rec.WriteString(spaces)
				rec.WriteString("<")
				rec.WriteString(tag)
				rec.WriteString(">")
				value = html.EscapeString(value)
				rec.WriteString(value)
				rec.WriteString("</")
				rec.WriteString(tag)
				rec.WriteString(">\n")
			}

			// each section will exit with the next line ready to process

			if strings.HasPrefix(line, "LOCUS") {

				// start of record

				// do not break if given artificial multi-line LOCUS
				str := readContinuationLines(line)

				cols := strings.Fields(str)
				ln := len(cols)
				if ln == 8 {
					division = cols[6]
				} else if ln == 7 {
					division = cols[5]
				}

				// read next line and continue - handled by readContinuationLines above
				// line = nextLine()
				// row++
			}

			if strings.HasPrefix(line, "DEFINITION") {
				readContinuationLines(line)
			}

			// record accession
			accn := ""

			if strings.HasPrefix(line, "ACCESSION") {

				txt := strings.TrimPrefix(line, "ACCESSION")
				str := readContinuationLines(txt)
				accessions := strings.Fields(str)
				ln := len(accessions)
				if ln > 0 {
					accn = accessions[0]
				}
			}

			if strings.HasPrefix(line, "VERSION") {

				cols := strings.Fields(line)
				if len(cols) == 2 || len(cols) == 3 {
					accn = cols[1]
				}

				// read next line and continue
				line = nextLine()
				row++
			}

			if strings.HasPrefix(line, "DBLINK") {
				readContinuationLines(line)
			}

			if strings.HasPrefix(line, "DBSOURCE") {
				readContinuationLines(line)
			}

			if strings.HasPrefix(line, "KEYWORDS") {
				readContinuationLines(line)
			}

			if strings.HasPrefix(line, "SOURCE") {
				readContinuationLines(line)
			}

			if strings.HasPrefix(line, "  ORGANISM") {

				line = nextLine()
				row++
				if strings.HasPrefix(line, twelvespaces) {
					readContinuationLines(line)
				}
			}

			// beginning of reference section
			for {
				if !strings.HasPrefix(line, "REFERENCE") {
					// exit out of reference section
					break
				}

				rec.Reset()

				rec.WriteString("  <CITATION>\n")

				if accn != "" {
					writeOneElement("    ", "ACCN", accn)
				}

				if division != "" {
					writeOneElement("    ", "DIV", division)
				}

				txt := strings.TrimPrefix(line, "REFERENCE")
				str := readContinuationLines(txt)
				str = CompressRunsOfSpaces(str)
				str = strings.TrimSpace(str)
				idx := strings.Index(str, "(")
				ref := ""
				if idx > 0 {
					ref = strings.TrimSpace(str[:idx])
				} else {
					ref = strings.TrimSpace(str)
				}
				// reference number
				writeOneElement("    ", "REF", ref)
				row++

				if strings.HasPrefix(line, "  AUTHORS") {

					faut := ""
					laut := ""
					count := 0

					txt := strings.TrimPrefix(line, "  AUTHORS")
					auths := readContinuationLines(txt)

					authors := strings.Split(auths, ", ")
					for _, auth := range authors {
						auth = strings.TrimSpace(auth)
						if auth == "" {
							continue
						}
						pair := strings.Split(auth, " and ")
						for _, name := range pair {

							// convert GenBank author to searchable form
							name = GenBankToMedlineAuthors(name)

							if faut == "" {
								faut = name
							}
							laut = name

							writeOneElement("    ", "AUTH", name)

							count++
						}
					}

					if faut != "" {
						writeOneElement("    ", "FAUT", faut)
						firstAuth = faut
					}
					if laut != "" {
						writeOneElement("    ", "LAUT", laut)
					}

					if count > 0 {
						// anum := strconv.Itoa(count)
						// writeOneElement("    ", "ANUM", anum)
					}
				}

				if strings.HasPrefix(line, "  CONSRTM") {

					txt := strings.TrimPrefix(line, "  CONSRTM")
					cons := readContinuationLines(txt)

					writeOneElement("    ", "CSRT", cons)
				}

				inPress := false

				if strings.HasPrefix(line, "  TITLE") {

					txt := strings.TrimPrefix(line, "  TITLE")
					titl := readContinuationLines(txt)

					if titl != "" {
						writeOneElement("    ", "TITL", titl)

						// ln := strconv.Itoa(len(titl))
						// writeOneElement("    ", "TLEN", ln)

						tnum := 0
						words := strings.FieldsFunc(titl, func(c rune) bool {
							return !unicode.IsLetter(c) && !unicode.IsDigit(c)
						})
						for _, item := range words {
							item = strings.ToLower(item)
							if deStop {
								// exclude stop words from count
								if IsStopWord(item) {
									continue
								}
							}
							if doStem {
								item = porter2.Stem(item)
								item = strings.TrimSpace(item)
							}
							if item == "" {
								continue
							}
							tnum++
						}
						// nm := strconv.Itoa(tnum)
						// writeOneElement("    ", "TNUM", nm)
					}

					if titl == "Direct Submission" {
						dirSub = true
					}

					currTitle = titl
				}

				if strings.HasPrefix(line, "  JOURNAL") {

					txt := strings.TrimPrefix(line, "  JOURNAL")
					jour := readContinuationLines(txt)

					if strings.HasSuffix(jour, " In press") {
						inPress = true
					}

					if dirSub {
						if len(jour) > 23 && strings.HasPrefix(jour, "Submitted (") {
							year := jour[18:22]
							year = strings.TrimSpace(year)
							if year != "" {
								writeOneElement("    ", "YEAR", year)
							}
						}
					} else if inPress {
						jour = strings.TrimSuffix(jour, " In press")
						if strings.HasSuffix(jour, ")") {
							idx0 := strings.LastIndex(jour, "(")
							if idx0 >= 0 {
								doi := ""
								year := jour[idx0:]
								jour := jour[:idx0]
								year = strings.TrimPrefix(year, "(")
								year = strings.TrimSuffix(year, ")")
								year = strings.TrimSpace(year)
								jour = strings.TrimSpace(jour)
								idxd := strings.Index(jour, ", doi:")
								if idxd >= 0 {
									doi = jour[idxd:]
									jour = jour[:idxd]
									doi = strings.TrimPrefix(doi, ", doi:")
									jour = strings.TrimSpace(jour)
									doi = strings.TrimSpace(doi)
								}
								// truncate journal at comma
								cmma := strings.Index(jour, ",")
								if cmma >= 0 {
									jour = jour[:cmma]
								}
								// truncate journal at first digit
								for i, r := range jour {
									if unicode.IsDigit(r) {
										jour = jour[:i]
										break
									}
								}
								jour = strings.Replace(jour, ".", " ", -1)
								jour = strings.TrimSpace(jour)
								jour = CompressRunsOfSpaces(jour)
								if jour != "" {
									writeOneElement("    ", "JOUR", jour)
								}
								if doi != "" {
									writeOneElement("    ", "DOI", doi)
								}
								if year != "" {
									writeOneElement("    ", "YEAR", year)
								}
							}
						}
					} else {
						journal := ""
						volume := ""
						issue := ""
						pages := ""
						year := ""
						lft, rgt := SplitInTwoLeft(jour, ",")
						if lft != "" && rgt != "" {
							if strings.HasSuffix(lft, ")") {
								idx1 := strings.LastIndex(lft, "(")
								if idx1 >= 0 {
									issue = lft[idx1:]
									lft = lft[:idx1]
									issue = strings.TrimPrefix(issue, "(")
									issue = strings.TrimSuffix(issue, ")")
									issue = strings.TrimSpace(issue)
									lft = strings.TrimSpace(lft)
								}
							}
							idx2 := strings.LastIndex(lft, " ")
							if idx2 >= 0 {
								volume = lft[idx2:]
								lft = lft[:idx2]
								volume = strings.TrimSpace(volume)
								lft = strings.TrimSpace(lft)
							}
							journal = lft
							if strings.HasSuffix(rgt, ")") {
								idx3 := strings.LastIndex(rgt, "(")
								if idx3 >= 0 {
									year = rgt[idx3:]
									rgt = rgt[:idx3]
									year = strings.TrimPrefix(year, "(")
									year = strings.TrimSuffix(year, ")")
									year = strings.TrimSpace(year)
									rgt = strings.TrimSpace(rgt)
								}
							}
							pages = rgt
							journal = strings.TrimSpace(journal)
							pages = strings.TrimSpace(pages)
							journal = strings.Replace(journal, ".", " ", -1)
							journal = strings.TrimSpace(journal)
							journal = CompressRunsOfSpaces(journal)
							if journal != "" {
								writeOneElement("    ", "JOUR", journal)
							}
							if volume != "" {
								writeOneElement("    ", "VOL", volume)
							}
							if issue != "" {
								writeOneElement("    ", "ISS", issue)
							}
							if pages != "" {
								writeOneElement("    ", "PAGE", pages)
							}
							if year != "" {
								writeOneElement("    ", "YEAR", year)
							}
						}
					}

				}

				pmid := ""
				if strings.HasPrefix(line, "   PUBMED") {

					txt := strings.TrimPrefix(line, "   PUBMED")
					pmid = readContinuationLines(txt)

					pmid = strings.TrimSpace(pmid)
				}

				stat := ""
				if dirSub {
					stat = "dirsub"
				} else if inPress {
					stat = "inpress"
				} else if pmid != "" {
					stat = "published"
				} else {
					stat = "unpub"
				}

				isDuplicate := false
				prv := firstAuth + "|" + currTitle + "|" + stat
				if prv != "" {
					if prevTitles[prv] {
						isDuplicate = true
					} else {
						prevTitles[prv] = true
					}
				}

				writeOneElement("    ", "STAT", stat)

				if strings.HasPrefix(line, "  MEDLINE") {
					// old MEDLINE uid not supported
					readContinuationLines(line)
				}

				if strings.HasPrefix(line, "  REMARK") {
					readContinuationLines(line)
				}

				if pmid != "" {
					// write PMID immediately before end of CITATION object
					writeOneElement("    ", "PMID", pmid)
				}

				rec.WriteString("  </CITATION>")

				txt = rec.String()

				if !isDuplicate {

					// record nonredundant reference
					arry = append(arry, txt)
				}

				// reset working buffer
				rec.Reset()
				// continue to next reference
			}
			// end of reference section

			if strings.HasPrefix(line, "COMMENT") {
				readContinuationLines(line)
			}

			if strings.HasPrefix(line, "PRIMARY") {
				readContinuationLines(line)
			}

			if strings.HasPrefix(line, "FEATURES") {

				line = nextLine()
				row++

				for {
					if !strings.HasPrefix(line, "     ") {
						// exit out of features section
						break
					}
					if len(line) < 22 {
						line = nextLine()
						row++
						continue
					}

					for {
						line = nextLine()
						row++
						if !strings.HasPrefix(line, twentyonespaces) {
							break
						}
						txt := strings.TrimPrefix(line, twentyonespaces)
						if strings.HasPrefix(txt, "/") {
							// if not continuation of location, break out of loop
							break
						}
						// append subsequent line and continue with loop
					}

					for {
						if !strings.HasPrefix(line, twentyonespaces) {
							// if not qualifier line, break out of loop
							break
						}
						txt := strings.TrimPrefix(line, twentyonespaces)
						if strings.HasPrefix(txt, "/") {
							// read new qualifier and start of value

							for {
								line = nextLine()
								row++
								if !strings.HasPrefix(line, twentyonespaces) {
									break
								}
								txt := strings.TrimPrefix(line, twentyonespaces)
								if strings.HasPrefix(txt, "/") {
									// if not continuation of qualifier, break out of loop
									break
								}
								// append subsequent line to value and continue with loop
							}
						}
					}
					// end of this feature
					// continue to next feature
				}
			}
			// TSA, TLS, WGS, or CONTIG lines may be next

			altName := ""

			if strings.HasPrefix(line, "TSA") ||
				strings.HasPrefix(line, "TLS") ||
				strings.HasPrefix(line, "WGS") {

				line = line[3:]
			}

			if strings.HasPrefix(line, "WGS_CONTIG") ||
				strings.HasPrefix(line, "WGS_SCAFLD") {

				line = line[10:]
			}

			if altName != "" {

				for {
					// read next line
					line = nextLine()
					row++
					if !strings.HasPrefix(line, twelvespaces) {
						// if not continuation of contig, break out of loop
						break
					}
					// append subsequent line and continue with loop
				}
			}

			if strings.HasPrefix(line, "CONTIG") {

				for {
					// read next line
					line = nextLine()
					row++
					if !strings.HasPrefix(line, twelvespaces) {
						// if not continuation of contig, break out of loop
						break
					}
					// append subsequent line and continue with loop
				}
			}

			if strings.HasPrefix(line, "BASE COUNT") {
				readContinuationLines(line)
				// not supported
			}

			if strings.HasPrefix(line, "ORIGIN") {

				line = nextLine()
				row++
			}

			// remainder should be sequence

			for line != "" {

				if strings.HasPrefix(line, "//") {

					for _, txt := range arry {
						// send formatted record down channel
						out <- txt
					}

					// go to top of loop for next record
					break
				}

				// read next line and continue
				line = nextLine()
				row++

			}

			// continue to next record
		}
	}

	// launch single indexer goroutine
	go indexGenBank(inp, out)

	return out
}

// CreateCitMatchers reads CITATION XML and returns matching PMIDs
func CreateCitMatchers(inp <-chan XMLRecord, options []string, deStop, doStem bool) <-chan XMLRecord {

	if inp == nil {
		return nil
	}

	out := make(chan XMLRecord, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create citmatch channel\n")
		os.Exit(1)
	}

	// obtain path from environment variable
	base := os.Getenv("EDIRECT_PUBMED_MASTER")
	if base != "" {
		if !strings.HasSuffix(base, "/") {
			base += "/"
		}
	}

	postingsBase := base + "Postings"
	archiveBase := base + "Archive"

	// check to make sure local archive is mounted
	_, err := os.Stat(archiveBase)
	if err != nil && os.IsNotExist(err) {
		fmt.Fprintf(os.Stderr, "\nERROR: Local archive and search index is not mounted\n\n")
		os.Exit(1)
	}

	debug := false
	strict := false

	for _, rgs := range options {
		opts := strings.Split(rgs, ",")
		for _, opt := range opts {
			if opt == "" {
				continue
			}
			switch opt {
			case "debug":
				debug = true
			case "strict":
				strict = true
			default:
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized -options choice '%s'\n", opt)
				os.Exit(1)
			}
		}
	}

	// gbCitMatch reads partitioned XML from channel and looks up candidate PMIDs
	gbCitMatch := func(wg *sync.WaitGroup, inp <-chan XMLRecord, out chan<- XMLRecord) {

		// report when this matcher has no more records to process
		defer wg.Done()

		sortedWords := func(titl string) []string {

			temp := make(map[string]bool)

			titl = strings.ToLower(titl)

			// break phrases into individual words
			words := strings.FieldsFunc(titl, func(c rune) bool {
				return !unicode.IsLetter(c) && !unicode.IsDigit(c)
			})

			for _, item := range words {
				if IsStopWord(item) {
					continue
				}
				temp[item] = true
			}

			var arry []string

			for item := range temp {
				if temp[item] {
					arry = append(arry, item)
				}
			}

			sort.Slice(arry, func(i, j int) bool { return arry[i] < arry[j] })

			return arry
		}

		fetchRecord := func(file string) string {

			var buf bytes.Buffer

			var arry [132]rune
			trie := MakeArchiveTrie(file, arry)

			if file == "" || trie == "" {
				return ""
			}

			zipp := true
			sfx := ".xml"
			if zipp {
				sfx += ".gz"
			}

			fpath := path.Join(archiveBase, trie, file+sfx)
			if fpath == "" {
				return ""
			}

			iszip := zipp

			inFile, err := os.Open(fpath)

			// if failed to find ".xml" or ".e2x" file, try ".xml.gz" or ".e2x.gz" without requiring -gzip
			if err != nil && os.IsNotExist(err) && !zipp {
				iszip = true
				fpath := path.Join(archiveBase, trie, file+sfx+".gz")
				if fpath == "" {
					return ""
				}
				inFile, err = os.Open(fpath)
			}
			if err != nil {
				msg := err.Error()
				if !strings.HasSuffix(msg, "no such file or directory") && !strings.HasSuffix(msg, "cannot find the path specified.") {
					fmt.Fprintf(os.Stderr, "%s\n", msg)
				}
				return ""
			}

			defer inFile.Close()

			brd := bufio.NewReader(inFile)

			if iszip {

				// using parallel pgzip for large release files
				zpr, err := pgzip.NewReader(brd)

				defer zpr.Close()

				if err == nil {
					// copy and decompress cached file contents
					buf.ReadFrom(zpr)
				}

			} else {

				// copy cached file contents
				buf.ReadFrom(brd)
			}

			str := buf.String()

			return str
		}

		getTitle := func(uid int32) string {

			if uid < 1 {
				return ""
			}

			val := strconv.Itoa(int(uid))
			pma := fetchRecord(val)
			pma = strings.TrimSpace(pma)
			if pma == "" {
				return ""
			}

			refFields := make(map[string]string)

			StreamValues(pma[:], "PubmedArticle", func(tag, attr, content string) { refFields[tag] = content })

			titl := refFields["ArticleTitle"]

			return titl
		}

		uniqueWords := func(strs []string) []string {

			rs := make([]string, 0, len(strs))
			mp := make(map[string]bool)

			for _, val := range strs {
				_, ok := mp[val]
				if !ok {
					mp[val] = true
					rs = append(rs, val)
				}
			}

			return rs
		}

		intersectWords := func(a, b []string) int {

			temp := make(map[string]bool)
			num := 0

			for _, item := range a {
				temp[item] = true
			}

			for _, item := range b {
				if temp[item] {
					num++
				}
			}

			return num
		}

		unionWords := func(a, b []string) int {

			temp := make(map[string]bool)
			num := 0

			for _, item := range a {
				if !temp[item] {
					num++
				}
				temp[item] = true
			}

			for _, item := range b {
				if !temp[item] {
					num++
				}
				temp[item] = true
			}

			return num
		}

		cleanTitle := func(str string) string {

			if str == "" {
				return str
			}

			str = FixMisusedLetters(str, true, false, true)
			str = TransformAccents(str, false, false)
			if HasUnicodeMarkup(str) {
				str = RepairUnicodeMarkup(str, SPACE)
			}
			if HasAngleBracket(str) {
				str = RepairTableMarkup(str, SPACE)
				// if wrp // str = EncodeAngleBrackets(str)
			}
			if HasAmpOrNotASCII(str) {
				str = html.UnescapeString(str)
			}
			if HasBadSpace(str) {
				str = CleanupBadSpaces(str)
			}
			if HasAdjacentSpaces(str) {
				str = CompressRunsOfSpaces(str)
			}
			if HasExtraSpaces(str) {
				str = RemoveExtraSpaces(str)
			}

			return str
		}

		cleanAuthor := func(str string) string {

			if str == "" {
				return str
			}

			str = FixMisusedLetters(str, true, true, false)
			str = TransformAccents(str, false, false)

			// convert numeric encoding to apostrophe
			str = strings.Replace(str, "&#39;", "'", -1)
			// remove space following apostrophe
			str = strings.Replace(str, "' ", "'", -1)
			// then remove apostrophe to match indexing logic
			str = strings.Replace(str, "'", "", -1)

			if HasUnicodeMarkup(str) {
				str = RepairUnicodeMarkup(str, SPACE)
			}
			if HasAngleBracket(str) {
				str = RepairTableMarkup(str, SPACE)
				// if wrp // str = EncodeAngleBrackets(str)
			}
			if HasBadSpace(str) {
				str = CleanupBadSpaces(str)
			}
			str = strings.ToLower(str)
			str = strings.Replace(str, "(", " ", -1)
			str = strings.Replace(str, ")", " ", -1)
			if HasAdjacentSpaces(str) {
				str = CompressRunsOfSpaces(str)
			}
			if HasExtraSpaces(str) {
				str = RemoveExtraSpaces(str)
			}

			// remove jr and sr suffix
			if strings.Index(str, " ") != strings.LastIndex(str, " ") {
				str = strings.TrimSuffix(str, " jr")
				str = strings.TrimSuffix(str, " sr")
			}

			// truncate to single initial
			pos := strings.LastIndex(str, " ")
			if pos > 0 && len(str) > pos+2 {
				str = str[:pos+2]
			}

			return str
		}

		cleanParens := func(str string) string {

			if str == "" {
				return str
			}

			// remove parentheses and other special characters
			str = strings.ToLower(str)

			str = strings.Replace(str, "(", " ", -1)
			str = strings.Replace(str, ")", " ", -1)

			if HasHyphenOrApostrophe(str) {
				str = FixSpecialCases(str)
			}

			str = strings.Replace(str, "_", " ", -1)
			str = strings.Replace(str, "-", " ", -1)
			str = strings.Replace(str, "+", " ", -1)
			str = strings.Replace(str, "~", " ", -1)

			str = CompressRunsOfSpaces(str)
			str = strings.TrimSpace(str)

			return str
		}

		// look for closest match to actual title among candidate PMIDs
		jaccard := func(titl string, ids []int32) int32 {

			if len(ids) < 1 {
				return 0
			}

			titl = cleanTitle(titl)
			titl = cleanParens(titl)
			titleWords := sortedWords(titl)
			titleWords = uniqueWords(titleWords)

			bestScore := 0
			bestPMID := int32(0)

			if debug {
				fmt.Fprintf(os.Stderr, "             %s\n", titl)
			}

			for _, uid := range ids {
				ttl := getTitle(uid)
				if ttl != "" {
					ttl = cleanTitle(ttl)
					ttl = cleanParens(ttl)
					ttlWords := sortedWords(ttl)
					ttlWords = uniqueWords(ttlWords)

					intrs := intersectWords(titleWords, ttlWords)
					union := unionWords(titleWords, ttlWords)
					score := intrs * 100 / union

					if debug {
						fmt.Fprintf(os.Stderr, "%8d %3d %s\n", uid, score, ttl)
					}

					// highest score should prefer original paper over errata and corrigenda
					if score > bestScore {
						bestScore = score
						bestPMID = uid
					}
				}
			}

			// require score of at least 60 to filter out false positives
			if bestScore < 60 {
				return 0
			}

			return bestPMID
		}

		intersectMatches := func(a, b []int32) []int32 {

			temp := make(map[int32]bool)
			var res []int32

			for _, item := range a {
				temp[item] = true
			}

			for _, item := range b {
				if temp[item] {
					res = append(res, item)
				}
			}

			return res
		}

		wordPairs := func(titl string) []string {

			var arry []string

			titl = strings.ToLower(titl)

			// break phrases into individual words
			words := strings.FieldsFunc(titl, func(c rune) bool {
				return !unicode.IsLetter(c) && !unicode.IsDigit(c)
			})

			// word pairs (or isolated singletons) separated by stop words
			if len(words) > 1 {
				past := ""
				run := 0
				for _, item := range words {
					if IsStopWord(item) {
						if run == 1 && past != "" {
							arry = append(arry, past)
						}
						past = ""
						run = 0
						continue
					}
					if item == "" {
						past = ""
						continue
					}
					if past != "" {
						arry = append(arry, past+" "+item)
					}
					past = item
					run++
				}
				if run == 1 && past != "" {
					arry = append(arry, past)
				}
			}

			return arry
		}

		/*
			wordRuns := func(titl string) []string {

				var arry []string

				titl = strings.ToLower(titl)

				words := strings.FieldsFunc(titl, func(c rune) bool {
					return !unicode.IsLetter(c) && !unicode.IsDigit(c)
				})

				// phrases (or isolated singletons) separated by stop words
				if len(words) > 1 {
					i := 0
					j := 0
					item := ""
					for j, item = range words {
						if IsStopWord(item) {
							if j > i {
								term := strings.Join(words[i:j], " ")
								arry = append(arry, term)
							}
							i = j + 1
						}
					}
					if j > i {
						term := strings.Join(words[i:j], " ")
						arry = append(arry, term)
					}
				}

				return arry
			}
		*/

		searchTitleParts := func(terms []string) []int32 {

			var candidates []int32

			if len(terms) < 1 {
				return candidates
			}

			// histogram of counts for PMIDs matching one or more word pair queries
			pmatch := make(map[int32]int)

			for _, item := range terms {
				arry := ProcessQuery(postingsBase, item+" [PAIR]", false, false, false, deStop)
				for _, uid := range arry {
					val := pmatch[uid]
					val++
					pmatch[uid] = val
				}
			}

			// find maximum number of adjacent overlapping word pair matches
			max := 0
			for _, vl := range pmatch {
				if vl > max {
					max = vl
				}
			}

			// require at least 4 matching pairs to avoid false positives
			if max < 4 {
				return candidates
			}

			// collect up to 25 PMIDs with maximum adjacent overlapping word pair count
			num := 0
			for ky, vl := range pmatch {
				if vl == max {
					candidates = append(candidates, ky)
					num++
					if num >= 25 {
						break
					}
				}
			}

			// histogram showing distribution of partial matches
			if debug {
				if num > 0 {
					top, mid, low := 0, 0, 0
					for _, vl := range pmatch {
						if vl == max {
							top++
						} else if vl == max-1 {
							mid++
						} else if vl == max-2 {
							low++
						}
					}
					sep := "matches: "
					if top > 0 {
						fmt.Fprintf(os.Stderr, "%s%d @ [%d]", sep, top, max)
						sep = ", "
					}
					if mid > 0 {
						fmt.Fprintf(os.Stderr, "%s%d @ [%d]", sep, mid, max-1)
						sep = ", "
					}
					if low > 0 {
						fmt.Fprintf(os.Stderr, "%s%d @ [%d]", sep, low, max-2)
					}
					fmt.Fprintf(os.Stderr, "\n")
				}

			}

			return candidates
		}

		matchByTitle := func(titl string) ([]int32, string) {

			var byTitle []int32

			titl = cleanTitle(titl)
			titl = cleanParens(titl)
			if debug {
				fmt.Fprintf(os.Stderr, "title:   %s\n", titl)
			}

			if titl == "" {
				return byTitle, "missing title"
			}

			// adjacent overlapping word pairs, plus single words between stop words or at ends
			pairs := wordPairs(titl)
			if len(pairs) < 1 {
				return byTitle, "title has no matching word pairs in search index"
			}

			byTitle = searchTitleParts(pairs)
			if len(byTitle) < 1 {
				return byTitle, "no candidates for matching title"
			}

			return byTitle, ""
		}

		matchByAuthor := func(faut, laut string) ([]int32, string, string) {

			var byAuthor []int32

			faut = cleanAuthor(faut)
			faut = cleanParens(faut)

			laut = cleanAuthor(laut)
			laut = cleanParens(laut)

			if faut == "" && laut == "" {
				return byAuthor, "empty authors", ""
			}

			query := faut
			if strings.Index(faut, " ") < 0 {
				// if just last name, space plus asterisk to wildcard on initials
				query += " "
			}
			// otherwise, if space between last name and initials, immediate asterisk for term truncation
			if strict {
				query += "* [FAUT]"
			} else {
				query += "* [AUTH]"
			}
			if faut != laut && laut != "" {
				query += " OR " + laut
				if strings.Index(laut, " ") < 0 {
					query += " "
				}
				if strict {
					query += "* [LAUT]"
				} else {
					query += "* [AUTH]"
				}
			}
			if debug {
				fmt.Fprintf(os.Stderr, "authors: %s\n", query)
			}

			names := faut
			if faut != laut && laut != "" {
				names += " OR " + laut
			}

			// find PMIDs indexed under first or last author, use wildcard after truncating to single initial
			byAuthor = ProcessQuery(postingsBase, query, false, false, false, deStop)
			if len(byAuthor) < 1 {
				return byAuthor, "unrecognized author '" + names + "'", names
			}

			return byAuthor, "", names
		}

		matchByJournal := func(jour string) ([]int32, string, string) {

			var byJournal []int32

			jour = cleanParens(jour)

			if jour == "" {
				return byJournal, "empty journal", ""
			}

			query := jour + " [JOUR]"
			if debug {
				fmt.Fprintf(os.Stderr, "journal: %s\n", query)
			}

			byJournal = ProcessQuery(postingsBase, query, false, false, false, deStop)
			if len(byJournal) < 1 {
				return byJournal, "unrecognized journal '" + jour + "'", jour
			}

			return byJournal, "", jour
		}

		matchByYear := func(year string) ([]int32, string, string) {

			var byYear []int32

			if year == "" {
				return byYear, "", ""
			}

			yr, err := strconv.Atoi(year)
			if err != nil {
				return byYear, "unrecognized year '" + year + "'", year
			}
			lst := strconv.Itoa(yr - 1)
			nxt := strconv.Itoa(yr + 1)
			if strict {
				lst = strconv.Itoa(yr)
			}

			query := lst + ":" + nxt + " [YEAR]"
			if debug {
				fmt.Fprintf(os.Stderr, "year:    %s\n", query)
			}

			span := lst + ":" + nxt

			byYear = ProcessQuery(postingsBase, query, false, false, false, deStop)
			if len(byYear) < 1 {
				return byYear, "unrecognized year range '" + span + "'", span
			}

			return byYear, "", span
		}

		// citFind returns PMID and optional message containing reason for failure
		citFind := func(citFields map[string]string) (int32, string) {

			if citFields == nil {
				return 0, "map missing"
			}

			note := ""
			between := ""

			// initial candidates based on most matches to overlapping word pairs in title

			titl := citFields["TITL"]

			byTitle, reasonT := matchByTitle(titl)
			if reasonT != "" {
				return 0, reasonT
			}

			// prepare postings subsets to filter candidates by author, journal, and year

			faut := citFields["FAUT"]
			laut := citFields["LAUT"]

			byAuthor, reasonA, labelA := matchByAuthor(faut, laut)
			if reasonA != "" {
				if strict {
					return 0, reasonA
				}
				note += between + reasonA
				between = ", "
			}

			jour := citFields["JOUR"]

			byJournal, reasonJ, labelJ := matchByJournal(jour)
			if reasonJ != "" {
				if strict {
					return 0, reasonJ
				}
				note += between + reasonJ
				between = ", "
			}

			year := citFields["YEAR"]

			byYear, reasonY, labelY := matchByYear(year)
			if reasonY != "" {
				if strict {
					return 0, reasonY
				}
				note += between + reasonY
				between = ", "
			}

			// interesections

			working := byTitle

			// restrict by author name
			if len(byAuthor) > 0 {
				temp := intersectMatches(working, byAuthor)
				if len(temp) < 1 {
					if strict {
						return 0, "author does not match title"
					}
					note += between + "title does not match author '" + labelA + "'"
					return 0, note + ", exiting"
				}
				working = temp
			} else if strict {
				return 0, "no author match"
			}

			// restrict by journal name, but ignore if no match
			if len(byJournal) > 0 {
				temp := intersectMatches(working, byJournal)
				if len(temp) < 1 {
					if strict {
						return 0, "journal does not match title"
					}
					note += between + "title does not match journal '" + labelJ + "'"
					between = ", "
				} else {
					working = temp
				}
			} else if strict {
				return 0, "no journal match"
			}

			// restrict by year range, but ignore if no match
			if len(byYear) > 0 {
				temp := intersectMatches(working, byYear)
				if len(temp) < 1 {
					if strict {
						return 0, "year range does not match title"
					}
					note += between + "title does not match year range '" + labelY + "'"
					between = ", "
				} else {
					working = temp
				}
			} else if strict {
				return 0, "no year match"
			}

			if len(working) < 1 {
				return 0, "match not found"
			}

			// get best matching candidate
			pmid := jaccard(titl, working)
			if pmid != 0 {
				return pmid, note
			}

			note += between + "jaccard failed"
			return pmid, note
		}

		citMatch := func(text string) (int32, string, string) {

			citFields := make(map[string]string)

			// stream tokens, collecting XML values in map
			StreamValues(text[:], "CITATION", func(tag, attr, content string) { citFields[tag] = content })

			orig := citFields["ORIG"]

			pmid, reason := citFind(citFields)

			return pmid, orig, reason
		}

		re, _ := regexp.Compile(">[ \n\r\t]+<")

		// read partitioned XML from producer channel
		for ext := range inp {

			text := ext.Text

			// rename original PMID
			text = strings.Replace(text, "<PMID>", "<ORIG>", -1)
			text = strings.Replace(text, "</PMID>", "</ORIG>", -1)

			pmid, orig, reason := citMatch(text[:])

			if debug {
				fmt.Fprintf(os.Stderr, "pmid %d, orig %s, reason %s\n", pmid, orig, reason)
			}

			/*
				if pmid == 0 && orig == "" && reason == "" {
					// if no match, and no ORIG value, send empty placeholder for unshuffler
					out <- XMLRecord{Index: ext.Index}
					continue
				}
			*/

			lft, rgt := SplitInTwoLeft(text, "</CITATION")
			if lft != "" && rgt != "" {
				pm := ""
				if pmid > 0 {
					pm = "<PMID>" + strconv.Itoa(int(pmid)) + "</PMID>"
				}
				nt := ""
				if reason != "" {
					nt = "<NOTE>" + reason + "</NOTE>"
				}
				text = lft + pm + nt + "</CITATION" + rgt
			}

			if re != nil {
				text = re.ReplaceAllString(text, "><")
			}

			// send record with newly-matched PMID
			out <- XMLRecord{Index: ext.Index, Text: text}
		}
	}

	var wg sync.WaitGroup

	// launch multiple citmatch goroutines
	for i := 0; i < NumServe(); i++ {
		wg.Add(1)
		go gbCitMatch(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all citation matchers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}
