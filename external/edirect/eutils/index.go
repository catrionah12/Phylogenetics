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
// File Name:  index.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package eutils

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"github.com/klauspost/pgzip"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"runtime/debug"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

// ENTREZ2INDEX COMMAND GENERATOR

// ProcessE2Index generates extraction commands to create input for Entrez2Index
func ProcessE2Index(tform string, isPipe, xtras bool) []string {

	currentYear := strconv.Itoa(time.Now().Year())

	var acc []string

	if !isPipe {
		if !deStop {
			acc = append(acc, "-stops")
		}
		if doStem {
			acc = append(acc, "-stems")
		}
	}

	acc = append(acc, "-set", "IdxDocumentSet", "-rec", "IdxDocument")
	acc = append(acc, "-pattern", "PubmedArticle")
	acc = append(acc, "-wrp", "IdxUid", "-element", "MedlineCitation/PMID", "-clr", "-rst", "-tab", "")

	acc = append(acc, "-group", "PubmedArticle", "-pkg", "IdxSearchFields")

	// index YEAR of publication
	acc = append(acc, "-block", "PubmedArticle", "-wrp", "YEAR", "-year", "PubDate/*")

	// citation fields - JOUR, VOL, ISS, PAGE, and LANG

	acc = append(acc, "-block", "MedlineJournalInfo", "-wrp", "JOUR", "-element", "MedlineTA", "NlmUniqueID", "ISSNLinking")
	acc = append(acc, "-block", "Article/Journal", "-wrp", "JOUR", "-element", "Title", "ISOAbbreviation", "ISSN")
	acc = append(acc, "-wrp", "VOL", "-element", "Volume", "-wrp", "ISS", "-element", "Issue")
	acc = append(acc, "-block", "Article/Pagination", "-wrp", "PAGE", "-page", "MedlinePgn")
	acc = append(acc, "-block", "Article/Language", "-wrp", "LANG", "-element", "Language")

	// author fields - ANUM, AUTH, FAUT, LAUT, CSRT, and INVR

	// only count human authors, not consortia
	acc = append(acc, "-block", "AuthorList", "-wrp", "ANUM", "-num", "Author/LastName")
	// use -position to get first author
	acc = append(acc, "-block", "AuthorList/Author", "-position", "first")
	acc = append(acc, "-wrp", "FAUT", "-sep", " ", "-author", "LastName,Initials")
	// expect consortium to be last in the author list, so explore each author, and if last name is present,
	// overwrite the LAST variable with the current person's name
	acc = append(acc, "-block", "AuthorList/Author", "-if", "LastName", "-sep", " ", "-LAST", "LastName,Initials")
	// then explore on (single instance) PubmedArticle to print one copy of variable containing the last author's name
	acc = append(acc, "-block", "PubmedArticle", "-if", "&LAST", "-wrp", "LAUT", "-author", "&LAST")
	// separate field for consortia
	acc = append(acc, "-block", "AuthorList/Author", "-wrp", "CSRT", "-prose", "CollectiveName")
	// now get all human authors and investigators
	acc = append(acc, "-block", "AuthorList/Author", "-wrp", "AUTH", "-sep", " ", "-author", "LastName,Initials")
	acc = append(acc, "-block", "InvestigatorList/Investigator", "-wrp", "INVR", "-sep", " ", "-author", "LastName,Initials")

	// title and abstract fields - TIAB, TITL, and PAIR

	// positional indices for TITL and TIAB fields
	acc = append(acc, "-block", "PubmedArticle", "-article", "ArticleTitle")
	acc = append(acc, "-block", "PubmedArticle", "-indices", "ArticleTitle,Abstract/AbstractText")
	// overlapping adjacent word pairs (or isolated singletons) separated by stop words
	acc = append(acc, "-block", "PubmedArticle", "-wrp", "PAIR", "-pairx", "ArticleTitle")

	// property fields - PROP

	acc = append(acc, "-block", "CommentsCorrections", "-wrp", "PROP", "-prop", "@RefType")
	acc = append(acc, "-block", "PublicationStatus", "-wrp", "PROP", "-prop", "PublicationStatus")
	acc = append(acc, "-block", "PublicationType", "-wrp", "PROP", "-element", "PublicationType")
	acc = append(acc, "-block", "Abstract", "-if", "AbstractText", "-wrp", "PROP", "-lbl", "Has Abstract")
	// dates
	acc = append(acc, "-block", "Journal", "-if", "MedlineDate", "-wrp", "PROP", "-lbl", "Medline Date")
	acc = append(acc, "-subset", "MedlineDate", "-if", "%MedlineDate", "-lt", "4", "-wrp", "PROP", "-lbl", "Bad Date")
	acc = append(acc, "-block", "PubMedPubDate", "-if", "%Year", "-lt", "4", "-wrp", "PROP", "-lbl", "Bad Date")
	acc = append(acc, "-block", "PubDate", "-if", "Year", "-gt", currentYear, "-wrp", "PROP", "-lbl", "Future Date")
	acc = append(acc, "-block", "PubMedPubDate", "-if", "PubMedPubDate@PubStatus", "-is-not", "pmc-release")
	acc = append(acc, "-and", "Year", "-gt", currentYear, "-wrp", "PROP", "-lbl", "Future Date")
	// version
	acc = append(acc, "-block", "MedlineCitation", "-if", "PMID@Version", "-gt", "1", "-wrp", "PROP", "-lbl", "Versioned")

	// if Extras/meshtree.txt is available, index CODE, TREE, and SUBH fields, and MESH for term list
	if tform != "" {
		acc = append(acc, "-block", "PubmedArticle", "-meshcode")
		acc = append(acc, "MeshHeading/DescriptorName@UI,Chemical/NameOfSubstance@UI,SupplMeshName@UI")
		acc = append(acc, "-block", "MeshHeading/QualifierName", "-wrp", "SUBH", "-element", "QualifierName")
		// only populating MESH for live term list, since query will redirect to wildcard on TREE
		acc = append(acc, "-block", "MeshHeading/DescriptorName", "-wrp", "MESH", "-element", "DescriptorName")
	}

	// optionally index MAJR fields - not (yet) expanded to lower elements in hierarchy tree
	if xtras {
		acc = append(acc, "-block", "MeshHeading/DescriptorName", "-if", "@MajorTopicYN")
		acc = append(acc, "-equals", "Y", "-wrp", "MAJR", "-element", "DescriptorName")
	}

	return acc
}

// UPDATE CACHED INVERTED FILES IN LOCAL ARCHIVE FOLDERS

// remove cached inverted index components with:
//
//   find /Volumes/cachet/Archive -name \*.inv.gz -type f -delete

// examineFolder collects two-digit subdirectories, xml files, and inv files
func examineFolder(base, path string) ([]string, []string, []string) {

	dir := filepath.Join(base, path)

	contents, err := os.ReadDir(dir)
	if err != nil {
		return nil, nil, nil
	}

	isTwoDigits := func(str string) bool {

		if len(str) != 2 {
			return false
		}

		ch := str[0]
		if ch < '0' || ch > '9' {
			return false
		}

		ch = str[1]
		if ch < '0' || ch > '9' {
			return false
		}

		return true
	}

	var dirs []string
	var xmls []string
	var invs []string

	for _, item := range contents {
		name := item.Name()
		if name == "" {
			continue
		}
		if item.IsDir() {
			if isTwoDigits(name) {
				dirs = append(dirs, name)
			}
		} else if strings.HasSuffix(name, ".xml.gz") {
			xmls = append(xmls, name)
		} else if strings.HasSuffix(name, ".inv.gz") {
			invs = append(invs, name)
		}
	}

	return dirs, xmls, invs
}

// gzFileToString reads selected gzipped file, uncompress, save contents as string
func gzFileToString(fpath string, parallelGz bool) string {

	file, err := os.Open(fpath)
	if err != nil {
		return ""
	}
	defer file.Close()

	var rdr io.Reader

	if parallelGz {

		// for large files, use parallel pgzip
		pgz, err := pgzip.NewReader(file)
		if err != nil {
			return ""
		}
		defer pgz.Close()
		rdr = pgz

	} else {

		// for small files, use regular gzip
		gz, err := gzip.NewReader(file)
		if err != nil {
			return ""
		}
		defer gz.Close()
		rdr = gz
	}

	byt, err := io.ReadAll(rdr)
	if err != nil {
		return ""
	}

	str := string(byt)
	if str == "" {
		return ""
	}

	if !strings.HasSuffix(str, "\n") {
		str += "\n"
	}

	return str
}

// readGzFiles reads and uncompresses a set of gzip-compressed files
func readGzFiles(names []string, parallelGz bool, numServ int) <-chan string {

	out := make(chan string, ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create read gz files channel\n")
		os.Exit(1)
	}

	readCompressedFiles := func(wg *sync.WaitGroup, parallelGz bool, inp <-chan string, out chan<- string) {

		defer wg.Done()

		for fpath := range inp {

			str := gzFileToString(fpath, parallelGz)

			if str == "" {
				continue
			}

			out <- str

			runtime.Gosched()
		}
	}

	inp := SliceToChan(names)
	if inp == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create file name channel\n")
		os.Exit(1)
	}

	var wg sync.WaitGroup

	// launch multiple file reader goroutines
	for i := 0; i < /* NumServe() */ numServ; i++ {
		wg.Add(1)
		go readCompressedFiles(&wg, parallelGz, inp, out)
	}

	// launch separate anonymous goroutine to wait until all file readers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// e2IndexConsumer callbacks have access to application-specific data as closures
type e2IndexConsumer func(inp <-chan XMLRecord) <-chan XMLRecord

// IndexAndInvert processes local archive, updating leaf-folder cached .inv.gz indexed and inverted file
func IndexAndInvert(archiveBase string, csmr e2IndexConsumer) <-chan string {

	if csmr == nil {
		return nil
	}

	if archiveBase == "" {

		// if not passed as an argument, obtain archive base path from environment variable
		base := os.Getenv("EDIRECT_PUBMED_MASTER")
		if base != "" {
			if !strings.HasSuffix(base, "/") {
				base += "/"
			}

			archiveBase = base + "Archive"
		}

		if archiveBase == "" {
			fmt.Fprintf(os.Stderr, "\nERROR: EDIRECT_PUBMED_MASTER environment variable is not set\n\n")
			os.Exit(1)
		}
	}

	// check to make sure local archive is mounted
	_, err := os.Stat(archiveBase)
	if err != nil && os.IsNotExist(err) {
		fmt.Fprintf(os.Stderr, "\nERROR: Local archive and search index is not mounted\n\n")
		os.Exit(1)
	}

	stringsToGzFile := func(base, path, name string, proc func(*bufio.Writer)) {

		if proc == nil {
			return
		}

		fpath := filepath.Join(base, path, name)

		// overwrites and truncates existing file
		fl, err := os.Create(fpath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return
		}

		// for small files, use regular gzip
		zpr, err := gzip.NewWriterLevel(fl, gzip.BestSpeed)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create compressor\n")
			os.Exit(1)
		}

		wrtr := bufio.NewWriter(zpr)

		// write contents
		proc(wrtr)

		err = wrtr.Flush()
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return
		}

		err = zpr.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return
		}

		// fl.Sync()

		err = fl.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return
		}
	}

	visitArchiveFolders := func(archiveBase string) <-chan string {

		out := make(chan string, ChanDepth())
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create archive explorer channel\n")
			os.Exit(1)
		}

		// recursive definition
		var visitSubFolders func(base, path string, out chan<- string)

		// visitSubFolders recurses to leaf directories
		visitSubFolders = func(base, path string, out chan<- string) {

			dirs, xmls, invs := examineFolder(base, path)

			// recursively explore subdirectories
			if dirs != nil {
				for _, dr := range dirs {
					sub := filepath.Join(path, dr)
					visitSubFolders(base, sub, out)
				}
				return
			}

			// looking for directory with at least one *.xml.gz archive file but no combined .inv.gz inverted index file
			if invs != nil || xmls == nil {
				return
			}

			// record directory that needs updating
			out <- path
		}

		visitArchiveSubset := func(wg *sync.WaitGroup, base, path string, out chan<- string) {

			defer wg.Done()

			visitSubFolders(base, path, out)
		}

		var wg sync.WaitGroup

		// launch multiple archive visitor goroutines
		dirs, _, _ := examineFolder(archiveBase, "")
		// reverse the order to put more recent folders first
		sort.Slice(dirs, func(i, j int) bool { return dirs[i] > dirs[j] })
		for _, dr := range dirs {
			wg.Add(1)
			go visitArchiveSubset(&wg, archiveBase, dr, out)
		}

		// launch separate anonymous goroutine to wait until all archive visitors are done
		go func() {
			wg.Wait()
			close(out)
		}()

		return out
	}

	e2indexAndInvert := func(archiveBase string, csmr e2IndexConsumer, inp <-chan string) <-chan string {

		if csmr == nil || inp == nil {
			return nil
		}

		out := make(chan string, ChanDepth())
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create e2invert channel\n")
			os.Exit(1)
		}

		// mutex to protect access to rollingCount and rollingColumn variables
		var vlock sync.Mutex

		rollingCount := 0
		rollingColumn := 0

		countSuccess := func() {

			vlock.Lock()

			rollingCount++
			// show a dot every 200 .inv files, generated from up to 20000 .xml files
			if rollingCount >= 200 {
				rollingCount = 0
				// print dot (progress monitor)
				fmt.Fprintf(os.Stderr, ".")
				rollingColumn++
				if rollingColumn > 49 {
					// print newline after 50 dots
					fmt.Fprintf(os.Stderr, "\n")
					rollingColumn = 0
				}
			}

			vlock.Unlock()
		}

		indexerInverter := func(num int, wg *sync.WaitGroup, base string, csmr e2IndexConsumer, inp <-chan string, out chan<- string) {

			defer wg.Done()

			for path := range inp {

				dirs, xmls, invs := examineFolder(base, path)
				// input was already filtered for folders with just .xml.gz files
				if dirs != nil || xmls == nil || invs != nil {
					continue
				}

				// collect PubMed files in single archive directory
				var filenames []string
				for _, file := range xmls {
					fpath := filepath.Join(base, path, file)
					filenames = append(filenames, fpath)
				}
				if len(filenames) < 1 {
					continue
				}

				var recs []string
				for _, fpath := range filenames {
					// for small files, use regular gzip
					txt := gzFileToString(fpath, false)
					recs = append(recs, txt)
				}
				if len(recs) < 1 {
					continue
				}

				dfls := strings.Join(recs, "\n")
				rrdr := strings.NewReader(dfls)
				if dfls == "" || rrdr == nil {
					return
				}

				// pass (xml) files through -e2index
				trdr := CreateXMLStreamer(rrdr)
				xmlq := CreateXMLProducer("PubmedArticle", "", false, trdr)
				// callback passes additional closure data to xtract createConsumers
				tblq := csmr(xmlq)
				if trdr == nil || xmlq == nil || tblq == nil {
					return
				}

				concatIndexedFiles := func() string {

					var e2xs []string

					// concatenate *.e2x intermediates
					for curr := range tblq {

						str := curr.Text

						str = strings.TrimSpace(str)
						if str == "" {
							continue
						}

						e2xs = append(e2xs, str)

						runtime.Gosched()
					}

					if len(e2xs) < 0 {
						return ""
					}

					res := strings.Join(e2xs, "\n")
					res = strings.TrimSpace(res)

					return res
				}

				idxs := concatIndexedFiles()
				irdr := strings.NewReader(idxs)
				if idxs == "" || irdr == nil {
					return
				}

				// pass (e2x) indexed product through -invert
				erdr := CreateXMLStreamer(irdr)
				colq := CreateXMLProducer("IdxDocument", "", false, erdr)
				dspq := CreateDispensers(colq)
				invq := CreateInverters(dspq)
				rslq := CreateResolver(invq)
				if erdr == nil || colq == nil || dspq == nil || invq == nil || rslq == nil {
					return
				}

				concatInvertedFiles := func() string {

					var invs []string

					// concatenate *.inv components
					for str := range rslq {

						str = strings.TrimSpace(str)
						if str == "" {
							continue
						}

						invs = append(invs, str)

						runtime.Gosched()
					}

					if len(invs) < 0 {
						return ""
					}

					res := strings.Join(invs, "\n  ")
					res = strings.TrimSpace(res)

					return "<InvDocumentSet>\n  " + res + "\n</InvDocumentSet>\n"
				}

				ivss := concatInvertedFiles()
				if ivss == "" {
					return
				}

				// derive 6-digit .inv file name from first .xml file name
				file := xmls[0]
				file = strings.TrimSuffix(file, ".xml.gz")

				if IsAllDigits(file) {

					// pad numeric identifier to 8 characters with leading zeros
					ln := len(file)
					if ln < 8 {
						zeros := "00000000"
						file = zeros[ln:] + file
					}
				}
				if IsAllDigitsOrPeriod(file) {

					// limit trie to first 6 characters
					if len(file) > 6 {
						file = file[:6]
					}
				}

				writeIvss := func(wrtr *bufio.Writer) {

					wrtr.WriteString(ivss)
					if !strings.HasSuffix(ivss, "\n") {
						wrtr.WriteString("\n")
					}
				}

				// write a single .inv.gz inverted index for all .xml.gz files in folder
				stringsToGzFile(base, path, file+".inv.gz", writeIvss)

				// progress monitor prints dot every 1000 records
				countSuccess()

				out <- file + ".inv"
			}
		}

		var wg sync.WaitGroup

		// launch multiple archive indexer inverter goroutines
		for i := 0; i < NumServe(); i++ {
			wg.Add(1)
			go indexerInverter(i, &wg, archiveBase, csmr, inp, out)
		}

		// launch separate anonymous goroutine to wait until all indexer inverters are done
		go func() {
			wg.Wait()
			close(out)
		}()

		return out
	}

	vrfq := visitArchiveFolders(archiveBase)
	iniq := e2indexAndInvert(archiveBase, csmr, vrfq)

	if vrfq == nil || iniq == nil {
		return nil
	}

	return iniq
}

// JOIN CACHED INVERTED FILES FROM LOCAL ARCHIVE FOLDERS

// JoinInverts combines and sorts local archive .inv.gz inverted files
func JoinInverts(archiveBase, invertBase string) <-chan string {

	if archiveBase == "" || invertBase == "" {

		// if not passed as an argument, obtain archive base path from environment variable
		base := os.Getenv("EDIRECT_PUBMED_MASTER")
		if base != "" {
			if !strings.HasSuffix(base, "/") {
				base += "/"
			}

			archiveBase = base + "Archive"
			invertBase = base + "Invert"
		}

		if archiveBase == "" {
			fmt.Fprintf(os.Stderr, "\nERROR: EDIRECT_PUBMED_MASTER environment variable is not set\n\n")
			os.Exit(1)
		}
	}

	// check to make sure local archive is mounted
	_, err := os.Stat(archiveBase)
	if err != nil && os.IsNotExist(err) {
		fmt.Fprintf(os.Stderr, "\nERROR: Local archive and search index is not mounted\n\n")
		os.Exit(1)
	}

	outp := make(chan string, ChanDepth())
	if outp == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create joiner channel\n")
		os.Exit(1)
	}

	writeGzFile := func(base, path, name, sfx string, inp <-chan string) <-chan string {

		if inp == nil {
			return nil
		}

		out := make(chan string, ChanDepth())
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create gz file writer channel\n")
			os.Exit(1)
		}

		saveOneGzFile := func(base, path, name, sfx string, inp <-chan string, out chan<- string) {

			// close channel when all records have been processed
			defer close(out)

			fpath := filepath.Join(base, path, name+sfx)

			// overwrites and truncates existing file
			fl, err := os.Create(fpath)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			// for large files, use parallel pgzip
			zpr, err := pgzip.NewWriterLevel(fl, pgzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create compressor\n")
				os.Exit(1)
			}

			wrtr := bufio.NewWriter(zpr)

			last := ""
			for str := range inp {
				if str == "" {
					continue
				}
				wrtr.WriteString(str)
				last = str
			}
			if !strings.HasSuffix(last, "\n") {
				wrtr.WriteString("\n")
			}

			err = wrtr.Flush()
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			err = zpr.Close()
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			// fl.Sync()

			err = fl.Close()
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			out <- name
		}

		// launch single combiner goroutine
		go saveOneGzFile(base, path, name, sfx, inp, out)

		return out
	}

	rejoinInverts := func(inp <-chan string) <-chan string {

		if inp == nil {
			return nil
		}

		out := make(chan string, ChanDepth())
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create invert splitter channel\n")
			os.Exit(1)
		}

		var jlock sync.Mutex

		master := make(map[string][]string)

		// add single posting
		addPost := func(key, str string) {

			// protect map with mutex
			jlock.Lock()

			lst, ok := master[key]
			if !ok {
				lst = make([]string, 0, 1)
			}
			lst = append(lst, str)
			master[key] = lst

			jlock.Unlock()
		}

		invSplitter := func(wg *sync.WaitGroup, inp <-chan string, out chan<- string) {

			defer wg.Done()

			currKey := ""

			addInv := func(tag, attr, content string) {

				if tag == "InvKey" {
					currKey = content
				} else if attr != "" {
					str := "<" + tag + " " + attr + ">" + content + "</" + tag + ">"
					addPost(currKey, str)
				} else {
					str := "<" + tag + ">" + content + "</" + tag + ">"
					addPost(currKey, str)
				}
			}

			for txt := range inp {

				StreamValues(txt[:], "InvDocument", addInv)

				runtime.Gosched()
			}
		}

		var wg sync.WaitGroup

		// launch multiple splitter goroutines
		for i := 0; i < NumServe(); i++ {
			wg.Add(1)
			go invSplitter(&wg, inp, out)
		}

		// launch separate anonymous goroutine to wait until all splitter are done
		go func() {
			wg.Wait()

			var keys []string
			for ky := range master {
				keys = append(keys, ky)
			}
			sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

			for _, key := range keys {

				lst, ok := master[key]
				if !ok {
					continue
				}

				itms := strings.Join(lst, "\n")
				if itms == "" {
					continue
				}

				out <- "<InvDocument>\n<InvKey>" + key + "</InvKey>\n<InvIDs>\n" + itms + "\n</InvIDs>\n</InvDocument>\n"
			}

			close(out)
		}()

		return out
	}

	combineInverts := func(archiveBase, invertBase string, out chan<- string) {

		// close channel when all records have been sent
		defer close(out)

		num := 0

		getSubFolders := func(path string, mids []string, fr, to string) {

			// increment output file number
			num++

			sfx := fmt.Sprintf("%03d", num)

			target := filepath.Join(invertBase, "pubmed"+sfx+".inv.gz")

			// previous joined index files in range of updated records were removed by pm-stash
			_, err := os.Stat(target)
			if err == nil {
				// if joined index file exists, no need to recreate
				return
			}

			var folders []string

			// filter second-level directories by range
			for _, val := range mids {
				if val >= fr && val <= to {
					folders = append(folders, val)
				}
			}

			var filenames []string

			// collect precomputed inverted index files
			for _, fld := range folders {
				mid := filepath.Join(path, fld)
				dirs, _, _ := examineFolder(archiveBase, mid)
				for _, lf := range dirs {
					sub := filepath.Join(mid, lf)
					_, _, invs := examineFolder(archiveBase, sub)
					for _, file := range invs {
						fpath := filepath.Join(archiveBase, sub, file)
						filenames = append(filenames, fpath)
					}
				}
			}

			if len(filenames) < 1 {
				return
			}

			// for large files, use parallel pgzip
			invq := readGzFiles(filenames, true, 10)
			splq := rejoinInverts(invq)
			wivq := writeGzFile(invertBase, "", "pubmed"+sfx, ".inv.gz", splq)

			for fname := range wivq {
				out <- fname + ".inv"
			}

			// force garbage collection
			runtime.GC()
			debug.FreeOSMemory()

			runtime.Gosched()
		}

		// get top-level directories
		dirs, _, _ := examineFolder(archiveBase, "")

		// collect inverted index files in groups derived from 250000 pubmed .xml files
		for _, dr := range dirs {
			mids, _, _ := examineFolder(archiveBase, dr)
			getSubFolders(dr, mids, "00", "24")
			getSubFolders(dr, mids, "25", "49")
			getSubFolders(dr, mids, "50", "74")
			getSubFolders(dr, mids, "75", "99")
		}
	}

	// launch single combiner goroutine
	go combineInverts(archiveBase, invertBase, outp)

	return outp
}

// deleteCacheFiles visits leaf folders in local archive, conditionally removing
// cached .xml.gz archive files and .inv.gz indexed and inverted files
func deleteCacheFiles(archiveBase string, removeInv, removeXML bool) <-chan string {

	if archiveBase == "" {

		// if not passed as an argument, obtain archive base path from environment variable
		base := os.Getenv("EDIRECT_PUBMED_MASTER")
		if base != "" {
			if !strings.HasSuffix(base, "/") {
				base += "/"
			}

			archiveBase = base + "Archive"
		}

		if archiveBase == "" {
			fmt.Fprintf(os.Stderr, "\nERROR: EDIRECT_PUBMED_MASTER environment variable is not set\n\n")
			os.Exit(1)
		}
	}

	// check to make sure local archive is mounted
	_, err := os.Stat(archiveBase)
	if err != nil && os.IsNotExist(err) {
		fmt.Fprintf(os.Stderr, "\nERROR: Local archive and search index is not mounted\n\n")
		os.Exit(1)
	}

	inspectArchiveFolders := func(base string) <-chan string {

		out := make(chan string, ChanDepth())
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create archive inspector channel\n")
			os.Exit(1)
		}

		// recursive definition
		var inspectSubFolders func(base, path string, out chan<- string)

		// inspectSubFolders recurses to leaf directories
		inspectSubFolders = func(base, path string, out chan<- string) {

			dirs, xmls, invs := examineFolder(base, path)

			// recursively explore subdirectories
			if dirs != nil {
				for _, dr := range dirs {
					sub := filepath.Join(path, dr)
					inspectSubFolders(base, sub, out)
				}
				return
			}

			// if clearing indices, report (single) inv file
			if removeInv && invs != nil {
				for _, iv := range invs {
					// .inv.gz inverted index file
					out <- filepath.Join(base, path, iv)
				}
			}

			// if clearing archive for repopulation, report all xml files
			if removeXML && xmls != nil {
				for _, xm := range xmls {
					// .xml.gz archive file
					out <- filepath.Join(base, path, xm)
				}
			}
		}

		inspectArchiveSubset := func(wg *sync.WaitGroup, base, path string, out chan<- string) {

			defer wg.Done()

			inspectSubFolders(base, path, out)
		}

		var wg sync.WaitGroup

		// launch multiple archive inspector goroutines
		dirs, _, _ := examineFolder(archiveBase, "")
		for _, dr := range dirs {
			wg.Add(1)
			go inspectArchiveSubset(&wg, archiveBase, dr, out)
		}

		// launch separate anonymous goroutine to wait until all archive inspector are done
		go func() {
			wg.Wait()
			close(out)
		}()

		return out
	}

	deleteInverts := func(base string, inp <-chan string) <-chan string {

		if inp == nil {
			return nil
		}

		out := make(chan string, ChanDepth())
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create deleter channel\n")
			os.Exit(1)
		}

		indexDeleter := func(wg *sync.WaitGroup, inp <-chan string, out chan<- string) {

			defer wg.Done()

			for fpath := range inp {

				// delete indicated file
				os.Remove(fpath)

				// report for file count
				out <- fpath
			}
		}

		var wg sync.WaitGroup

		// launch multiple archive deleter goroutines
		for i := 0; i < NumServe(); i++ {
			wg.Add(1)
			go indexDeleter(&wg, inp, out)
		}

		// launch separate anonymous goroutine to wait until all deleters are done
		go func() {
			wg.Wait()
			close(out)
		}()

		return out
	}

	irfq := inspectArchiveFolders(archiveBase)
	dinq := deleteInverts(archiveBase, irfq)

	if irfq == nil || dinq == nil {
		return nil
	}

	return dinq
}

// DeleteInvertCache clears cached .inv.gz inverted index files (use if Entrez indexing logic has changed)
func DeleteInvertCache(archiveBase string) <-chan string {

	return deleteCacheFiles(archiveBase, true, false)
}

// DeleteArchiveCache clears all archived .xml.gz and .inv.gz files (use if cleanup/normalize logic has changed)
func DeleteArchiveCache(archiveBase string) <-chan string {

	return deleteCacheFiles(archiveBase, true, true)
}
