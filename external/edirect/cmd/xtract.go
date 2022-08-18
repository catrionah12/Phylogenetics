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
// File Name:  xtract.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package main

import (
	"bufio"
	"eutils"
	"fmt"
	"io"
	"math"
	"os"
	"runtime"
	"runtime/debug"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
	"unicode"
)

// GLOBAL VARIABLES

var (
	doStem bool
	deStop bool
)

// TYPED CONSTANTS

// LevelType is the integer type for exploration arguments
type LevelType int

// LevelType keys for exploration arguments
const (
	_ LevelType = iota
	UNIT
	SUBSET
	SECTION
	BLOCK
	BRANCH
	GROUP
	DIVISION
	PATH
	PATTERN
)

// IndentType is the integer type for XML formatting
type IndentType int

// IndentType keys for XML formatting
const (
	SINGULARITY IndentType = iota
	COMPACT
	FLUSH
	INDENT
	SUBTREE
	WRAPPED
)

// OpType is the integer type for operations
type OpType int

// OpType keys for operations
const (
	UNSET OpType = iota
	ELEMENT
	FIRST
	LAST
	BACKWARD
	ENCODE
	DECODE
	UPPER
	LOWER
	CHAIN
	TITLE
	BASIC
	PLAIN
	SIMPLE
	AUTHOR
	PROSE
	ORDER
	YEAR
	MONTH
	PAGE
	AUTH
	PROP
	TRIM
	WCT
	DOI
	TRANSLATE
	REPLACE
	TERMS
	WORDS
	PAIRS
	PAIRX
	REVERSE
	LETTERS
	CLAUSES
	INDICES
	ARTICLE
	MESHCODE
	MATRIX
	HISTOGRAM
	ACCENTED
	SCAN
	PFX
	SFX
	SEP
	TAB
	RET
	LBL
	CLR
	PFC
	DEQ
	PLG
	ELG
	FWD
	AWD
	WRP
	ENC
	PKG
	RST
	DEF
	REG
	EXP
	COLOR
	POSITION
	SELECT
	IF
	UNLESS
	MATCH
	AVOID
	AND
	OR
	EQUALS
	CONTAINS
	INCLUDES
	ISWITHIN
	STARTSWITH
	ENDSWITH
	ISNOT
	ISBEFORE
	ISAFTER
	MATCHES
	RESEMBLES
	ISEQUALTO
	DIFFERSFROM
	GT
	GE
	LT
	LE
	EQ
	NE
	NUM
	LEN
	SUM
	ACC
	MIN
	MAX
	INC
	DEC
	SUB
	AVG
	DEV
	MED
	MUL
	DIV
	MOD
	BIN
	OCT
	HEX
	BIT
	RAW
	ZEROBASED
	ONEBASED
	UCSCBASED
	REVCOMP
	NUCLEIC
	FASTA
	NCBI2NA
	NCBI4NA
	MOLWT
	HGVS
	ELSE
	VARIABLE
	ACCUMULATOR
	VALUE
	QUESTION
	TILDE
	STAR
	DOT
	PRCNT
	DOLLAR
	ATSIGN
	COUNT
	LENGTH
	DEPTH
	INDEX
	UNRECOGNIZED
)

// ArgumentType is the integer type for argument classification
type ArgumentType int

// ArgumentType keys for argument classification
const (
	_ ArgumentType = iota
	EXPLORATION
	CONDITIONAL
	EXTRACTION
	CUSTOMIZATION
)

// RangeType is the integer type for element range choices
type RangeType int

// RangeType keys for element range choices
const (
	NORANGE RangeType = iota
	STRINGRANGE
	VARIABLERANGE
	INTEGERRANGE
)

// SeqEndType is used for -ucsc-based decisions
type SeqEndType int

// SeqEndType keys for -ucsc-based decisions
const (
	_ SeqEndType = iota
	ISSTART
	ISSTOP
	ISPOS
)

// SequenceType is used to record XML tag and position for -ucsc-based
type SequenceType struct {
	Based int
	Which SeqEndType
}

// MUTEXES

var hlock sync.Mutex

var slock sync.RWMutex

// ARGUMENT MAPS

var argTypeIs = map[string]ArgumentType{
	"-unit":         EXPLORATION,
	"-Unit":         EXPLORATION,
	"-subset":       EXPLORATION,
	"-Subset":       EXPLORATION,
	"-section":      EXPLORATION,
	"-Section":      EXPLORATION,
	"-block":        EXPLORATION,
	"-Block":        EXPLORATION,
	"-branch":       EXPLORATION,
	"-Branch":       EXPLORATION,
	"-group":        EXPLORATION,
	"-Group":        EXPLORATION,
	"-division":     EXPLORATION,
	"-Division":     EXPLORATION,
	"-path":         EXPLORATION,
	"-Path":         EXPLORATION,
	"-pattern":      EXPLORATION,
	"-Pattern":      EXPLORATION,
	"-position":     CONDITIONAL,
	"-select":       CONDITIONAL,
	"-if":           CONDITIONAL,
	"-unless":       CONDITIONAL,
	"-match":        CONDITIONAL,
	"-avoid":        CONDITIONAL,
	"-and":          CONDITIONAL,
	"-or":           CONDITIONAL,
	"-equals":       CONDITIONAL,
	"-contains":     CONDITIONAL,
	"-includes":     CONDITIONAL,
	"-is-within":    CONDITIONAL,
	"-starts-with":  CONDITIONAL,
	"-ends-with":    CONDITIONAL,
	"-is-not":       CONDITIONAL,
	"-is-before":    CONDITIONAL,
	"-is-after":     CONDITIONAL,
	"-matches":      CONDITIONAL,
	"-resembles":    CONDITIONAL,
	"-is-equal-to":  CONDITIONAL,
	"-differs-from": CONDITIONAL,
	"-gt":           CONDITIONAL,
	"-ge":           CONDITIONAL,
	"-lt":           CONDITIONAL,
	"-le":           CONDITIONAL,
	"-eq":           CONDITIONAL,
	"-ne":           CONDITIONAL,
	"-element":      EXTRACTION,
	"-first":        EXTRACTION,
	"-last":         EXTRACTION,
	"-backward":     EXTRACTION,
	"-encode":       EXTRACTION,
	"-decode":       EXTRACTION,
	"-decode64":     EXTRACTION,
	"-upper":        EXTRACTION,
	"-lower":        EXTRACTION,
	"-chain":        EXTRACTION,
	"-title":        EXTRACTION,
	"-basic":        EXTRACTION,
	"-plain":        EXTRACTION,
	"-simple":       EXTRACTION,
	"-author":       EXTRACTION,
	"-prose":        EXTRACTION,
	"-order":        EXTRACTION,
	"-year":         EXTRACTION,
	"-month":        EXTRACTION,
	"-page":         EXTRACTION,
	"-auth":         EXTRACTION,
	"-prop":         EXTRACTION,
	"-trim":         EXTRACTION,
	"-wct":          EXTRACTION,
	"-doi":          EXTRACTION,
	"-translate":    EXTRACTION,
	"-replace":      EXTRACTION,
	"-terms":        EXTRACTION,
	"-words":        EXTRACTION,
	"-pairs":        EXTRACTION,
	"-pairx":        EXTRACTION,
	"-reverse":      EXTRACTION,
	"-letters":      EXTRACTION,
	"-clauses":      EXTRACTION,
	"-indices":      EXTRACTION,
	"-article":      EXTRACTION,
	"-meshcode":     EXTRACTION,
	"-matrix":       EXTRACTION,
	"-histogram":    EXTRACTION,
	"-accented":     EXTRACTION,
	"-scan":         EXTRACTION,
	"-num":          EXTRACTION,
	"-len":          EXTRACTION,
	"-sum":          EXTRACTION,
	"-acc":          EXTRACTION,
	"-min":          EXTRACTION,
	"-max":          EXTRACTION,
	"-inc":          EXTRACTION,
	"-dec":          EXTRACTION,
	"-sub":          EXTRACTION,
	"-avg":          EXTRACTION,
	"-dev":          EXTRACTION,
	"-med":          EXTRACTION,
	"-mul":          EXTRACTION,
	"-div":          EXTRACTION,
	"-mod":          EXTRACTION,
	"-bin":          EXTRACTION,
	"-oct":          EXTRACTION,
	"-hex":          EXTRACTION,
	"-bit":          EXTRACTION,
	"-raw":          EXTRACTION,
	"-0-based":      EXTRACTION,
	"-zero-based":   EXTRACTION,
	"-1-based":      EXTRACTION,
	"-one-based":    EXTRACTION,
	"-ucsc":         EXTRACTION,
	"-ucsc-based":   EXTRACTION,
	"-ucsc-coords":  EXTRACTION,
	"-bed-based":    EXTRACTION,
	"-bed-coords":   EXTRACTION,
	"-revcomp":      EXTRACTION,
	"-nucleic":      EXTRACTION,
	"-fasta":        EXTRACTION,
	"-ncbi2na":      EXTRACTION,
	"-ncbi4na":      EXTRACTION,
	"-molwt":        EXTRACTION,
	"-hgvs":         EXTRACTION,
	"-else":         EXTRACTION,
	"-pfx":          CUSTOMIZATION,
	"-sfx":          CUSTOMIZATION,
	"-sep":          CUSTOMIZATION,
	"-tab":          CUSTOMIZATION,
	"-ret":          CUSTOMIZATION,
	"-lbl":          CUSTOMIZATION,
	"-clr":          CUSTOMIZATION,
	"-pfc":          CUSTOMIZATION,
	"-deq":          CUSTOMIZATION,
	"-plg":          CUSTOMIZATION,
	"-elg":          CUSTOMIZATION,
	"-fwd":          CUSTOMIZATION,
	"-awd":          CUSTOMIZATION,
	"-wrp":          CUSTOMIZATION,
	"-enc":          CUSTOMIZATION,
	"-pkg":          CUSTOMIZATION,
	"-rst":          CUSTOMIZATION,
	"-def":          CUSTOMIZATION,
	"-reg":          CUSTOMIZATION,
	"-exp":          CUSTOMIZATION,
	"-color":        CUSTOMIZATION,
}

var opTypeIs = map[string]OpType{
	"-element":      ELEMENT,
	"-first":        FIRST,
	"-last":         LAST,
	"-backward":     BACKWARD,
	"-encode":       ENCODE,
	"-decode":       DECODE,
	"-decode64":     DECODE,
	"-upper":        UPPER,
	"-lower":        LOWER,
	"-chain":        CHAIN,
	"-title":        TITLE,
	"-basic":        BASIC,
	"-plain":        PLAIN,
	"-simple":       SIMPLE,
	"-author":       AUTHOR,
	"-prose":        PROSE,
	"-order":        ORDER,
	"-year":         YEAR,
	"-month":        MONTH,
	"-page":         PAGE,
	"-auth":         AUTH,
	"-prop":         PROP,
	"-trim":         TRIM,
	"-wct":          WCT,
	"-doi":          DOI,
	"-translate":    TRANSLATE,
	"-replace":      REPLACE,
	"-terms":        TERMS,
	"-words":        WORDS,
	"-pairs":        PAIRS,
	"-pairx":        PAIRX,
	"-reverse":      REVERSE,
	"-letters":      LETTERS,
	"-clauses":      CLAUSES,
	"-indices":      INDICES,
	"-article":      ARTICLE,
	"-meshcode":     MESHCODE,
	"-matrix":       MATRIX,
	"-histogram":    HISTOGRAM,
	"-accented":     ACCENTED,
	"-scan":         SCAN,
	"-pfx":          PFX,
	"-sfx":          SFX,
	"-sep":          SEP,
	"-tab":          TAB,
	"-ret":          RET,
	"-lbl":          LBL,
	"-clr":          CLR,
	"-pfc":          PFC,
	"-deq":          DEQ,
	"-plg":          PLG,
	"-elg":          ELG,
	"-fwd":          FWD,
	"-awd":          AWD,
	"-wrp":          WRP,
	"-enc":          ENC,
	"-pkg":          PKG,
	"-rst":          RST,
	"-def":          DEF,
	"-reg":          REG,
	"-exp":          EXP,
	"-color":        COLOR,
	"-position":     POSITION,
	"-select":       SELECT,
	"-if":           IF,
	"-unless":       UNLESS,
	"-match":        MATCH,
	"-avoid":        AVOID,
	"-and":          AND,
	"-or":           OR,
	"-equals":       EQUALS,
	"-contains":     CONTAINS,
	"-includes":     INCLUDES,
	"-is-within":    ISWITHIN,
	"-starts-with":  STARTSWITH,
	"-ends-with":    ENDSWITH,
	"-is-not":       ISNOT,
	"-is-before":    ISBEFORE,
	"-is-after":     ISAFTER,
	"-matches":      MATCHES,
	"-resembles":    RESEMBLES,
	"-is-equal-to":  ISEQUALTO,
	"-differs-from": DIFFERSFROM,
	"-gt":           GT,
	"-ge":           GE,
	"-lt":           LT,
	"-le":           LE,
	"-eq":           EQ,
	"-ne":           NE,
	"-num":          NUM,
	"-len":          LEN,
	"-sum":          SUM,
	"-acc":          ACC,
	"-min":          MIN,
	"-max":          MAX,
	"-inc":          INC,
	"-dec":          DEC,
	"-sub":          SUB,
	"-avg":          AVG,
	"-dev":          DEV,
	"-med":          MED,
	"-mul":          MUL,
	"-div":          DIV,
	"-mod":          MOD,
	"-bin":          BIN,
	"-oct":          OCT,
	"-hex":          HEX,
	"-bit":          BIT,
	"-raw":          RAW,
	"-0-based":      ZEROBASED,
	"-zero-based":   ZEROBASED,
	"-1-based":      ONEBASED,
	"-one-based":    ONEBASED,
	"-ucsc":         UCSCBASED,
	"-ucsc-based":   UCSCBASED,
	"-ucsc-coords":  UCSCBASED,
	"-bed-based":    UCSCBASED,
	"-bed-coords":   UCSCBASED,
	"-revcomp":      REVCOMP,
	"-nucleic":      NUCLEIC,
	"-fasta":        FASTA,
	"-ncbi2na":      NCBI2NA,
	"-ncbi4na":      NCBI4NA,
	"-molwt":        MOLWT,
	"-hgvs":         HGVS,
	"-else":         ELSE,
}

var sequenceTypeIs = map[string]SequenceType{
	"INSDSeq:INSDInterval_from":       {1, ISSTART},
	"INSDSeq:INSDInterval_to":         {1, ISSTOP},
	"DocumentSummary:ChrStart":        {0, ISSTART},
	"DocumentSummary:ChrStop":         {0, ISSTOP},
	"DocumentSummary:Chr_start":       {1, ISSTART},
	"DocumentSummary:Chr_end":         {1, ISSTOP},
	"DocumentSummary:Chr_inner_start": {1, ISSTART},
	"DocumentSummary:Chr_inner_end":   {1, ISSTOP},
	"DocumentSummary:Chr_outer_start": {1, ISSTART},
	"DocumentSummary:Chr_outer_end":   {1, ISSTOP},
	"DocumentSummary:start":           {1, ISSTART},
	"DocumentSummary:stop":            {1, ISSTOP},
	"DocumentSummary:display_start":   {1, ISSTART},
	"DocumentSummary:display_stop":    {1, ISSTOP},
	"Entrezgene:Seq-interval_from":    {0, ISSTART},
	"Entrezgene:Seq-interval_to":      {0, ISSTOP},
	"GenomicInfoType:ChrStart":        {0, ISSTART},
	"GenomicInfoType:ChrStop":         {0, ISSTOP},
	"RS:position":                     {0, ISPOS},
	"RS:@asnFrom":                     {0, ISSTART},
	"RS:@asnTo":                       {0, ISSTOP},
	"RS:@end":                         {0, ISSTOP},
	"RS:@leftContigNeighborPos":       {0, ISSTART},
	"RS:@physMapInt":                  {0, ISPOS},
	"RS:@protLoc":                     {0, ISPOS},
	"RS:@rightContigNeighborPos":      {0, ISSTOP},
	"RS:@start":                       {0, ISSTART},
	"RS:@structLoc":                   {0, ISPOS},
}

var monthTable = map[string]int{
	"jan":       1,
	"january":   1,
	"feb":       2,
	"february":  2,
	"mar":       3,
	"march":     3,
	"apr":       4,
	"april":     4,
	"may":       5,
	"jun":       6,
	"june":      6,
	"jul":       7,
	"julu":      7,
	"aug":       8,
	"august":    8,
	"sep":       9,
	"september": 9,
	"oct":       10,
	"october":   10,
	"nov":       11,
	"november":  11,
	"dec":       12,
	"december":  12,
}

var propertyTable = map[string]string{
	"AssociatedDataset":           "Associated Dataset",
	"AssociatedPublication":       "Associated Publication",
	"CommentIn":                   "Comment In",
	"CommentOn":                   "Comment On",
	"ErratumFor":                  "Erratum For",
	"ErratumIn":                   "Erratum In",
	"ExpressionOfConcernFor":      "Expression Of Concern For",
	"ExpressionOfConcernIn":       "Expression Of Concern In",
	"OriginalReportIn":            "Original Report In",
	"ReprintIn":                   "Reprint In",
	"ReprintOf":                   "Reprint Of",
	"RepublishedFrom":             "Republished From",
	"RepublishedIn":               "Republished In",
	"RetractedandRepublishedFrom": "Retracted And Republished From",
	"RetractedandRepublishedIn":   "Retracted And Republished In",
	"RetractionIn":                "Retraction In",
	"RetractionOf":                "Retraction Of",
	"SummaryForPatientsIn":        "Summary For Patients In",
	"UpdateIn":                    "Update In",
	"UpdateOf":                    "Update Of",
	"aheadofprint":                "Ahead Of Print",
	"epublish":                    "Electronically Published",
	"ppublish":                    "Published In Print",
}

// DATA OBJECTS

// Step contains parameters for executing a single command step
type Step struct {
	Type   OpType
	Value  string
	Parent string
	Match  string
	Attrib string
	TypL   RangeType
	StrL   string
	IntL   int
	TypR   RangeType
	StrR   string
	IntR   int
	Norm   bool
	Wild   bool
	Unesc  bool
}

// Operation breaks commands into sequential steps
type Operation struct {
	Type   OpType
	Value  string
	Stages []*Step
}

// Block contains nested instructions for executing commands
type Block struct {
	Visit      string
	Parent     string
	Match      string
	Path       []string
	Working    []string
	Parsed     []string
	Position   string
	Foreword   string
	Afterword  string
	Conditions []*Operation
	Commands   []*Operation
	Failure    []*Operation
	Subtasks   []*Block
}

// Limiter is used for collecting specific nodes (e.g., first and last)
type Limiter struct {
	Obj *eutils.XMLNode
	Idx int
	Lvl int
}

// UTILITIES

func hasSpaceOrHyphen(str string) bool {

	for _, ch := range str {
		if ch == ' ' || ch == '-' {
			return true
		}
	}

	return false
}

func isAllCapsOrDigits(str string) bool {

	for _, ch := range str {
		if !unicode.IsUpper(ch) && !unicode.IsDigit(ch) {
			return false
		}
	}

	return true
}

// hasCommaOrSemicolon reports on comma, semicolon, or hyphen
func hasCommaOrSemicolon(str string) bool {

	for _, ch := range str {
		if ch == ',' || ch == ';' || ch == '-' {
			return true
		}
	}

	return false
}

// removeCommaOrSemicolon replaces comma or semicolon with space
func removeCommaOrSemicolon(str string) string {

	str = strings.ToLower(str)

	if hasCommaOrSemicolon(str) {
		str = strings.Replace(str, ",", " ", -1)
		str = strings.Replace(str, ";", " ", -1)
		str = eutils.CompressRunsOfSpaces(str)
	}
	str = strings.TrimSpace(str)
	str = strings.TrimRight(str, ".?:")

	return str
}

// sortStringByWords sorts the individual words in a string
func sortStringByWords(str string) string {

	str = removeCommaOrSemicolon(str)

	// check for multiple words
	if hasSpaceOrHyphen(str) {
		flds := strings.Fields(str)
		sort.Slice(flds, func(i, j int) bool { return flds[i] < flds[j] })
		str = strings.Join(flds, " ")
		str = strings.Replace(str, "-", " ", -1)
		str = eutils.CompressRunsOfSpaces(str)
		str = strings.TrimRight(str, ".?:")
	}

	return str
}

/*
func parseMarkup(str, cmd string) int {

	switch str {
	case "fuse", "fused":
		return eutils.FUSE
	case "space", "spaces":
		return eutils.SPACE
	case "period", "periods":
		return eutils.PERIOD
	case "concise":
		return eutils.CONCISE
	case "bracket", "brackets":
		return eutils.BRACKETS
	case "markdown":
		return eutils.MARKDOWN
	case "slash":
		return eutils.SLASH
	case "tag", "tags":
		return eutils.TAGS
	case "terse":
		return eutils.TERSE
	default:
		if str != "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized %s value '%s'\n", cmd, str)
			os.Exit(1)
		}
	}
	return eutils.NOMARKUP
}
*/

// INSDSEQ EXTRACTION COMMAND GENERATOR

// e.g., xtract -insd complete mat_peptide "%peptide" product peptide

// processINSD generates extraction commands for GenBank/RefSeq records in INSDSet format
func processINSD(args []string, isPipe, addDash, doIndex bool) []string {

	// legal GenBank / GenPept / RefSeq features

	features := []string{
		"-10_signal",
		"-35_signal",
		"3'clip",
		"3'UTR",
		"5'clip",
		"5'UTR",
		"allele",
		"assembly_gap",
		"attenuator",
		"Bond",
		"C_region",
		"CAAT_signal",
		"CDS",
		"centromere",
		"conflict",
		"D_segment",
		"D-loop",
		"enhancer",
		"exon",
		"gap",
		"GC_signal",
		"gene",
		"iDNA",
		"intron",
		"J_segment",
		"LTR",
		"mat_peptide",
		"misc_binding",
		"misc_difference",
		"misc_feature",
		"misc_recomb",
		"misc_RNA",
		"misc_signal",
		"misc_structure",
		"mobile_element",
		"modified_base",
		"mRNA",
		"mutation",
		"N_region",
		"ncRNA",
		"old_sequence",
		"operon",
		"oriT",
		"polyA_signal",
		"polyA_site",
		"precursor_RNA",
		"prim_transcript",
		"primer_bind",
		"promoter",
		"propeptide",
		"protein_bind",
		"Protein",
		"RBS",
		"Region",
		"regulatory",
		"rep_origin",
		"repeat_region",
		"repeat_unit",
		"rRNA",
		"S_region",
		"satellite",
		"scRNA",
		"sig_peptide",
		"Site",
		"snoRNA",
		"snRNA",
		"source",
		"stem_loop",
		"STS",
		"TATA_signal",
		"telomere",
		"terminator",
		"tmRNA",
		"transit_peptide",
		"tRNA",
		"unsure",
		"V_region",
		"V_segment",
		"variation",
	}

	// legal GenBank / GenPept / RefSeq qualifiers

	qualifiers := []string{
		"allele",
		"altitude",
		"anticodon",
		"artificial_location",
		"bio_material",
		"bond_type",
		"bound_moiety",
		"breed",
		"calculated_mol_wt",
		"cell_line",
		"cell_type",
		"chloroplast",
		"chromoplast",
		"chromosome",
		"circular_RNA",
		"citation",
		"clone_lib",
		"clone",
		"coded_by",
		"codon_start",
		"codon",
		"collected_by",
		"collection_date",
		"compare",
		"cons_splice",
		"country",
		"cultivar",
		"culture_collection",
		"cyanelle",
		"db_xref",
		"derived_from",
		"dev_stage",
		"direction",
		"EC_number",
		"ecotype",
		"encodes",
		"endogenous_virus",
		"environmental_sample",
		"estimated_length",
		"evidence",
		"exception",
		"experiment",
		"focus",
		"frequency",
		"function",
		"gap_type",
		"gdb_xref",
		"gene_synonym",
		"gene",
		"germline",
		"haplogroup",
		"haplotype",
		"host",
		"identified_by",
		"inference",
		"insertion_seq",
		"isolate",
		"isolation_source",
		"kinetoplast",
		"lab_host",
		"label",
		"lat_lon",
		"linkage_evidence",
		"locus_tag",
		"macronuclear",
		"map",
		"mating_type",
		"metagenome_source",
		"metagenomic",
		"mitochondrion",
		"mobile_element_type",
		"mobile_element",
		"mod_base",
		"mol_type",
		"name",
		"nat_host",
		"ncRNA_class",
		"non_functional",
		"note",
		"number",
		"old_locus_tag",
		"operon",
		"organelle",
		"organism",
		"partial",
		"PCR_conditions",
		"PCR_primers",
		"peptide",
		"phenotype",
		"plasmid",
		"pop_variant",
		"product",
		"protein_id",
		"proviral",
		"pseudo",
		"pseudogene",
		"rearranged",
		"recombination_class",
		"region_name",
		"regulatory_class",
		"replace",
		"ribosomal_slippage",
		"rpt_family",
		"rpt_type",
		"rpt_unit_range",
		"rpt_unit_seq",
		"rpt_unit",
		"satellite",
		"segment",
		"sequenced_mol",
		"serotype",
		"serovar",
		"sex",
		"site_type",
		"specific_host",
		"specimen_voucher",
		"standard_name",
		"strain",
		"structural_class",
		"sub_clone",
		"sub_species",
		"sub_strain",
		"submitter_seqid",
		"tag_peptide",
		"tissue_lib",
		"tissue_type",
		"trans_splicing",
		"transcript_id",
		"transcription",
		"transgenic",
		"transl_except",
		"transl_table",
		"translation",
		"transposon",
		"type_material",
		"UniProtKB_evidence",
		"usedin",
		"variety",
		"virion",
	}

	// legal INSDSeq XML fields

	insdtags := []string{
		"INSDAltSeqData_items",
		"INSDAltSeqData",
		"INSDAltSeqItem_first-accn",
		"INSDAltSeqItem_gap-comment",
		"INSDAltSeqItem_gap-length",
		"INSDAltSeqItem_gap-linkage",
		"INSDAltSeqItem_gap-type",
		"INSDAltSeqItem_interval",
		"INSDAltSeqItem_isgap",
		"INSDAltSeqItem_isgap@value",
		"INSDAltSeqItem_last-accn",
		"INSDAltSeqItem_value",
		"INSDAltSeqItem",
		"INSDAuthor",
		"INSDComment_paragraphs",
		"INSDComment_type",
		"INSDComment",
		"INSDCommentParagraph",
		"INSDFeature_intervals",
		"INSDFeature_key",
		"INSDFeature_location",
		"INSDFeature_operator",
		"INSDFeature_partial3",
		"INSDFeature_partial3@value",
		"INSDFeature_partial5",
		"INSDFeature_partial5@value",
		"INSDFeature_quals",
		"INSDFeature_xrefs",
		"INSDFeature",
		"INSDFeatureSet_annot-source",
		"INSDFeatureSet_features",
		"INSDFeatureSet",
		"INSDInterval_accession",
		"INSDInterval_from",
		"INSDInterval_interbp",
		"INSDInterval_interbp@value",
		"INSDInterval_iscomp",
		"INSDInterval_iscomp@value",
		"INSDInterval_point",
		"INSDInterval_to",
		"INSDInterval",
		"INSDKeyword",
		"INSDQualifier_name",
		"INSDQualifier_value",
		"INSDQualifier",
		"INSDReference_authors",
		"INSDReference_consortium",
		"INSDReference_journal",
		"INSDReference_position",
		"INSDReference_pubmed",
		"INSDReference_reference",
		"INSDReference_remark",
		"INSDReference_title",
		"INSDReference_xref",
		"INSDReference",
		"INSDSecondary-accn",
		"INSDSeq_accession-version",
		"INSDSeq_alt-seq",
		"INSDSeq_comment-set",
		"INSDSeq_comment",
		"INSDSeq_contig",
		"INSDSeq_create-date",
		"INSDSeq_create-release",
		"INSDSeq_database-reference",
		"INSDSeq_definition",
		"INSDSeq_division",
		"INSDSeq_entry-version",
		"INSDSeq_feature-set",
		"INSDSeq_feature-table",
		"INSDSeq_keywords",
		"INSDSeq_length",
		"INSDSeq_locus",
		"INSDSeq_moltype",
		"INSDSeq_organism",
		"INSDSeq_other-seqids",
		"INSDSeq_primary-accession",
		"INSDSeq_primary",
		"INSDSeq_project",
		"INSDSeq_references",
		"INSDSeq_secondary-accessions",
		"INSDSeq_segment",
		"INSDSeq_sequence",
		"INSDSeq_source-db",
		"INSDSeq_source",
		"INSDSeq_strandedness",
		"INSDSeq_struc-comments",
		"INSDSeq_taxonomy",
		"INSDSeq_topology",
		"INSDSeq_update-date",
		"INSDSeq_update-release",
		"INSDSeq_xrefs",
		"INSDSeq",
		"INSDSeqid",
		"INSDSet",
		"INSDStrucComment_items",
		"INSDStrucComment_name",
		"INSDStrucComment",
		"INSDStrucCommentItem_tag",
		"INSDStrucCommentItem_url",
		"INSDStrucCommentItem_value",
		"INSDStrucCommentItem",
		"INSDXref_dbname",
		"INSDXref_id",
		"INSDXref",
	}

	checkAgainstVocabulary := func(str, objtype string, arry []string) {

		if str == "" || arry == nil {
			return
		}

		// skip past pound, percent, or caret character at beginning of string
		if len(str) > 1 {
			switch str[0] {
			case '#', '%', '^':
				str = str[1:]
			default:
			}
		}

		for _, txt := range arry {
			if str == txt {
				return
			}
			if strings.ToUpper(str) == strings.ToUpper(txt) {
				fmt.Fprintf(os.Stderr, "\nERROR: Incorrect capitalization of '%s' %s, change to '%s'\n", str, objtype, txt)
				os.Exit(1)
			}
		}

		fmt.Fprintf(os.Stderr, "\nERROR: Item '%s' is not a legal -insd %s\n", str, objtype)
		os.Exit(1)
	}

	var acc []string

	max := len(args)
	if max < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract -insd\n")
		os.Exit(1)
	}

	// record accession and sequence

	if doIndex {
		if isPipe {
			acc = append(acc, "-head", "<IdxDocumentSet>", "-tail", "</IdxDocumentSet>")
			acc = append(acc, "-hd", "  <IdxDocument>\n", "-tl", "  </IdxDocument>")
			acc = append(acc, "-pattern", "INSDSeq", "-pfx", "    <IdxUid>", "-sfx", "</IdxUid>\n")
			acc = append(acc, "-element", "INSDSeq_accession-version", "-clr", "-rst", "-tab", "\n")
		} else {
			acc = append(acc, "-head", "\"<IdxDocumentSet>\"", "-tail", "\"</IdxDocumentSet>\"")
			acc = append(acc, "-hd", "\"  <IdxDocument>\\n\"", "-tl", "\"  </IdxDocument>\"")
			acc = append(acc, "-pattern", "INSDSeq", "-pfx", "\"    <IdxUid>\"", "-sfx", "\"</IdxUid>\\n\"")
			acc = append(acc, "-element", "INSDSeq_accession-version", "-clr", "-rst", "-tab", "\\n")
		}
	} else {
		acc = append(acc, "-pattern", "INSDSeq", "-ACCN", "INSDSeq_accession-version")
		acc = append(acc, "-LCUS", "INSDSeq_locus", "-SEQ", "INSDSeq_sequence")
	}

	if doIndex {
		if isPipe {
			acc = append(acc, "-group", "INSDSeq", "-lbl", "    <IdxSearchFields>\n")
		} else {
			acc = append(acc, "-group", "INSDSeq", "-lbl", "\"    <IdxSearchFields>\\n\"")
		}
	}

	printAccn := true

	// collect descriptors

	if strings.HasPrefix(args[0], "INSD") {

		if doIndex {
			acc = append(acc, "-clr", "-indices")
		} else {
			if isPipe {
				acc = append(acc, "-clr", "-pfx", "\\n", "-element", "&ACCN")
				acc = append(acc, "-group", "INSDSeq", "-sep", "|", "-element")
			} else {
				acc = append(acc, "-clr", "-pfx", "\"\\n\"", "-element", "\"&ACCN\"")
				acc = append(acc, "-group", "INSDSeq", "-sep", "\"|\"", "-element")
			}
			printAccn = false
		}

		for {
			if len(args) < 1 {
				return acc
			}
			str := args[0]
			if !strings.HasPrefix(args[0], "INSD") {
				break
			}
			checkAgainstVocabulary(str, "element", insdtags)
			acc = append(acc, str)
			args = args[1:]
		}

	} else if strings.HasPrefix(strings.ToUpper(args[0]), "INSD") {

		// report capitalization or vocabulary failure
		checkAgainstVocabulary(args[0], "element", insdtags)

		// program should not get to this point, but warn and exit anyway
		fmt.Fprintf(os.Stderr, "\nERROR: Item '%s' is not a legal -insd %s\n", args[0], "element")
		os.Exit(1)
	}

	processOneFeature := func(ftargs []string) {

		// skip past -insd feature clause separator

		if ftargs[0] == "-insd" {
			ftargs = ftargs[1:]
		}

		// collect qualifiers

		partial := false
		complete := false

		if ftargs[0] == "+" || ftargs[0] == "complete" {
			complete = true
			ftargs = ftargs[1:]
			max--
		} else if ftargs[0] == "-" || ftargs[0] == "partial" {
			partial = true
			ftargs = ftargs[1:]
			max--
		}

		if max < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: No feature key supplied to xtract -insd\n")
			os.Exit(1)
		}

		acc = append(acc, "-group", "INSDFeature")

		// limit to designated features

		feature := ftargs[0]

		fcmd := "-if"

		// can specify multiple features separated by plus sign (e.g., CDS+mRNA) or comma (e.g., CDS,mRNA)
		plus := strings.Split(feature, "+")
		for _, pls := range plus {
			comma := strings.Split(pls, ",")
			for _, cma := range comma {

				checkAgainstVocabulary(cma, "feature", features)
				acc = append(acc, fcmd, "INSDFeature_key", "-equals", cma)

				fcmd = "-or"
			}
		}

		if max < 2 {
			// still need at least one qualifier even on legal feature
			fmt.Fprintf(os.Stderr, "\nERROR: Feature '%s' must be followed by at least one qualifier\n", feature)
			os.Exit(1)
		}

		ftargs = ftargs[1:]

		if complete {
			acc = append(acc, "-unless", "INSDFeature_partial5", "-or", "INSDFeature_partial3")
		} else if partial {
			acc = append(acc, "-if", "INSDFeature_partial5", "-or", "INSDFeature_partial3")
		}

		if printAccn {
			if doIndex {
			} else {
				if isPipe {
					acc = append(acc, "-clr", "-pfx", "\\n", "-first", "&ACCN,&LCUS")
				} else {
					acc = append(acc, "-clr", "-pfx", "\"\\n\"", "-first", "\"&ACCN,&LCUS\"")
				}
				printAccn = false
			}
		}

		for _, str := range ftargs {

			if str == "mol_wt" {
				str = "calculated_mol_wt"
			}

			if strings.HasPrefix(str, "INSD") {

				checkAgainstVocabulary(str, "element", insdtags)
				if doIndex {
					acc = append(acc, "-block", "INSDFeature", "-clr", "-indices")
				} else {
					if isPipe {
						acc = append(acc, "-block", "INSDFeature", "-sep", "|", "-element")
					} else {
						acc = append(acc, "-block", "INSDFeature", "-sep", "\"|\"", "-element")
					}
				}
				acc = append(acc, str)
				if addDash {
					acc = append(acc, "-block", "INSDFeature", "-unless", str)
					if strings.HasSuffix(str, "@value") {
						if isPipe {
							acc = append(acc, "-lbl", "false")
						} else {
							acc = append(acc, "-lbl", "\"false\"")
						}
					} else {
						if isPipe {
							acc = append(acc, "-lbl", "\\-")
						} else {
							acc = append(acc, "-lbl", "\"\\-\"")
						}
					}
				}

			} else if strings.HasPrefix(str, "#INSD") {

				checkAgainstVocabulary(str, "element", insdtags)
				if doIndex {
					acc = append(acc, "-block", "INSDFeature", "-clr", "-indices")
				} else {
					if isPipe {
						acc = append(acc, "-block", "INSDFeature", "-sep", "|", "-element")
						acc = append(acc, str)
					} else {
						acc = append(acc, "-block", "INSDFeature", "-sep", "\"|\"", "-element")
						ql := fmt.Sprintf("\"%s\"", str)
						acc = append(acc, ql)
					}
				}

			} else if strings.HasPrefix(strings.ToUpper(str), "#INSD") {

				// report capitalization or vocabulary failure
				checkAgainstVocabulary(str, "element", insdtags)

			} else if str == "sub_sequence" {

				// special sub_sequence qualifier shows sequence under feature intervals
				acc = append(acc, "-block", "INSDFeature_intervals")

				acc = append(acc, "-subset", "INSDInterval", "-FR", "INSDInterval_from", "-TO", "INSDInterval_to")
				if isPipe {
					acc = append(acc, "-pfx", "", "-tab", "", "-nucleic", "&SEQ[&FR:&TO]")
				} else {
					acc = append(acc, "-pfx", "\"\"", "-tab", "\"\"", "-nucleic", "\"&SEQ[&FR:&TO]\"")
				}

				acc = append(acc, "-subset", "INSDFeature_intervals")
				if isPipe {
					acc = append(acc, "-deq", "\\t")
				} else {
					acc = append(acc, "-deq", "\"\\t\"")
				}

			} else if str == "feat_location" {

				// special feat_location qualifier shows feature intervals, in 1-based GenBank convention
				acc = append(acc, "-block", "INSDFeature_intervals")

				acc = append(acc, "-subset", "INSDInterval", "-FR", "INSDInterval_from", "-TO", "INSDInterval_to")
				if isPipe {
					acc = append(acc, "-pfx", "", "-tab", "..", "-element", "&FR")
					acc = append(acc, "-pfx", "", "-tab", ",", "-element", "&TO")
				} else {
					acc = append(acc, "-pfx", "\"\"", "-tab", "\"..\"", "-element", "\"&FR\"")
					acc = append(acc, "-pfx", "\"\"", "-tab", "\",\"", "-element", "\"&TO\"")
				}

				acc = append(acc, "-subset", "INSDFeature_intervals")
				if isPipe {
					acc = append(acc, "-deq", "\\t")
				} else {
					acc = append(acc, "-deq", "\"\\t\"")
				}

			} else if str == "feat_intervals" {

				// special feat_intervals qualifier shows feature intervals, decremented to 0-based
				acc = append(acc, "-block", "INSDFeature_intervals")

				acc = append(acc, "-subset", "INSDInterval")
				if isPipe {
					acc = append(acc, "-pfx", "", "-tab", "..", "-dec", "INSDInterval_from")
					acc = append(acc, "-pfx", "", "-tab", ",", "-dec", "INSDInterval_to")
				} else {
					acc = append(acc, "-pfx", "\"\"", "-tab", "\"..\"", "-dec", "\"INSDInterval_from\"")
					acc = append(acc, "-pfx", "\"\"", "-tab", "\",\"", "-dec", "\"INSDInterval_to\"")
				}

				acc = append(acc, "-subset", "INSDFeature_intervals")
				if isPipe {
					acc = append(acc, "-deq", "\\t")
				} else {
					acc = append(acc, "-deq", "\"\\t\"")
				}

			} else if str == "chloroplast" ||
				str == "chromoplast" ||
				str == "cyanelle" ||
				str == "environmental_sample" ||
				str == "focus" ||
				str == "germline" ||
				str == "kinetoplast" ||
				str == "macronuclear" ||
				str == "metagenomic" ||
				str == "mitochondrion" ||
				str == "partial" ||
				str == "proviral" ||
				str == "pseudo" ||
				str == "rearranged" ||
				str == "ribosomal_slippage" ||
				str == "trans_splicing" ||
				str == "transgenic" ||
				str == "virion" {

				acc = append(acc, "-block", "INSDQualifier")

				checkAgainstVocabulary(str, "qualifier", qualifiers)
				if doIndex {
					acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
					acc = append(acc, "-clr", "-indices", "INSDQualifier_name")
				} else {
					acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
					acc = append(acc, "-lbl", str)
				}
				if addDash {
					acc = append(acc, "-block", "INSDFeature", "-unless", "INSDQualifier_name", "-equals", str)
					if isPipe {
						acc = append(acc, "-lbl", "\\-")
					} else {
						acc = append(acc, "-lbl", "\"\\-\"")
					}
				}

			} else {

				acc = append(acc, "-block", "INSDQualifier")

				isTaxID := false
				if feature == "source" && (str == "taxon" || str == "taxid") {
					// special taxid qualifier extracts number from taxon db_xref
					isTaxID = true
					str = "db_xref"
				} else {
					checkAgainstVocabulary(str, "qualifier", qualifiers)
				}

				if len(str) > 2 && str[0] == '%' {
					acc = append(acc, "-if", "INSDQualifier_name", "-equals", str[1:])
					if doIndex {
						if isPipe {
							acc = append(acc, "-clr", "-indices", "%INSDQualifier_value")
						} else {
							acc = append(acc, "-clr", "-indices", "\"%INSDQualifier_value\"")
						}
					} else {
						if isPipe {
							acc = append(acc, "-element", "%INSDQualifier_value")
						} else {
							acc = append(acc, "-element", "\"%INSDQualifier_value\"")
						}
					}
					if addDash {
						acc = append(acc, "-block", "INSDFeature", "-unless", "INSDQualifier_name", "-equals", str[1:])
						if isPipe {
							acc = append(acc, "-lbl", "\\-")
						} else {
							acc = append(acc, "-lbl", "\"\\-\"")
						}
					}
				} else {
					if doIndex {
						acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
						acc = append(acc, "-clr", "-indices", "INSDQualifier_value")
					} else if isTaxID {
						acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
						acc = append(acc, "-and", "INSDQualifier_value", "-starts-with", "taxon:")
						acc = append(acc, "-element", "INSDQualifier_value[taxon:|]")
					} else {
						acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
						acc = append(acc, "-element", "INSDQualifier_value")
					}
					if addDash {
						if isTaxID {
							acc = append(acc, "-block", "INSDFeature", "-unless", "INSDQualifier_value", "-starts-with", "taxon:")
						} else {
							acc = append(acc, "-block", "INSDFeature", "-unless", "INSDQualifier_name", "-equals", str)
						}
						if isPipe {
							acc = append(acc, "-lbl", "\\-")
						} else {
							acc = append(acc, "-lbl", "\"\\-\"")
						}
					}
				}
			}
		}
	}

	// multiple feature clauses are separated by additional -insd arguments

	last := 0
	curr := 0
	nxt := ""

	for curr, nxt = range args {
		if nxt == "-insd" {
			if last < curr {
				processOneFeature(args[last:curr])
				last = curr
			}
		}
	}

	if last < curr {
		processOneFeature(args[last:])
	}

	if doIndex {
		if isPipe {
			acc = append(acc, "-group", "INSDSeq", "-clr", "-lbl", "    </IdxSearchFields>\n")
		} else {
			acc = append(acc, "-group", "INSDSeq", "-clr", "-lbl", "\"    </IdxSearchFields>\\n\"")
		}
	}

	return acc
}

// BIOTHINGS EXTRACTION COMMAND GENERATOR

// processBiopath generates extraction commands for BioThings resources (undocumented)
func processBiopath(args []string, isPipe bool) []string {

	// nquire -get "http://myvariant.info/v1/variant/chr6:g.26093141G>A" \
	//   -fields clinvar.rcv.conditions.identifiers \
	//   -always_list clinvar.rcv.conditions.identifiers |
	// transmute -j2x |
	// xtract -biopath opt clinvar.rcv.conditions.identifiers.omim

	var acc []string

	max := len(args)
	if max < 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract -biopath\n")
		os.Exit(1)
	}

	obj := args[0]
	args = args[1:]

	acc = append(acc, "-pattern", obj)

	paths := args[0]

	items := strings.Split(paths, ",")

	for _, path := range items {

		dirs := strings.Split(path, ".")
		max = len(dirs)
		if max < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient path arguments supplied to xtract -biopath\n")
			os.Exit(1)
		}
		if max > 7 {
			fmt.Fprintf(os.Stderr, "\nERROR: Too many nodes in argument supplied to xtract -biopath\n")
			os.Exit(1)
		}

		str := dirs[max-1]

		acc = append(acc, "-path")
		if isPipe {
			acc = append(acc, path)
			acc = append(acc, "-tab", "\\n")
			acc = append(acc, "-element", str)
		} else {
			acc = append(acc, "\""+path+"\"")
			acc = append(acc, "-tab", "\"\\n\"")
			acc = append(acc, "-element", "\""+str+"\"")
		}
	}

	return acc
}

// -select SUPPORT FUNCTIONS

func createSelectors(parent, indx string, order map[string]bool, inp <-chan eutils.XMLRecord) <-chan eutils.XMLRecord {

	if parent == "" || indx == "" || order == nil || inp == nil {
		return nil
	}

	find := eutils.ParseIndex(indx)

	out := make(chan eutils.XMLRecord, eutils.ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create selector channel\n")
		os.Exit(1)
	}

	// xmlSelector reads partitioned XML from channel and matches identifiers of records to keep
	xmlSelector := func(wg *sync.WaitGroup, inp <-chan eutils.XMLRecord, out chan<- eutils.XMLRecord) {

		// report when this selector has no more records to process
		defer wg.Done()

		// read partitioned XML from producer channel
		for ext := range inp {

			text := ext.Text

			found := false

			eutils.FindIdentifiers(text[:], parent, find,
				func(id string) {
					id = sortStringByWords(id)
					_, ok := order[id]
					if ok {
						found = true
					}
				})

			if !found {
				// identifier field not found or not in identifier list, send empty placeholder for unshuffler
				out <- eutils.XMLRecord{Index: ext.Index}
				continue
			}

			// send selected record
			out <- eutils.XMLRecord{Index: ext.Index, Text: text}
		}
	}

	var wg sync.WaitGroup

	// launch multiple selector goroutines
	for i := 0; i < eutils.NumServe(); i++ {
		wg.Add(1)
		go xmlSelector(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all selectors are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func createUnicoders(inp <-chan eutils.XMLRecord) <-chan eutils.XMLRecord {

	if inp == nil {
		return nil
	}

	out := make(chan eutils.XMLRecord, eutils.ChanDepth())
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create selector channel\n")
		os.Exit(1)
	}

	// xmlUnicoder reads partitioned XML from channel and keeps records with non-ASCII characters
	xmlUnicoder := func(wg *sync.WaitGroup, inp <-chan eutils.XMLRecord, out chan<- eutils.XMLRecord) {

		// report when this selector has no more records to process
		defer wg.Done()

		// read partitioned XML from producer channel
		for ext := range inp {

			text := ext.Text

			if !eutils.IsNotASCII(text) {
				// if only ASCII, send empty placeholder for unshuffler
				out <- eutils.XMLRecord{Index: ext.Index}
				continue
			}

			// send selected record
			out <- eutils.XMLRecord{Index: ext.Index, Text: text}
		}
	}

	var wg sync.WaitGroup

	// launch multiple unicoder goroutines
	for i := 0; i < eutils.NumServe(); i++ {
		wg.Add(1)
		go xmlUnicoder(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all unicoders are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// MAIN FUNCTION

// e.g., xtract -pattern PubmedArticle -element MedlineCitation/PMID -block Author -sep " " -element Initials,LastName

func main() {

	// skip past executable name
	args := os.Args[1:]

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// performance arguments
	chanDepth := 0
	farmSize := 0
	heapSize := 0
	numServe := 0
	goGc := 0

	// processing option arguments
	doCompress := false
	doCleanup := false
	doStrict := false
	doMixed := false
	doSelf := false
	deAccent := false
	deSymbol := false
	doASCII := false
	doStem = false
	deStop = true

	/*
		doUnicode := false
		doScript := false
		doMathML := false
	*/

	// CONCURRENCY, CLEANUP, AND DEBUGGING FLAGS

	// do these first because -defcpu and -maxcpu can be sent from wrapper before other arguments

	ncpu := runtime.NumCPU()
	if ncpu < 1 {
		ncpu = 1
	}

	// wrapper can limit maximum number of processors to use (undocumented)
	maxProcs := ncpu
	defProcs := 0

	// concurrent performance tuning parameters, can be overridden by -proc and -cons
	numProcs := 0
	serverRatio := 4

	// -flag sets -strict or -mixed cleanup flags from argument
	flgs := ""

	/*
		unicodePolicy := ""
		scriptPolicy := ""
		mathmlPolicy := ""
	*/

	// read data from file instead of stdin
	fileName := ""

	// flag for indexed input file
	turbo := false

	// debugging
	mpty := false
	idnt := false
	stts := false
	timr := false

	// profiling
	prfl := false

	// repeat the specified extraction 5 times for each -proc from 1 to nCPU
	trial := false

	inSwitch := true

	// get concurrency, cleanup, and debugging flags in any order
	for {

		inSwitch = true

		switch args[0] {
		// concurrency override arguments can be passed in by local wrapper script (undocumented)
		case "-maxcpu":
			maxProcs = eutils.GetNumericArg(args, "Maximum number of processors", 1, 1, ncpu)
			args = args[1:]
		case "-defcpu":
			defProcs = eutils.GetNumericArg(args, "Default number of processors", ncpu, 1, ncpu)
			args = args[1:]
		// performance tuning flags
		case "-proc":
			numProcs = eutils.GetNumericArg(args, "Number of processors", ncpu, 1, ncpu)
			args = args[1:]
		case "-cons":
			serverRatio = eutils.GetNumericArg(args, "Parser to processor ratio", 4, 1, 32)
			args = args[1:]
		case "-serv":
			numServe = eutils.GetNumericArg(args, "Concurrent parser count", 0, 1, 128)
			args = args[1:]
		case "-chan":
			chanDepth = eutils.GetNumericArg(args, "Communication channel depth", 0, ncpu, 128)
			args = args[1:]
		case "-heap":
			heapSize = eutils.GetNumericArg(args, "Unshuffler heap size", 8, 8, 64)
			args = args[1:]
		case "-farm":
			farmSize = eutils.GetNumericArg(args, "Node buffer length", 4, 4, 2048)
			args = args[1:]
		case "-gogc":
			goGc = eutils.GetNumericArg(args, "Garbage collection percentage", 0, 50, 1000)
			args = args[1:]

		// read data from file
		case "-input":
			fileName = eutils.GetStringArg(args, "Input file name")
			args = args[1:]

		// input is indexed with <NEXT_RECORD_SIZE> objects
		case "-turbo":
			turbo = true

		// data cleanup flags
		case "-compress", "-compressed":
			doCompress = true
		case "-spaces", "-cleanup":
			doCleanup = true
		case "-strict":
			doStrict = true
		case "-mixed":
			doMixed = true
		case "-self":
			doSelf = true
		case "-accent":
			deAccent = true
		case "-symbol":
			deSymbol = true
		case "-ascii":
			doASCII = true

		// previously visible processing flags (undocumented)
		case "-stems", "-stem":
			doStem = true
		case "-stops", "-stop":
			deStop = false

		// allow setting of unicode, script, and mathml flags (undocumented)
		case "-unicode":
			// unicodePolicy = GetStringArg(args, "Unicode argument")
			args = args[1:]
		case "-script":
			// scriptPolicy = GetStringArg(args, "Script argument")
			args = args[1:]
		case "-mathml":
			// mathmlPolicy = GetStringArg(args, "MathML argument")
			args = args[1:]

		case "-flag", "-flags":
			flgs = eutils.GetStringArg(args, "Flags argument")
			args = args[1:]

		// debugging flags
		case "-debug":
			// dbug = true
		case "-empty":
			mpty = true
		case "-ident":
			idnt = true
		case "-stats", "-stat":
			stts = true
		case "-timer":
			timr = true
		case "-profile":
			prfl = true
		case "-trial", "-trials":
			trial = true

		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past argument
		args = args[1:]

		if len(args) < 1 {
			break
		}
	}

	// -flag allows script to set -strict or -mixed (or -stems, or -stops) from argument
	switch flgs {
	case "strict":
		doStrict = true
	case "mixed":
		doMixed = true
	case "stems", "stem":
		doStem = true
	case "stops", "stop":
		deStop = false
	case "none", "default":
	default:
		if flgs != "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized -flag value '%s'\n", flgs)
			os.Exit(1)
		}
	}

	/*
		UnicodeFix = parseMarkup(unicodePolicy, "-unicode")
		ScriptFix = parseMarkup(scriptPolicy, "-script")
		MathMLFix = parseMarkup(mathmlPolicy, "-mathml")

		if UnicodeFix != NOMARKUP {
			doUnicode = true
		}

		if ScriptFix != NOMARKUP {
			doScript = true
		}

		if MathMLFix != NOMARKUP {
			doMathML = true
		}
	*/

	if numProcs == 0 {
		if defProcs > 0 {
			numProcs = defProcs
		} else if maxProcs > 0 {
			numProcs = maxProcs
		}
	}
	if numProcs > ncpu {
		numProcs = ncpu
	}
	if numProcs > maxProcs {
		numProcs = maxProcs
	}

	eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, turbo)

	eutils.SetOptions(doStrict, doMixed, doSelf, deAccent, deSymbol, doASCII, doCompress, doCleanup, doStem, deStop)

	// -stats prints number of CPUs and performance tuning values if no other arguments (undocumented)
	if stts && len(args) < 1 {

		eutils.PrintStats()

		return
	}

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// DOCUMENTATION COMMANDS

	inSwitch = true

	switch args[0] {
	case "-version":
		fmt.Printf("%s\n", eutils.EDirectVersion)
	case "-help", "help":
		eutils.PrintHelp("xtract", "xtract-help.txt")
	case "-examples", "-example":
		eutils.PrintHelp("xtract", "xtract-examples.txt")
	case "-extras", "-extra", "-advanced":
		fmt.Printf("Please run rchive -help for local record indexing information\n")
	case "-internal", "-internals":
		eutils.PrintHelp("xtract", "xtract-internal.txt")
	case "-keys":
		eutils.PrintHelp("xtract", "xtract-keys.txt")
	case "-unix":
		eutils.PrintHelp("xtract", "xtract-unix.txt")
	default:
		// if not any of the documentation commands, keep going
		inSwitch = false
	}

	if inSwitch {
		return
	}

	// FILE NAME CAN BE SUPPLIED WITH -input COMMAND

	in := os.Stdin

	// check for data being piped into stdin
	isPipe := false
	fi, err := os.Stdin.Stat()
	if err == nil {
		isPipe = bool((fi.Mode() & os.ModeNamedPipe) != 0)
	}

	usingFile := false

	if fileName != "" {

		inFile, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		defer inFile.Close()

		// use indicated file instead of stdin
		in = inFile
		usingFile = true

		if isPipe && runtime.GOOS != "windows" {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: Input data from both stdin and file '%s', mode is '%s'\n", fileName, mode)
			os.Exit(1)
		}
	}

	// check for -input command after extraction arguments
	for _, str := range args {
		if str == "-input" {
			fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -input command\n")
			os.Exit(1)
		}
	}

	// START PROFILING IF REQUESTED

	if prfl {

		f, err := os.Create("cpu.pprof")
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create profile output file\n")
			os.Exit(1)
		}

		pprof.StartCPUProfile(f)

		defer pprof.StopCPUProfile()
	}

	// INITIALIZE RECORD COUNT

	recordCount := 0
	byteCount := 0

	// print processing rate and program duration
	printDuration := func(name string) {

		eutils.PrintDuration(name, recordCount, byteCount)
	}

	// NAME OF OUTPUT STRING TRANSFORMATION FILE

	tform := ""
	transform := make(map[string]string)

	populateTx := func(tf string) {

		inFile, err := os.Open(tf)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Unable to open transformation file %s\n", err.Error())
			os.Exit(1)
		}
		defer inFile.Close()

		scanr := bufio.NewScanner(inFile)

		// populate transformation map for -translate (and -matrix) output
		for scanr.Scan() {

			line := scanr.Text()
			frst, scnd := eutils.SplitInTwoLeft(line, "\t")

			transform[frst] = scnd
		}
	}

	if len(args) > 2 && args[0] == "-transform" {
		tform = args[1]
		args = args[2:]
		if tform != "" {
			populateTx(tform)
		}
	}

	// SEQUENCE RECORD EXTRACTION COMMAND GENERATOR

	// -insd simplifies extraction of INSDSeq qualifiers
	if args[0] == "-insd" || args[0] == "-insd-" || args[0] == "-insd-idx" {

		addDash := true
		doIndex := false
		// -insd- variant suppresses use of dash as placeholder for missing qualifiers (undocumented)
		if args[0] == "-insd-" {
			addDash = false
		}
		// -insd-idx variant creates word index using -indices command (undocumented)
		if args[0] == "-insd-idx" {
			doIndex = true
			addDash = false
		}

		args = args[1:]

		insd := processINSD(args, isPipe || usingFile, addDash, doIndex)

		if !isPipe && !usingFile {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range insd {
				fmt.Printf(" %s", str)
			}
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = insd
	}

	// CITATION MATCHER EXTRACTION COMMAND GENERATOR

	// -citmatch extracts PMIDs from nquire -citmatch output (undocumented)
	if args[0] == "-citmatch" {

		var acc []string

		acc = append(acc, "-pattern", "opt")
		if isPipe {
			acc = append(acc, "-sep", "\n")
		} else {
			acc = append(acc, "-sep", "\"\\n\"")
		}
		acc = append(acc, "-element", "uids/pubmed")

		if !isPipe && !usingFile {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range acc {
				fmt.Printf(" %s", str)
			}
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = acc
	}

	// BIOTHINGS EXTRACTION COMMAND GENERATOR

	// -biopath takes a parent object and a dotted exploration path for BioThings resources (undocumented)
	if args[0] == "-biopath" {

		args = args[1:]

		biopath := processBiopath(args, isPipe || usingFile)

		if !isPipe && !usingFile {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range biopath {
				fmt.Printf(" %s", str)
			}
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = biopath
	}

	// SPECIFY STRINGS TO GO BEFORE AND AFTER ENTIRE OUTPUT OR EACH RECORD

	head := ""
	tail := ""

	hd := ""
	tl := ""

	for {

		inSwitch = true

		switch args[0] {
		case "-head":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -head command\n")
				os.Exit(1)
			}
			head = eutils.ConvertSlash(args[1])
			// allow splitting of -head argument, keep appending until next command (undocumented)
			ofs, nxt := 0, args[2:]
			for {
				if len(nxt) < 1 {
					break
				}
				tmp := nxt[0]
				if strings.HasPrefix(tmp, "-") {
					break
				}
				ofs++
				txt := eutils.ConvertSlash(tmp)
				if head != "" && !strings.HasSuffix(head, "\t") {
					head += "\t"
				}
				head += txt
				nxt = nxt[1:]
			}
			if ofs > 0 {
				args = args[ofs:]
			}
		case "-tail":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -tail command\n")
				os.Exit(1)
			}
			tail = eutils.ConvertSlash(args[1])
		case "-hd":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -hd command\n")
				os.Exit(1)
			}
			hd = eutils.ConvertSlash(args[1])
		case "-tl":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -tl command\n")
				os.Exit(1)
			}
			tl = eutils.ConvertSlash(args[1])
		case "-wrp":
			// shortcut to wrap records in XML tags
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -wrp command\n")
				os.Exit(1)
			}
			tmp := eutils.ConvertSlash(args[1])
			lft, rgt := eutils.SplitInTwoLeft(tmp, ",")
			if lft != "" {
				head = "<" + lft + ">"
				tail = "</" + lft + ">"
			}
			if rgt != "" {
				hd = "<" + rgt + ">"
				tl = "</" + rgt + ">"
			}
		case "-set":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -set command\n")
				os.Exit(1)
			}
			tmp := eutils.ConvertSlash(args[1])
			if tmp != "" {
				head = "<" + tmp + ">"
				tail = "</" + tmp + ">"
			}
		case "-rec":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -rec command\n")
				os.Exit(1)
			}
			tmp := eutils.ConvertSlash(args[1])
			if tmp != "" {
				hd = "<" + tmp + ">"
				tl = "</" + tmp + ">"
			}
		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past arguments
		args = args[2:]

		if len(args) < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract\n")
			os.Exit(1)
		}
	}

	// CREATE XML BLOCK READER FROM STDIN OR FILE

	const FirstBuffSize = 4096

	getFirstBlock := func() string {

		buffer := make([]byte, FirstBuffSize)
		n, err := in.Read(buffer)
		if err != nil && err != io.EOF {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to read first block: %s\n", err.Error())
			// os.Exit(1)
		}
		bufr := buffer[:n]
		return string(bufr)
	}

	first := getFirstBlock()

	mlt := io.MultiReader(strings.NewReader(first), in)

	isJsn := false
	isAsn := false
	isGbf := false
	matched := 0

	// auto-detect XML, JSON, or ASN.1 format
	if first != "" {
		posJ1 := strings.Index(first, "{")
		posJ2 := strings.Index(first, "\":")
		if posJ1 >= 0 && posJ2 >= 0 && posJ1 < posJ2 {
			isJsn = true
			matched++
		} else {
			posJ1 = FirstBuffSize
			posJ2 = FirstBuffSize
		}
		posA1 := strings.Index(first, "::=")
		posA2 := strings.Index(first, "{")
		if posA1 >= 0 && posA2 >= 0 && posA1 < posA2 {
			isAsn = true
			matched++
		} else {
			posA1 = FirstBuffSize
			posA2 = FirstBuffSize
		}
		posG1 := strings.Index(first, "LOCUS")
		posG2 := strings.Index(first, "DEFINITION")
		if posG1 >= 0 && posG2 >= 0 && posG1 < posG2 {
			isGbf = true
			matched++
		} else {
			posG1 = FirstBuffSize
			posG2 = FirstBuffSize
		}
		posX1 := strings.Index(first, "<")
		posX2 := strings.Index(first, ">")
		if posX1 >= 0 && posX2 >= 0 && posX1 < posX2 {
			matched++
		} else {
			posX1 = FirstBuffSize
			posX2 = FirstBuffSize
		}
		if matched > 1 {
			if posX1 < posJ1 && posX1 < posA1 && posX1 < posG1 {
				isJsn = false
				isAsn = false
				isGbf = false
			} else if posJ1 < posA1 && posJ1 < posG1 {
				isAsn = false
				isGbf = false
			} else if posA1 < posJ1 && posA1 < posG1 {
				isJsn = false
				isGbf = false
			} else if posG1 < posJ1 && posG1 < posA1 {
				isJsn = false
				isAsn = false
			}
		}
	}

	if isJsn {
		jrdr := eutils.JSONConverter(mlt, "root", "", "element")
		mlt = eutils.ChanToReader(jrdr)
	} else if isAsn {
		ardr := eutils.ASN1Converter(mlt, "", "")
		mlt = eutils.ChanToReader(ardr)
	} else if isGbf {
		grdr := eutils.GenBankConverter(mlt)
		mlt = eutils.ChanToReader(grdr)
	}

	rdr := eutils.CreateXMLStreamer(mlt)
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
		os.Exit(1)
	}

	// CONFIRM INPUT DATA AVAILABILITY AFTER RUNNING COMMAND GENERATORS

	if fileName == "" && runtime.GOOS != "windows" {

		fromStdin := bool((fi.Mode() & os.ModeCharDevice) == 0)
		if !isPipe || !fromStdin {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: No data supplied to xtract from stdin or file, mode is '%s'\n", mode)
			os.Exit(1)
		}
	}

	if !usingFile && !isPipe {

		fmt.Fprintf(os.Stderr, "\nERROR: No XML input data supplied to xtract\n")
		os.Exit(1)
	}

	// XML VALIDATION

	nextArg := func() (string, bool) {

		if len(args) < 1 {
			return "", false
		}

		// remove next token from slice
		nxt := args[0]
		args = args[1:]

		return nxt, true
	}

	if args[0] == "-verify" || args[0] == "-validate" {

		// skip past command name
		args = args[1:]

		find := ""
		html := false

		// look for optional arguments
		for {
			arg, ok := nextArg()
			if !ok {
				break
			}

			switch arg {
			case "-find":
				// override set wrapper
				find, ok = nextArg()
			case "-html":
				html = true
			}
		}

		recordCount = eutils.ValidateXML(rdr, find, html)

		debug.FreeOSMemory()

		// suppress printing of lines if not properly counted
		if recordCount == 1 {
			recordCount = 0
		}

		if timr {
			printDuration("lines")
		}

		return
	}

	// MISCELLANEOUS TIMING COMMANDS

	if args[0] == "-chunk" {

		for str := range rdr {
			recordCount++
			byteCount += len(str)
		}

		printDuration("blocks")

		return
	}

	if args[0] == "-split" {

		if len(args) > 1 {
			if args[1] == "-pattern" {
				// skip past -split if followed by -pattern
				args = args[1:]
			}
		}
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -split command\n")
			os.Exit(1)
		}
		pat := args[1]

		eutils.PartitionPattern(pat, "", turbo, rdr,
			func(str string) {
				recordCount++
				byteCount += len(str)
			})

		printDuration("patterns")

		return
	}

	if args[0] == "-token" {

		eutils.StreamTokens(rdr,
			func(tkn eutils.XMLToken) {
				recordCount++
				byteCount += len(tkn.Name) + len(tkn.Attr)
			})

		printDuration("tokens")

		return
	}

	// INDEXED XML FILE PREPARATION

	// cat carotene.xml | xtract -timer -index -pattern PubmedArticle > carindex.txt
	// xtract -timer -turbo -input carindex.txt -pattern PubmedArticle -element LastName
	if args[0] == "-index" {

		if len(args) > 1 {
			if args[1] == "-pattern" {
				// skip past -index if followed by -pattern
				args = args[1:]
			}
		}
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -index command\n")
			os.Exit(1)
		}
		pat := args[1]

		retlength := len("\n")

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		eutils.PartitionPattern(pat, "", false, rdr,
			func(str string) {
				recordCount++
				nxt := len(str)
				byteCount += nxt
				newln := false

				if !strings.HasSuffix(str, "\n") {
					nxt += retlength
					newln = true
				}

				os.Stdout.WriteString("<NEXT_RECORD_SIZE>")
				val := strconv.Itoa(nxt)
				os.Stdout.WriteString(val)
				os.Stdout.WriteString("</NEXT_RECORD_SIZE>\n")

				os.Stdout.WriteString(str)
				if newln {
					os.Stdout.WriteString("\n")
				}
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// ENSURE PRESENCE OF PATTERN ARGUMENT

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// allow -record as synonym of -pattern (undocumented)
	if args[0] == "-record" || args[0] == "-Record" {
		args[0] = "-pattern"
	}

	// make sure top-level -pattern command is next
	if args[0] != "-pattern" && args[0] != "-Pattern" {
		fmt.Fprintf(os.Stderr, "\nERROR: No -pattern in command-line arguments\n")
		os.Exit(1)
	}
	if len(args) < 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: Item missing after -pattern command\n")
		os.Exit(1)
	}

	topPat := args[1]
	if topPat == "" {
		fmt.Fprintf(os.Stderr, "\nERROR: Item missing after -pattern command\n")
		os.Exit(1)
	}
	if strings.HasPrefix(topPat, "-") {
		fmt.Fprintf(os.Stderr, "\nERROR: Misplaced %s command\n", topPat)
		os.Exit(1)
	}

	// look for -pattern Parent/* construct for heterogeneous data, e.g., -pattern PubmedArticleSet/*
	topPattern, star := eutils.SplitInTwoLeft(topPat, "/")
	if topPattern == "" {
		return
	}

	parent := ""
	if star == "*" {
		parent = topPattern
	} else if star != "" {
		fmt.Fprintf(os.Stderr, "\nERROR: -pattern Parent/Child construct is not supported\n")
		os.Exit(1)
	}

	// SAVE ONLY RECORDS WITH NON-ASCII CHARACTERS

	// -pattern record_name -select -nonascii
	if len(args) == 4 && args[2] == "-select" && args[3] == "-nonascii" {

		xmlq := eutils.CreateXMLProducer(topPattern, star, false, rdr)
		fchq := createUnicoders(xmlq)
		unsq := eutils.CreateXMLUnshuffler(fchq)

		if xmlq == nil || fchq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create selector\n")
			os.Exit(1)
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		// drain output channel
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			// send result to output
			os.Stdout.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				os.Stdout.WriteString("\n")
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}

			recordCount++
			runtime.Gosched()
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// READ FILE OF IDENTIFIERS AND CONCURRENTLY EXTRACT SELECTED RECORDS

	// -pattern record_name -select parent/element@attribute^version -in file_of_identifiers
	if len(args) == 6 && args[2] == "-select" && (args[4] == "-in" || args[4] == "-retaining") {

		indx := args[3]
		unqe := args[5]

		// read file of identifiers to use for filtering
		fl, err := os.Open(unqe)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open identifier file '%s'\n", unqe)
			os.Exit(1)
		}

		// create map that records each UID
		order := make(map[string]bool)

		scanr := bufio.NewScanner(fl)

		// read lines of identifiers
		for scanr.Scan() {

			line := scanr.Text()
			id, _ := eutils.SplitInTwoLeft(line, "\t")

			id = sortStringByWords(id)

			// add identifier to map
			order[id] = true
		}

		fl.Close()

		xmlq := eutils.CreateXMLProducer(topPattern, star, false, rdr)
		fchq := createSelectors(topPattern, indx, order, xmlq)
		unsq := eutils.CreateXMLUnshuffler(fchq)

		if xmlq == nil || fchq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create selector\n")
			os.Exit(1)
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		// drain output channel
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			// send result to output
			os.Stdout.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				os.Stdout.WriteString("\n")
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}

			recordCount++
			runtime.Gosched()
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// READ FILE OF IDENTIFIERS AND EXCLUDE SELECTED RECORDS

	// -pattern record_name -exclude element -excluding file_of_identifiers (undocumented)
	if len(args) == 6 && args[2] == "-select" && args[4] == "-excluding" {

		indx := args[3]
		unqe := args[5]

		// read file of identifiers to use for filtering
		fl, err := os.Open(unqe)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open identifier file '%s'\n", unqe)
			os.Exit(1)
		}

		// create map that records each UID
		order := make(map[string]bool)

		scanr := bufio.NewScanner(fl)

		// read lines of identifiers
		for scanr.Scan() {

			line := scanr.Text()
			id, _ := eutils.SplitInTwoLeft(line, "\t")
			id = strings.ToLower(id)

			// add identifier to map
			order[id] = true
		}

		fl.Close()

		find := eutils.ParseIndex(indx)

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				id := eutils.FindIdentifier(str[:], parent, find)
				if id != "" {
					id = strings.ToLower(id)
					_, ok := order[id]
					if ok {
						// in exclusion list, skip
						return
					}
				}

				if hd != "" {
					os.Stdout.WriteString(hd)
					os.Stdout.WriteString("\n")
				}

				// write selected record
				os.Stdout.WriteString(str[:])
				os.Stdout.WriteString("\n")

				if tl != "" {
					os.Stdout.WriteString(tl)
					os.Stdout.WriteString("\n")
				}
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// READ ORDERED FILE OF IDENTIFIERS AND XML STRINGS, APPEND XML JUST INSIDE CLOSING TAG OF APPROPRIATE RECORD

	// -pattern record_name -select element -appending file_of_identifiers_and_metadata (undocumented)
	if len(args) == 6 && args[2] == "-select" && args[4] == "-appending" {

		indx := args[3]
		apnd := args[5]

		fl, err := os.Open(apnd)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open transformation file '%s'\n", apnd)
			os.Exit(1)
		}

		scanr := bufio.NewScanner(fl)

		find := eutils.ParseIndex(indx)

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		rgt := "</" + topPattern + ">"

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				id := eutils.FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}
				id = strings.ToLower(id)

				for scanr.Scan() {

					line := scanr.Text()
					frst, scnd := eutils.SplitInTwoLeft(line, "\t")
					frst = strings.ToLower(frst)

					if id != frst {
						return
					}
					if !strings.HasSuffix(str, rgt) {
						return
					}

					lft := strings.TrimSuffix(str, rgt)
					str = lft + "  " + scnd + "\n" + rgt

					if hd != "" {
						os.Stdout.WriteString(hd)
						os.Stdout.WriteString("\n")
					}

					os.Stdout.WriteString(str[:])
					os.Stdout.WriteString("\n")

					if tl != "" {
						os.Stdout.WriteString(tl)
						os.Stdout.WriteString("\n")
					}

					break
				}
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		fl.Close()

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// SORT XML RECORDS BY IDENTIFIER

	// -pattern record_name -sort parent/element@attribute^version
	if len(args) == 4 && args[2] == "-sort" {

		indx := args[3]

		// create map that records each UID
		order := make(map[string][]string)

		find := eutils.ParseIndex(indx)

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				id := eutils.FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				data, ok := order[id]
				if !ok {
					data = make([]string, 0, 1)
				}
				data = append(data, str)
				// always need to update order, since data may be reallocated
				order[id] = data
			})

		var keys []string
		for ky := range order {
			keys = append(keys, ky)
		}
		// sort fields in alphabetical or numeric order
		sort.Slice(keys, func(i, j int) bool {
			// numeric sort on strings checks lengths first
			if eutils.IsAllDigits(keys[i]) && eutils.IsAllDigits(keys[j]) {
				lni := len(keys[i])
				lnj := len(keys[j])
				// shorter string is numerically less, assuming no leading zeros
				if lni < lnj {
					return true
				}
				if lni > lnj {
					return false
				}
			}
			// same length or non-numeric, can now do string comparison on contents
			return keys[i] < keys[j]
		})

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		for _, id := range keys {

			strs := order[id]
			for _, str := range strs {
				os.Stdout.WriteString(str)
				os.Stdout.WriteString("\n")
			}
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// SPLIT FILE BY BY RECORD COUNT

	// split XML record into subfiles by count
	if len(args) == 8 && args[2] == "-split" && args[4] == "-prefix" && args[6] == "-suffix" {

		// e.g., -head "<IdxDocumentSet>" -tail "</IdxDocumentSet>" -pattern IdxDocument -split 250000 -prefix "biocon" -suffix "e2x"
		count := 0
		fnum := 0
		var (
			fl  *os.File
			err error
		)
		chunk, err := strconv.Atoi(args[3])
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return
		}
		prefix := args[5]
		suffix := args[7]

		eutils.PartitionPattern(topPattern, star, false, rdr,
			func(str string) {
				recordCount++

				if count >= chunk {
					if tail != "" {
						fl.WriteString(tail)
						fl.WriteString("\n")
					}
					fl.Close()
					count = 0
				}
				if count == 0 {
					fpath := fmt.Sprintf("%s%03d.%s", prefix, fnum, suffix)
					fl, err = os.Create(fpath)
					if err != nil {
						fmt.Fprintf(os.Stderr, "%s\n", err.Error())
						return
					}
					os.Stderr.WriteString(fpath + "\n")
					fnum++
					if head != "" {
						fl.WriteString(head)
						fl.WriteString("\n")
					}
				}
				count++

				fl.WriteString(str[:])
				fl.WriteString("\n")
			})

		if count >= chunk {
			if tail != "" {
				fl.WriteString(tail)
				fl.WriteString("\n")
			}
			fl.Close()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// PARSE AND VALIDATE EXTRACTION ARGUMENTS

	// parse nested exploration instruction from command-line arguments
	cmds := eutils.ParseArguments(args, topPattern)
	if cmds == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Problem parsing command-line arguments\n")
		os.Exit(1)
	}

	// GLOBAL MAP FOR SORT-UNIQ-COUNT HISTOGRAM ARGUMENT

	histogram := make(map[string]int)

	// PERFORMANCE TIMING COMMAND

	// -stats with an extraction command prints XML size and processing time for each record
	if stts {

		legend := "REC\tOFST\tSIZE\tTIME"

		rec := 0

		eutils.PartitionPattern(topPattern, star, turbo, rdr,
			func(str string) {
				rec++
				beginTime := time.Now()
				eutils.ProcessExtract(str[:], parent, rec, hd, tl, transform, histogram, cmds)
				endTime := time.Now()
				duration := endTime.Sub(beginTime)
				micro := int(float64(duration.Nanoseconds()) / 1e3)
				if legend != "" {
					fmt.Printf("%s\n", legend)
					legend = ""
				}
				fmt.Printf("%d\t%d\t%d\n", rec, len(str), micro)
			})

		return
	}

	// PERFORMANCE OPTIMIZATION FUNCTION

	// -trial -input fileName runs the specified extraction for each -proc from 1 to nCPU
	if trial && fileName != "" {

		legend := "CPU\tRATE\tDEV"

		for numServ := 1; numServ <= ncpu; numServ++ {

			numServe = numServ

			eutils.SetTunings(numProcs, numServe, serverRatio, chanDepth, farmSize, heapSize, goGc, turbo)

			runtime.GOMAXPROCS(numServ)

			sum := 0
			count := 0
			mean := 0.0
			m2 := 0.0

			// calculate mean and standard deviation of processing rate
			for trials := 0; trials < 5; trials++ {

				inFile, err := os.Open(fileName)
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
					os.Exit(1)
				}

				trdr := eutils.CreateXMLStreamer(inFile)
				if trdr == nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to read input file\n")
					os.Exit(1)
				}

				xmlq := eutils.CreateXMLProducer(topPattern, star, turbo, trdr)
				tblq := eutils.CreateConsumers(cmds, parent, hd, tl, transform, histogram, xmlq)

				if xmlq == nil || tblq == nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to create servers\n")
					os.Exit(1)
				}

				begTime := time.Now()
				recordCount = 0

				for range tblq {
					recordCount++
					runtime.Gosched()
				}

				inFile.Close()

				debug.FreeOSMemory()

				endTime := time.Now()
				expended := endTime.Sub(begTime)
				secs := float64(expended.Nanoseconds()) / 1e9

				if secs >= 0.000001 && recordCount > 0 {
					speed := int(float64(recordCount) / secs)
					sum += speed
					count++
					x := float64(speed)
					delta := x - mean
					mean += delta / float64(count)
					m2 += delta * (x - mean)
				}
			}

			if legend != "" {
				fmt.Printf("%s\n", legend)
				legend = ""
			}
			if count > 1 {
				vrc := m2 / float64(count-1)
				dev := int(math.Sqrt(vrc))
				fmt.Printf("%d\t%d\t%d\n", numServ, sum/count, dev)
			}
		}

		return
	}

	// PROCESS SINGLE SELECTED RECORD IF -pattern ARGUMENT IS IMMEDIATELY FOLLOWED BY -position COMMAND

	posn := ""
	if cmds.Visit == topPat {
		if cmds.Position == "outer" ||
			cmds.Position == "inner" ||
			cmds.Position == "even" ||
			cmds.Position == "odd" ||
			cmds.Position == "all" {
			// filter by record position when draining unshuffler channel
			posn = cmds.Position
			cmds.Position = ""
		}
	}

	if cmds.Visit == topPat && cmds.Position != "" && cmds.Position != "select" {

		qry := ""
		idx := 0
		rec := 0

		if cmds.Position == "first" {

			eutils.PartitionPattern(topPattern, star, turbo, rdr,
				func(str string) {
					rec++
					if rec == 1 {
						qry = str
						idx = rec
					}
				})

		} else if cmds.Position == "last" {

			eutils.PartitionPattern(topPattern, star, turbo, rdr,
				func(str string) {
					qry = str
					idx = rec
				})

		} else {

			// use numeric position
			number, err := strconv.Atoi(cmds.Position)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized position '%s'\n", cmds.Position)
				os.Exit(1)
			}

			eutils.PartitionPattern(topPattern, star, turbo, rdr,
				func(str string) {
					rec++
					if rec == number {
						qry = str
						idx = rec
					}
				})
		}

		if qry == "" {
			return
		}

		// clear position on top node to prevent condition test failure
		cmds.Position = ""

		// process single selected record
		res := eutils.ProcessExtract(qry[:], parent, idx, hd, tl, transform, histogram, cmds)

		if res != "" {
			fmt.Printf("%s", res)
		}

		return
	}

	// LAUNCH PRODUCER, CONSUMER, AND UNSHUFFLER GOROUTINES

	// launch producer goroutine to partition XML by pattern
	xmlq := eutils.CreateXMLProducer(topPattern, star, turbo, rdr)

	// launch consumer goroutines to parse and explore partitioned XML objects
	tblq := eutils.CreateConsumers(cmds, parent, hd, tl, transform, histogram, xmlq)

	// launch unshuffler goroutine to restore order of results
	unsq := eutils.CreateXMLUnshuffler(tblq)

	if xmlq == nil || tblq == nil || unsq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create servers\n")
		os.Exit(1)
	}

	// PERFORMANCE SUMMARY

	/*
		if dbug {

			// drain results, but suppress extraction output
			for ext := range unsq {
				byteCount += len(ext.Text)
				recordCount++
				runtime.Gosched()
			}

			// force garbage collection, return memory to operating system
			debug.FreeOSMemory()

			// print processing parameters as XML object
			stopTime := time.Now()
			duration := stopTime.Sub(StartTime)
			seconds := float64(duration.Nanoseconds()) / 1e9

			// Threads is a more easily explained concept than GOMAXPROCS
			fmt.Printf("<Xtract>\n")
			fmt.Printf("  <Threads>%d</Threads>\n", numProcs)
			fmt.Printf("  <Parsers>%d</Parsers>\n", NumServe)
			fmt.Printf("  <Time>%.3f</Time>\n", seconds)
			if seconds >= 0.001 && recordCount > 0 {
				rate := int(float64(recordCount) / seconds)
				fmt.Printf("  <Rate>%d</Rate>\n", rate)
			}
			fmt.Printf("</Xtract>\n")

			return
		}
	*/

	// DRAIN OUTPUT CHANNEL TO EXECUTE EXTRACTION COMMANDS, RESTORE OUTPUT ORDER WITH HEAP

	var buffer strings.Builder
	count := 0
	okay := false
	lastTime := time.Now()

	wrtr := bufio.NewWriter(os.Stdout)

	// printResult prints output for current pattern, handles -empty and -ident flags, and periodically flushes buffer
	printResult := func(curr eutils.XMLRecord) {

		str := curr.Text

		if mpty {

			if str == "" {

				okay = true

				idx := curr.Index
				val := strconv.Itoa(idx)
				buffer.WriteString(val[:])
				buffer.WriteString("\n")

				count++
			}

		} else if str != "" {

			okay = true

			if idnt {
				idx := curr.Index
				val := strconv.Itoa(idx)
				buffer.WriteString(val[:])
				buffer.WriteString("\t")
			}

			// save output to byte buffer
			buffer.WriteString(str[:])

			count++
		}

		thisTime := time.Now()
		duration := thisTime.Sub(lastTime)
		milliSeconds := duration.Milliseconds()

		if count > 1000 || milliSeconds > 4999 {
			count = 0
			lastTime = thisTime
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
			buffer.Reset()
		}
	}

	if head != "" {
		buffer.WriteString(head[:])
		buffer.WriteString("\n")
	}

	// drain unshuffler channel

	if posn == "outer" {

		// print only first and last records
		var beg *eutils.XMLRecord
		var end *eutils.XMLRecord

		for curr := range unsq {

			if beg == nil {
				beg = &eutils.XMLRecord{Index: curr.Index, Ident: curr.Ident, Text: curr.Text}
			} else {
				end = &eutils.XMLRecord{Index: curr.Index, Ident: curr.Ident, Text: curr.Text}
			}

			recordCount++
		}

		if beg != nil {
			printResult(*beg)
		}
		if end != nil {
			printResult(*end)
		}

	} else if posn == "inner" {

		// print all but first and last records
		var prev *eutils.XMLRecord
		var next *eutils.XMLRecord
		first := true

		for curr := range unsq {

			if first {
				first = false
			} else {
				prev = next
				next = &eutils.XMLRecord{Index: curr.Index, Ident: curr.Ident, Text: curr.Text}
			}

			if prev != nil {
				printResult(*prev)
			}

			recordCount++
		}

	} else if posn == "even" {

		even := false

		for curr := range unsq {

			if even {
				printResult(curr)
			}
			even = !even

			recordCount++
		}

	} else if posn == "odd" {

		odd := true

		for curr := range unsq {

			if odd {
				printResult(curr)
			}
			odd = !odd

			recordCount++
		}

	} else {

		// default or -position all
		for curr := range unsq {

			// send result to output
			printResult(curr)

			recordCount++
			runtime.Gosched()
		}
	}

	if tail != "" {
		buffer.WriteString(tail[:])
		buffer.WriteString("\n")
	}

	// do not print head or tail if no extraction output
	if okay {
		txt := buffer.String()
		if txt != "" {
			// print final buffer
			wrtr.WriteString(txt[:])
		}
	}
	buffer.Reset()

	wrtr.Flush()

	// print -histogram results, if populated
	var keys []string
	for ky := range histogram {
		keys = append(keys, ky)
	}
	if len(keys) > 0 {
		// sort fields in alphabetical or numeric order
		sort.Slice(keys, func(i, j int) bool {
			// numeric sort on strings checks lengths first
			if eutils.IsAllDigits(keys[i]) && eutils.IsAllDigits(keys[j]) {
				lni := len(keys[i])
				lnj := len(keys[j])
				// shorter string is numerically less, assuming no leading zeros
				if lni < lnj {
					return true
				}
				if lni > lnj {
					return false
				}
			}
			// same length or non-numeric, can now do string comparison on contents
			return keys[i] < keys[j]
		})

		for _, str := range keys {

			count := histogram[str]
			val := strconv.Itoa(count)
			os.Stdout.WriteString(val)
			os.Stdout.WriteString("\t")
			os.Stdout.WriteString(str)
			os.Stdout.WriteString("\n")
		}
	}

	// force garbage collection and return memory before calculating processing rate
	debug.FreeOSMemory()

	if timr {
		printDuration("records")
	}
}
