//----------------------------------------------------------------------------------------
// IVC: ivc.go
// Main program for variant calling module.
// Copyright 2015 Nam Sy Vo.
//----------------------------------------------------------------------------------------

package main

import (
	"flag"
	"github.com/namsyvo/IVC"
	"log"
	"path"
)

func main() {
	//Starting program----------------------------------------------------------//
	log.Printf("IVC - Integrated Variant Caller using next-generation sequencing data.")
	log.Printf("IVC-main: Calling variants based on alignment between reads and reference multi-genomes.")

	//Setting up all parameters------------------------------------------------//
	input_para_info := ReadParaInfo()
	ivc.Setup(input_para_info)
	//-------------------------------------------------------------------------//

	//Initializing indexes and parameters--------------------------------------//
	variant_caller := ivc.NewVariantCaller()
	//-------------------------------------------------------------------------//

	//Call variants from read-multigenome alignment----------------------------//
	variant_caller.CallVariants()
	//-------------------------------------------------------------------------//

	//Output variant calls-----------------------------------------------------//
	variant_caller.OutputVarCalls()
	//-------------------------------------------------------------------------//
	log.Printf("Finish whole variant calling process.")
}

func ReadParaInfo() *ivc.ParaInfo {
	var genome_file = flag.String("R", "", "reference genome file")
	var var_prof_file = flag.String("V", "", "variant profile file")
	var idx_dir = flag.String("I", "", "index directory")
	var read_file_1 = flag.String("1", "", "pairend read file, first end")
	var read_file_2 = flag.String("2", "", "pairend read file, second end")
	var var_call_file = flag.String("O", "", "variant call output file")
	var search_mode = flag.Int("mode", 0, "searching mode for finding seeds (1: random, 2: deterministic)")
	var start_pos = flag.Int("start", 0, "starting position on reads for finding seeds")
	var search_step = flag.Int("step", 0, "step for searching in deterministic mode")
	var max_snum = flag.Int("maxs", 0, "maximum number of seeds")
	var max_psnum = flag.Int("maxp", 0, "maximum number of paired-seeds")
	var min_slen = flag.Int("lmin", 0, "minimum length of seeds")
	var max_slen = flag.Int("lmax", 0, "maximum length of seeds")
	var dist_thres = flag.Float64("d", 0, "threshold of alignment distances")
	var iter_num = flag.Int("r", 0, "maximum number of iterations")
	var sub_cost = flag.Float64("s", 0, "substitution cost")
	var gap_open = flag.Float64("o", 0, "gap open cost")
	var gap_ext = flag.Float64("e", 0, "gap extension cost")
	var proc_num = flag.Int("t", 0, "maximum number of CPUs")
	var debug_mode = flag.Bool("debug", false, "turn on debug mode.")
	flag.Parse()

	_, genome_file_name := path.Split(*genome_file)
	multi_seq_file_name := path.Join(*idx_dir, genome_file_name) + "_mgf.fasta"
	rev_multi_seq_file_name := path.Join(*idx_dir, genome_file_name) + "_mgf_rev.fasta"
	_, var_prof_file_name := path.Split(*var_prof_file)
	var_prof_index_file_name := path.Join(*idx_dir, var_prof_file_name) + ".idx"

	para_info := new(ivc.ParaInfo)
	para_info.Ref_file = multi_seq_file_name
	para_info.Var_prof_file = var_prof_index_file_name
	para_info.Index_file = multi_seq_file_name + ".index/"
	para_info.Rev_index_file = rev_multi_seq_file_name + ".index/"
	para_info.Read_file_1 = *read_file_1
	para_info.Read_file_2 = *read_file_2
	para_info.Var_call_file = *var_call_file
	para_info.Search_mode = *search_mode
	para_info.Start_pos = *start_pos
	para_info.Search_step = *search_step
	para_info.Max_snum = *max_snum
	para_info.Max_psnum = *max_psnum
	para_info.Min_slen = *min_slen
	para_info.Max_slen = *max_slen
	para_info.Dist_thres = *dist_thres
	para_info.Iter_num = *iter_num
	para_info.Sub_cost = *sub_cost
	para_info.Gap_open = *gap_open
	para_info.Gap_ext = *gap_ext
	para_info.Proc_num = *proc_num
	para_info.Debug_mode = *debug_mode

	return para_info
}
