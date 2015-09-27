//----------------------------------------------------------------------------------------
// IVC: ivc.go - Main program.
// Copyright 2015 Nam Sy Vo.
//----------------------------------------------------------------------------------------

package main

import (
	"flag"
	"github.com/namsyvo/IVC"
	"log"
	"os"
	"path"
	"runtime"
	"runtime/pprof"
)

func main() {

	//Starting program----------------------------------------------------------//
	log.Printf("IVC - Integrated Variant Caller using next-generation sequencing data.")
	log.Printf("IVC-main: Calling variants based on alignment between reads and reference multi-genomes.")

	input_info := ReadInputInfo()
	ivc.MEM_STATS = new(runtime.MemStats)

	//Initializing indexes and parameters--------------------------------------//
	variant_caller := ivc.NewVariantCaller(input_info)
	//-------------------------------------------------------------------------//

	//Call variants from read-multigenome alignment----------------------------//
	variant_caller.CallVariants()
	//-------------------------------------------------------------------------//

	//Output variant calls-----------------------------------------------------//
	variant_caller.OutputVarCalls()
	//-------------------------------------------------------------------------//

	log.Printf("Finish whole variant calling process.")
}

//--------------------------------------------------------------------------------------------------
// Read input information and parameters
//--------------------------------------------------------------------------------------------------
func ReadInputInfo() *ivc.InputInfo {
	var genome_file = flag.String("R", "", "reference genome file")
	var var_prof_file = flag.String("V", "", "variant profile file")
	var idx_dir = flag.String("I", "", "index directory")
	var read_file_1 = flag.String("1", "", "pairend read file, first end")
	var read_file_2 = flag.String("2", "", "pairend read file, second end")
	var var_call_file = flag.String("O", "", "variant call output file")
	var search_mode = flag.Int("mode", 1, "searching mode for finding seeds (1: random, 2: deterministic)")
	var start_pos = flag.Int("start", 0, "starting position on reads for finding seeds")
	var search_step = flag.Int("step", 5, "step for searching in deterministic mode")
	var proc_num = flag.Int("t", 0, "maximum number of CPUs")
	var max_snum = flag.Int("maxs", 512, "maximum number of seeds")
	var max_psnum = flag.Int("maxp", 128, "maximum number of paired-seeds")
	var min_slen = flag.Int("lmin", 15, "minimum length of seeds")
	var max_slen = flag.Int("lmax", 25, "maximum length of seeds")
	var dist_thres = flag.Int("d", 36, "threshold of alignment distances")
	var prob_thres = flag.Float64("p", 36.0, "threshold of alignment probabilities")
	var iter_num = flag.Int("r", 12, "maximum number of iterations")
	var sub_cost = flag.Float64("s", 4.0, "substitution cost")
	var gap_open = flag.Float64("o", 4.1, "gap open cost")
	var gap_ext = flag.Float64("e", 1.0, "gap extension cost")
	var debug_mode = flag.Bool("debug", false, "turn on debug mode.")
	flag.Parse()

	_, genome_file_name := path.Split(*genome_file)
	multigenome_file_name := path.Join(*idx_dir, genome_file_name) + ".mgf"
	rev_multigenome_file_name := path.Join(*idx_dir, genome_file_name) + "_rev.mgf"
	_, var_prof_file_name := path.Split(*var_prof_file)
	var_prof_index_file_name := path.Join(*idx_dir, var_prof_file_name) + ".idx"

	input_info := new(ivc.InputInfo)
	input_info.Ref_file = multigenome_file_name
	input_info.Var_prof_file = var_prof_index_file_name
	input_info.Index_file = multigenome_file_name + ".index/"
	input_info.Rev_index_file = rev_multigenome_file_name + ".index/"
	input_info.Read_file_1 = *read_file_1
	input_info.Read_file_2 = *read_file_2
	input_info.Var_call_file = *var_call_file
	input_info.Search_mode = *search_mode
	input_info.Start_pos = *start_pos
	input_info.Search_step = *search_step
	input_info.Proc_num = *proc_num
	input_info.Max_snum = *max_snum
	input_info.Max_psnum = *max_psnum
	input_info.Min_slen = *min_slen
	input_info.Max_slen = *max_slen
	input_info.Dist_thres = *dist_thres
	input_info.Prob_thres = *prob_thres
	input_info.Iter_num = *iter_num
	input_info.Sub_cost = *sub_cost
	input_info.Gap_open = *gap_open
	input_info.Gap_ext = *gap_ext
	input_info.Debug_mode = *debug_mode

	if input_info.Proc_num == 0 {
		input_info.Proc_num = runtime.NumCPU()
		log.Printf("No input for number of threads, set to maximum number of CPUs of the current machine (=%d)", input_info.Proc_num)
	}
	if input_info.Debug_mode {
		var err error
		ivc.CPU_FILE, err = os.Create(input_info.Var_call_file + ".cprof")
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(ivc.CPU_FILE)
		defer pprof.StopCPUProfile()

		ivc.MEM_FILE, err = os.Create(input_info.Var_call_file + ".mprof")
		if err != nil {
			log.Fatal(err)
		}
		defer ivc.MEM_FILE.Close()
		log.Printf("Debug mode:\tCpu_prof_file=%s, Mem_prof_file=%s", input_info.Var_call_file+".cprof", input_info.Var_call_file+".mprof")
	}
	log.Printf("Input files:\tGenome_file: %s, Var_file: %s, Index_file=%s, Rev_index_file=%s,"+
		"Read_file_1=%s, Read_file_2=%s, Var_call_file=%s", input_info.Ref_file, input_info.Var_prof_file,
		input_info.Index_file, input_info.Rev_index_file, input_info.Read_file_1, input_info.Read_file_2, input_info.Var_call_file)
	log.Printf("Input paras:\tSearch_mode=%d, Start_pos=%d, Search_step=%d, Proc_num=%d, Max_snum=%d, Max_psnum=%d, "+
		"Min_slen=%d, Max_slen=%d, Dist_thres=%d, Prob_thres=%.5f, Iter_num=%d, Sub_cost=%.5f, Gap_open=%.5f, Gap_ext=%.5f, Debug_mode=%t",
		input_info.Search_mode, input_info.Start_pos, input_info.Search_step, input_info.Proc_num, input_info.Max_snum, input_info.Max_psnum,
		input_info.Min_slen, input_info.Max_slen, input_info.Dist_thres, input_info.Prob_thres, input_info.Iter_num,
		input_info.Sub_cost, input_info.Gap_open, input_info.Gap_ext, input_info.Debug_mode)

	return input_info
}
