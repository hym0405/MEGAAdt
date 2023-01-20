#!/usr/bin/env python3

### Patent Pending
### © [2022] The Trustees of Columbia University in the City of New York. 
### Use is subject to the terms of the License Agreement. 

import os
import sys
import argparse
from Bio.SeqUtils import MeltingTemp as mt

def main():

	parser = argparse.ArgumentParser(description = "Design tools for Mutagenesis by Template-guided Amplicon Assembly (MEGAA). " + \
						"See details in https://github.com/hym0405/MEGAAdt", add_help = False)

	group1 = parser.add_argument_group("Input and output [Required]")
	group2 = parser.add_argument_group("Oligo design [Optional]")
	group3 = parser.add_argument_group("Warning information [Optional]")

	group1.add_argument("-f", "--fasta", type = str,
					help = "Input template sequences in FASTA format (See example: ./test_input/input_seq.fasta)")

	group1.add_argument("-i", "--variant", type = str,
					help = "Input variants information in tab-separated table (See example: ./test_input/variant_info.tsv)")

	group1.add_argument("-o", "--output", type = str,
					help = "Prefix of output files; [output].variant.fasta and [output].oligo.tsv will be generated")

	group2.add_argument("-el", "--end_length", type = int, default = 10,
					help="Minimum number of perfectly matched bases (nt) at the end of mutagenesis oligo (both 3' and 5')' [default: 10]")

	group2.add_argument("-gc3", "--gc_clamp3", type = int, default = 1, choices = [0, 1],
					help="# of GC bases (nt) required at 3' end of mutagenesis oligo. The value should be either 1 or 0 [default: 1]")

	group2.add_argument("-Tm5", "--melting5", type = int, default = 35,
					help="Minimum melting temperature (°C) at 5' end of mutagenesis oligos [default: 35]")

	group2.add_argument("-Tm3h", "--melting3_high", type = int, default = 65,
					help="Upper limit of gradient melting temperature (°C) at 3' end of mutagenesis oligos [default: 65]")

	group2.add_argument("-Tm3l", "--melting3_low", type = int, default = 50,
					help="Lower limit of gradient melting temperature (°C) at 3' end of mutagenesis oligos [default: 50]")

	group3.add_argument("-wg", "--warning_gap", type = int, default = 5,
					help="Minimum number of gaps (nt) between oligo to print warnings [default: 5]")

	group3.add_argument("-wl", "--warning_length", type = int, default = 60,
					help="Minimum length of oligos to print warnings [default: 60]")

	group3.add_argument("-h", "--help", action="help", help="show this help message and exit")

	args = parser.parse_args(args = None if sys.argv[1:] else ['--help'])

	p5_Tm_cut = args.melting5
	p3_Tm_high = args.melting3_high
	p3_Tm_low = args.melting3_low
	p3_GCclamp = args.gc_clamp3
	end_extension = args.end_length
	warning_oligo_gap = args.warning_gap
	warining_primer_length = args.warning_length

	internal_paramters = {"Na": 80, "Tris": 50, "Mg": 20}

	reference_fasta_pool = readFASTAfile(args.fasta)
	variant_info_pool, variant_list = readInfoFile(args.variant)
	check_flag = sanityCheck_info(variant_info_pool, reference_fasta_pool)
	if check_flag == "pass":
		variant_primer_pool = determinePrimerCount(variant_info_pool, reference_fasta_pool, end_extension)
		variant_design_pool = designPrimer(variant_primer_pool, variant_info_pool, reference_fasta_pool, end_extension, p3_GCclamp, \
				            p5_Tm_cut, p3_Tm_high, p3_Tm_low, internal_paramters)
		output_primerInfo_pool, output_finalSeq_pool = generateOutput(variant_design_pool, variant_primer_pool, variant_info_pool, \
							reference_fasta_pool, warining_primer_length, warning_oligo_gap)

		f_fasta = open(args.output + ".variant.fasta", "w")
		for each in variant_list:
			f_fasta.writelines([">" + each + os.linesep])
			f_fasta.writelines([output_finalSeq_pool[each] + os.linesep])
		f_fasta.close()

		f_primer = open(args.output + ".oligo.tsv", "w")
		f_primer.writelines(["Variant\tTemplate\tPrimer\tStart_site\tSequence\tLength\tTm_5p\tTm_3p\tMutations\tNotes\tWarnings" + os.linesep])
		for each in variant_list:
			primer_output = output_primerInfo_pool[each]
			for eachP in primer_output:
				f_primer.writelines(["\t".join([str(e) for e in eachP]) + os.linesep])
		f_primer.close()


def readFASTAfile(file):
	pool = {}
	f = open(file, "r")
	data = f.readlines()
	f.close()
	tmpLabel = ""
	tmpSeq = ""
	flag = 0
	for each in data:
		each = each.replace("\n", "")
		each = each.replace("\r", "")
		if each[0] == ">":
			if flag == 1:
				pool[tmpLabel] = tmpSeq
			tmpLabel = each[1:]
			tmpSeq = ""
		else:
			tmpSeq += each.upper()
		flag = 1
	pool[tmpLabel] = tmpSeq
	return pool
def readInfoFile(file):
	pool = {}
	variant_list = []
	f = open(file, "r")
	data = f.readlines()
	f.close()
	for each in data[1:]:
		each = each.replace("\n", "")
		each = each.replace("\r", "")
		if len(each) == 0:
			continue
		if each[0] == "#":
			continue
		tmp = each.split("\t")
		pool[tmp[0]] = []
		if tmp[0] not in variant_list:
			variant_list.append(tmp[0])
	for each in data[1:]:
		each = each.replace("\n", "")
		each = each.replace("\r", "")
		if len(each) == 0:
			continue
		if each[0] == "#":
			continue
		tmp = each.split("\t")
		pool[tmp[0]].append([tmp[1], int(tmp[2]), tmp[3].upper(), tmp[4].upper()])
	return pool, variant_list

def sanityCheck_info(info_pool, seq_pool):
	for eachVariant in info_pool.keys():
		site_list = info_pool[eachVariant]
		tmpTemplate = site_list[0][0]
		for e in site_list:
			if e[0] != tmpTemplate:
				print("ERROR: multiple templates are found in variant: " + eachVariant)
				return "error"
			if e[0] not in seq_pool.keys():
				print("ERROR: template " + e[0] + " of variant " + eachVariant + " are not found in input FASTA file")
				return "error"
		site_used = {}
		for e in site_list:
			position = e[1]
			refbase_length = len(e[2])
			for i in range(refbase_length):
				site_used[position + i] = []
		for e in site_list:
			position = e[1]
			refbase_length = len(e[2])
			for i in range(refbase_length):
				site_used[position + i].append(e)
		for e in site_used:
			mutation_list = site_used[e]
			if len(mutation_list) > 1:
				overlap_base = []
				for p in mutation_list:
					overlap_base.append(str(p[1]) + "_" + p[2])
				print("ERROR: overlapped mutations are found in template " + tmpTemplate + " of variant " + eachVariant + ": " + \
					";".join(overlap_base))
				return "error"

		for eachSite in site_list:
			tmpPosition = eachSite[1]
			tmpRef = eachSite[2]
			tmpInputRef = seq_pool[tmpTemplate][(tmpPosition - 1):(tmpPosition - 1 + len(tmpRef))]
			if tmpInputRef != tmpRef:
				print("ERROR: provided reference sequence of position " + str(tmpPosition) + " in template " + tmpTemplate + \
					" is not matched in input FASTA file:\nProvided in info file: " + tmpRef + "\nFound in FASTA file: " + tmpInputRef)
				return "error"
	return "pass"
def determinePrimerCount(info_pool, seq_pool, end_extension):
	variant_primer_pool = {}
	for eachVariant in info_pool.keys():
		variant_primer_pool[eachVariant] = []
	for eachVariant in info_pool.keys():
		site_list = info_pool[eachVariant]
		tmpTemplate = site_list[0][0]
		reference_fasta = seq_pool[tmpTemplate]

		site_position_list = [[e[1], e[2], e[3]] for e in site_list]
		site_position_list.sort()
		previous_start = 1
		previous_reflength = 1
		previous_primer_index = 0
		variant_primer_pool[eachVariant].append([[previous_start,
						reference_fasta[(previous_start - 1):(previous_start - 1 + previous_reflength)],
						reference_fasta[(previous_start - 1):(previous_start - 1 + previous_reflength)]],])

		for i in range(len(site_position_list)):
			tmpSite, tmpRef, tmpAlt = site_position_list[i]
			if tmpSite - (previous_start + previous_reflength) - 1 <= end_extension * 2 + 2:
				variant_primer_pool[eachVariant][previous_primer_index].append([tmpSite, tmpRef, tmpAlt])
			else:
				previous_primer_index += 1
				variant_primer_pool[eachVariant].append([[tmpSite, tmpRef, tmpAlt],])
			previous_start = tmpSite
			previous_reflength = len(tmpRef)
	return variant_primer_pool

def Tm_calculation(oligo_seq, internal_paramters):
	return mt.Tm_NN(oligo_seq, Na = internal_paramters["Na"], \
				Tris = internal_paramters["Tris"], \
				Mg = internal_paramters["Mg"])

def designPrimer(variant_primer_pool, info_pool, seq_pool, end_extension, p3_GCclamp, \
			p5_Tm_cut, p3_Tm_high, p3_Tm_low, internal_paramters):

	variant_design_pool = {}
	for eachVariant in variant_primer_pool.keys():
		site_list = info_pool[eachVariant]
		tmpTemplate = site_list[0][0]
		reference_fasta = seq_pool[tmpTemplate]

		primer_list = variant_primer_pool[eachVariant]
		primer_num = len(primer_list)
		primer_5p_length_list = []
		primer_5p_Tm_list = []
		for eachPrimer in primer_list:
			site_5p = [e[0] for e in eachPrimer]
			site_5p.sort()
			site_5p_start = site_5p[0]
			if site_5p_start == 1:
				primer_5p_length = 0
				primer_5p_length_list.append(primer_5p_length)
				primer_5p_Tm_list.append(0)
			else:
				tmp_5p_length = end_extension
				while(True):
					tmp_5p_seq = reference_fasta[(site_5p_start - 1 - tmp_5p_length):(site_5p_start - 1)]
					tmp_5p_Tm = Tm_calculation(tmp_5p_seq, internal_paramters)
					if tmp_5p_Tm >= p5_Tm_cut:
						break
					else:
						tmp_5p_length += 1
						continue
				primer_5p_length = tmp_5p_length
				primer_5p_length_list.append(primer_5p_length)
				primer_5p_Tm_list.append(tmp_5p_Tm)

		tm_p3_target = [e * (p3_Tm_high - p3_Tm_low) / (primer_num - 1) + p3_Tm_low for e in range(primer_num)]
		primer_3p_length_list = []
		primer_3p_Tm_list = []
		primer_3p_flag_list = []
		for i in range(len(primer_list)):
			site_3p = [e[0] + len(e[1]) for e in primer_list[i]]
			site_3p.sort()
			site_3p.reverse()
			site_3p_start = site_3p[0] - 1
			if i == 0:
				previous_p3_Tm = 0
			else:
				previous_p3_Tm = primer_3p_Tm_list[i - 1]

			if i == (len(primer_list) - 1):
				site_3p_length_max_pos = len(reference_fasta)
			else:
				site_5p_next = [e[0] for e in primer_list[i + 1]]
				site_5p_next.sort()
				site_5p_next_start = site_5p_next[0]
				primer_5p_length_next = primer_5p_length_list[i + 1]
				site_3p_length_max_pos = site_5p_next_start - primer_5p_length_next - 1

			if site_3p_start + end_extension > site_3p_length_max_pos:
				primer_3p_length = site_3p_length_max_pos - site_3p_start
				primer_3p_length_list.append(primer_3p_length)
				primer_3p_Tm_list.append(Tm_calculation(reference_fasta[site_3p_start:(site_3p_start + primer_3p_length)], internal_paramters))
				primer_3p_flag_list.append("Just next to the next primer")
			else:
				tmp_3p_length = end_extension
				while(True):
					tmp_3p_seq = reference_fasta[site_3p_start:(site_3p_start + tmp_3p_length)]
					tmp_3p_Tm = Tm_calculation(tmp_3p_seq, internal_paramters)
					if tmp_3p_Tm >= tm_p3_target[i] and tmp_3p_Tm >= previous_p3_Tm:
						if p3_GCclamp == 1:
							if tmp_3p_seq[-1] in ["G", "C"]:
								primer_3p_flag_list.append("")
								break
							else:
								tmp_3p_length += 1
								continue
						else:
							primer_3p_flag_list.append("")
							break
					elif site_3p_start + tmp_3p_length >= site_3p_length_max_pos:
						primer_3p_flag_list.append("Just next to the next primer")
						break
					else:
						tmp_3p_length += 1
						continue
				primer_3p_length = tmp_3p_length
				primer_3p_length_list.append(primer_3p_length)
				primer_3p_Tm_list.append(tmp_3p_Tm)


		variant_design_pool[eachVariant] = [primer_5p_length_list, primer_5p_Tm_list, \
											primer_3p_length_list, primer_3p_Tm_list, primer_3p_flag_list]
	return variant_design_pool

def generateOutput(variant_design_pool, variant_primer_pool, info_pool, seq_pool, warining_primer_length, warning_oligo_gap):
	output_primerInfo_pool = {}
	for eachVariant in variant_primer_pool.keys():
		output_primerInfo_pool[eachVariant] = []
	output_finalSeq_pool = {}
	for eachVariant in variant_primer_pool.keys():
		site_list = info_pool[eachVariant]
		tmpTemplate = site_list[0][0]
		reference_fasta = seq_pool[tmpTemplate]
		primer_list = variant_primer_pool[eachVariant]
		primer_5p_length_list, primer_5p_Tm_list, primer_3p_length_list, primer_3p_Tm_list, primer_3p_flag_list = variant_design_pool[eachVariant]
		output_finalSeq = ""
		output_finalSeq_previousPos = 0
		previous_primer_end = -999999
		
		for i in range(len(primer_list)):
			mutations = primer_list[i]
			p5_length = primer_5p_length_list[i]
			p5_Tm = primer_5p_Tm_list[i]
			p3_length = primer_3p_length_list[i]
			p3_Tm = primer_3p_Tm_list[i]
			p3_flag = primer_3p_flag_list[i]
			mutations.sort()

			ifStartPrimer = 0
			site_5p = [e[0] for e in mutations]
			site_5p.sort()
			site_5p_start = site_5p[0]

			site_3p = [e[0] + len(e[1]) for e in mutations]
			site_3p.sort()
			site_3p.reverse()
			site_3p_start = site_3p[0] - 1
	
			primer_seq = ""
			primer_seq += reference_fasta[(site_5p_start - 1 - p5_length):(site_5p_start - 1)]
			primer_start_site = (site_5p_start - p5_length)
			previous_end = site_5p_start - 1

			for eachM in mutations:
				tmpSite, refBase, altBase = eachM
				output_finalSeq += reference_fasta[output_finalSeq_previousPos:(tmpSite - 1)]
				primer_seq += reference_fasta[previous_end: (tmpSite - 1)]
				if tmpSite == 1:
					output_finalSeq += altBase
					primer_seq += altBase
					ifStartPrimer = 1
				else:
					output_finalSeq += altBase.lower()
					primer_seq += altBase.lower()
				output_finalSeq_previousPos = tmpSite + len(refBase) - 1
				previous_end = tmpSite + len(refBase) - 1

			primer_seq += reference_fasta[site_3p_start:(site_3p_start + p3_length)]

			tmp_notes = []
			if len(mutations) >= 2:
				tmp_notes.append("Contains multiple mutations")
			if ifStartPrimer == 1:
				tmp_notes.append("Extension F-primer")

			tmp_warning = []
			if len(primer_seq) >= warining_primer_length:
				tmp_warning.append("Primer longer than " + str(warining_primer_length) + "nt")
			if p3_flag != "":
				tmp_warning.append(p3_flag)
			elif (site_5p_start - p5_length) - previous_primer_end - 1 <= warning_oligo_gap:
				tmp_warning.append("Close to the upstream primer: " + str((site_5p_start - p5_length) - previous_primer_end - 1) + "nt gaps")
			previous_primer_end = (site_3p_start + p3_length)

			tmp_primer_info = [eachVariant, tmpTemplate, eachVariant + "-P" + str(i), str(primer_start_site), primer_seq, \
							str(len(primer_seq)), round(p5_Tm, 1), round(p3_Tm, 1), \
							";".join([str(e[0]) + "_" + e[1] + "_" + e[2] for e in mutations if e[0] != 1]), ";".join(tmp_notes), ";".join(tmp_warning)]

			output_primerInfo_pool[eachVariant].append(tmp_primer_info)
		output_finalSeq += reference_fasta[output_finalSeq_previousPos:]
		output_finalSeq_pool[eachVariant] = output_finalSeq
	return output_primerInfo_pool, output_finalSeq_pool
if __name__ == "__main__":
	main()
