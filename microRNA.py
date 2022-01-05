#!/usr/bin/env python
#20190816 initial version, y:/home/mengf_lab/daimh/microRNA-standalone
import sys, argparse, subprocess, hashlib, os.path, glob, re, multiprocessing, piping

_header = [
	'miRNA_name',
	'miRNA_count_rpm',
	'miRNA_count_raw',
	'Sum_of_Aligned_Reads',
	'Read_count',
	'Reads_per_miRNA',
	'Num_mapped_miRNA_for_this_read',
	'Offset_from_miRBase_seq',
	'Most_Freq_Seq',
	'miRBase_Seq'
]
def print_list(fs, fout):
	print('\t'.join(['{:.1f}'.format(s) if isinstance(s, float) else str(s) for s in fs]), file=fout)
def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', "--fastq-directory", help="fastq files directory", required=True)
	parser.add_argument('-m', "--microRNA-fasta-file", help="microRNA fasta file", required=True)
	parser.add_argument('-o', "--output-directory", help="output directory", required=True)
	parser.add_argument("-f", "--force", help='override files under output directory', action='store_true')
	parser.add_argument("-a", "--adapter-sequence", required=False)
	parser.add_argument("-p", "--parallel", help='default is the number of cores', type=int, required=False)
	parser.add_argument("-b", "--buffer-size", help='parameter for "sort -S"', required=False)
	return parser.parse_args()
def test_env(args):
	if os.path.exists(args.output_directory) and not args.force:
		raise Exception(f'ERR-901: "{args.output_directory}" already exists')
	if args.adapter_sequence: 
		_versions = {
			'a4efab2dba4ada14ab36da141ae417be': '2.4'
		}
		try:
			ps = subprocess.run(['cutadapt', '--help'], capture_output=True)
		except FileNotFoundError:
			raise Exception('ERR-902: cutadapt is not in your PATH environment.')
		md5 = hashlib.md5(ps.stdout).hexdigest()
		if md5 not in _versions:
			raise Exception(f'ERR-903: microRNA.py supports these cutadapt versions only. {",".join(_versions.values())}')
		args.adapter_sequence = args.adapter_sequence.upper()
		if re.match('^[ACGT]+$', args.adapter_sequence) is None:
			raise Exception(f'ERR-904: adapter sequence has non acgt character')
	fqs = glob.glob(f'{args.fastq_directory}/*')
	if len(fqs) == 0:
		raise Exception(f'ERR-001: there is no files under the input fastq directory')
	fqs.sort()
	samples = []
	for fq in fqs:
		good = False
		for end in ['.fq', '.fastq', 'fq.gz', '.fastq.gz']:
			if fq.endswith(end):
				samples.append(os.path.basename(fq)[:-len(end)])
				good = True
				break
		if not good:
			raise Exception(f'ERR-002: file {fq} is ignored because it is not .fq, .fastq, or .fq.gz or .fastq.gz', file=sys.stderr)
	if len(samples) != len(set(samples)): 
			raise Exception(f'ERR-003: some files have both fastq/fq and its gz file there', file=sys.stderr)
	if not os.path.exists(args.output_directory): os.mkdir(args.output_directory)
	return fqs, samples

def get_microRNA(filename):
	fas = []
	sid, lst = None, []
	for ln in open(filename):
		if ln[-1] != '\n': raise Exception(f'ERR-800: microRNA file is corrupted at the end')
		if ln[0] == '>':
			if sid is not None: fas.append((''.join(lst), sid))	
			sid, lst = ln[1:ln.find(' ')], []
		else:
			for c in ln[:-1]:
				if c not in ['A', 'C', 'G', 'T']:
					raise Exception(f'ERR-801: microRNA file has non-ACGT sequence')
			lst.append(ln[:-1])
	if sid is not None: fas.append((''.join(lst), sid))	
	return fas	

class FuncCounter(piping.FuncBase):
	def __init__(self):
		self.cnt = 0
	def work(self, line, fout):
		self.cnt += 1
	def close(self, fout):
		if (self.cnt % 4 != 0): raise Exception('ERR-100: remainder of line number of one fastq file divided by 4 is not zero')
		print(int(self.cnt / 4), file=fout)

class FuncAligner(piping.FuncBase):
	def __init__(self, prefix, microRNA):
		self.microRNA = microRNA
		self.report = {'UNIQUE':0, 'MULTIPLE':0, 'miRBase':0, 'UNMAPPED':0}
		self.funmapped = open(prefix + '_unmapped_seq.txt', 'w')
		self.foverlong = open(prefix + '_overlong_seq.txt', 'w')
		print('REPEAT\tSEQ', file=self.funmapped)
		print('REPEAT\tSEQ', file=self.foverlong)
	def work(self, line, fout):
		fs = line[:-1].lstrip().split(' ')
		if len(fs) != 2 or line[-1] != '\n': raise Exception('ERR-001')
		self.align_to_all_microRNA(self.microRNA, fs[1], int(fs[0]), fout)
	def align_to_all_microRNA(self, microRNA, rseq, read_repeat, fout):
		if len(rseq) < 18: raise Exception('ERR-700')
		if len(rseq) > 30:
			print(read_repeat, rseq, sep='\t', file=self.foverlong)
			return 
		rev = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
		min_match_length = 17
		alns = []
		for sseq, sid in microRNA:
			aln = self.align_to_one_microRNA(min_match_length, sid, sseq, rseq)
			if aln is None: continue
			alns.append(aln)
		if len(alns) == 0: 
			print(read_repeat, rseq, sep='\t', file=self.funmapped)
			self.report['UNMAPPED'] += 1
			return
		alns.sort()
		rna_count = 1
		for idx in range(1, len(alns)):
			if alns[idx][:3] != alns[0][:3]: break
			rna_count += 1
		if rna_count == 1:
			self.report['UNIQUE'] += read_repeat
		else:
			self.report['MULTIPLE'] += read_repeat
		for idx in range(rna_count):
			aln = alns[idx] + [rna_count, rseq, read_repeat]
			aln[1] = -aln[1]
			if aln[3] >= -1 and aln[3] <= 1:
				self.report['miRBase'] += read_repeat
			print(','.join([str(i) for i in aln]), file=fout)
	def align_to_one_microRNA(self, min_match_length, subject_id, subject_seq, query_seq):
		slen, qlen = len(subject_seq), len(query_seq)
		min_rslt = None
		for idx in range(-qlen+min_match_length, slen - min_match_length + 1):
			pm, mm_all, mm19 = 0, 0, 0
			for bp in range(max(0, idx), min(slen, qlen+idx)):
				if subject_seq[bp] == query_seq[bp-idx] and subject_seq[bp] != 'N':
					pm += 1
				else:
					mm_all += 1
					if bp-idx < 19:
						mm19 += 1
						if mm19 > 1: break
			if pm < min_match_length: continue
			if mm19 > 1: continue
			cur_rslt = [mm19, 0-pm, mm_all, idx, subject_id, subject_seq]
			if min_rslt is None or min_rslt > cur_rslt:
				min_rslt = cur_rslt
		return min_rslt

class FuncSummary(piping.FuncBase):
	def __init__(self, prefix):
		self.dict_summary_one = {}
		self.prev = None
		self.lst = []
		self.mirna = {}
		self.fone = open(prefix + '_summary_one.txt', 'w')
		self.fall = open(prefix + '_summary_all.txt', 'w')
		print('\t'.join(_header), file=self.fone)
		print('\t'.join(_header), file=self.fall)
	def work(self, line, fout):
		fs = line[:-1].split(',')
		if self.prev != fs[4]:
			self.summarize()
			self.prev = fs[4]
			self.lst = []
		self.lst.append(fs)
#		print(fs, file=fout)
	def close(self, fout):
		self.summarize()
		self.printseq(self.mirna)
	def summarize(self):
		if len(self.lst) == 0: return
		sum_of_aligned_reads, miRNA_count_raw, offset_eq_1_read_count, offset_gt_1_read_count = 0, 0, 0, 0
		for fs in self.lst:
			rna_count = int(fs[6])
			read_count = float(fs[8]) / rna_count
			sum_of_aligned_reads += read_count
			offset = int(fs[3])
			if offset == 0:
				miRNA_count_raw += read_count
			elif offset in [1, -1]:
				miRNA_count_raw += read_count
				offset_eq_1_read_count += read_count
			else:
				offset_gt_1_read_count += read_count
			fs.append([int(fs[8]), read_count])
		fraction_5_not_mirbase = (offset_eq_1_read_count + offset_gt_1_read_count)/sum_of_aligned_reads
		for fs in self.lst:
			offset = int(fs[3])
			rseq, rna = fs[7], fs[5]
			if offset < 0:
				left = rseq[:-offset]
			else:
				left = ''
			diff = []
			lower_seq = []
			for i in range(max(0, offset), min(len(rseq)+offset, len(rna))):
				if rseq[i-offset] == rna[i]: 
					lower_seq.append(rna[i])
				else:
					lower_seq.append(rseq[i-offset].lower())
					diff.append('%d:%c' % (i+1, rseq[i-offset]))
			if len(diff) != int(fs[2]): raise Exception('Err3:%s' % str(fs))
			if len(rseq) + offset > len(rna):
				right = rseq[len(rna) - len(rseq) - offset:]
			else:
				right = ''
			lower_seq = ''.join(lower_seq)
			if (left + lower_seq + right).lower() != rseq.lower(): raise Exception('Err4')
			if offset > 0: lower_seq = '%s%s' % (''.join(['_' for i in range(offset)]), lower_seq)
			sequence = ''.join(lower_seq)
			if left != '': sequence = '%s<%s' % (left, sequence)
			if right != '': sequence = '%s>%s' % (sequence, right)
			if self.mirna.get('ID') != fs[4]:
				self.printseq(self.mirna)
				self.mirna['ID'] = fs[4]
				self.mirna['LST'] = []
			self.mirna['LST'].append((fs[4], 1000000 * miRNA_count_raw/ self.total, miRNA_count_raw, sum_of_aligned_reads, fs[-1][0], fs[-1][1], fs[6], offset, sequence, rna))
	def printseq(self, mirna):
		if len(mirna) == 0: return
		first = True
		for fs in sorted(mirna['LST'], key=lambda fs: fs[4], reverse=True):
			if first:
				print_list(fs, self.fone)
				self.dict_summary_one[fs[0]] = fs
				first = False
			print_list(fs, self.fall)

def worker_per_fq(args, microRNA, fq, sample):
	sort_cmd = ['sort']
	if args.buffer_size is not None: sort_cmd += ['-S', args.buffer_size]
	if fq.endswith('.gz'):
		p_fq = piping.Pipe(['gunzip', '-c'], open(fq))
	else:
		p_fq = piping.Pipe(['cat'], open(fq))
	p_counter = piping.Pipe(FuncCounter(), subprocess.PIPE)
	if args.adapter_sequence:
		p_align = piping.Pipe(['cutadapt', '-m', '18', '-a', args.adapter_sequence, '-'], subprocess.PIPE)
	else:
		p_align = piping.Pipe(['cat'], subprocess.PIPE)
	p_align.append(['sed', '-n', '2~4p'])
	p_align.append(sort_cmd)
	p_align.append(['uniq', '-c'])
	prefix = args.output_directory + '/' + sample
	f_align = FuncAligner(prefix, microRNA)
	p_align.append(f_align)
	p_align.append(sort_cmd + ['-t', ',', '-k', '5,5', '-k', '8,8'])
	f_summary = FuncSummary(prefix)
	p_align.append(f_summary)
	p_fq.tee([p_counter, p_align])
	f_summary.total = int(p_counter.stdout.read())
	p_fq.close()
	p_counter.close()
	p_align.close()
	f_align.report['TOTAL'] = f_summary.total
	return {
		'SUMMARY'	: f_summary.dict_summary_one, 
		'REPORT'	: f_align.report,
	}

def report_all_samples_column(samples, summaries, args, column):
	fout = open(args.output_directory + '/ALL_SAMPLES_' + _header[column] + '.txt', 'w')
	print('miRNA_name\t' + '\t'.join(samples), file=fout)
	for rna in sorted(set(sum([list(s.keys()) for s in summaries], []))):
		row = [rna]
		for i in range(len(samples)):
			fs = summaries[i].get(rna)
			if fs is None:
				row.append('0')
			else:
				row.append(fs[column])
		print_list(row, fout)

def report_all_samples_MFS(samples, summaries, args, column):
	fout = open(args.output_directory + '/ALL_SAMPLES_' + _header[column] + '.txt', 'w')
	print('miRNA_name\tavg_miRNA_count_rpm\tsame_as_MFS\tsame_as_miRBase\t' + '\t'.join(samples), file=fout)
	tab = []
	for rna in sorted(set(sum([list(s.keys()) for s in summaries], []))):
		row = []
		avg = 0.0
		for i in range(len(samples)):
			fs = summaries[i].get(rna)
			if fs is None:
				row.append('0')
			else:
				row.append(fs[column])
				avg += float(fs[1])
		avg /= len(samples)
		same_across_all_sample = len(set(row)) == 1
		if same_across_all_sample:
			identical_to_mirbase = row[0] == fs[9]
		else:
			identical_to_mirbase = False
		same_across_all_sample = str(same_across_all_sample)[0]
		identical_to_mirbase = str(identical_to_mirbase)[0]
		tab.append((avg, '%s\t%.1f\t%s\t%s\t%s' % (rna, avg, same_across_all_sample, identical_to_mirbase, '\t'.join(row))))
	for avg, ln in sorted(tab, reverse=True):
		print(ln, file=fout)

def report_all_samples_count(samples, reports, args):
	fout = open(args.output_directory + '/ALL_SAMPLES_REPORT.txt', 'w')
	print('Sample\tTotal\tCount_of_reads_mapped_to_a_miRBase_annotation\tmapped_to_one_annotation_only\tmapped_to_multiple_annotations\tmiRNA_count\tmiRNA_percentage\tunmapped_read_count', file=fout)
	for i in range(len(samples)):
		rpt = reports[i]
		row = (samples[i], rpt['TOTAL'], rpt['UNIQUE'] + rpt['MULTIPLE'], rpt['UNIQUE'], rpt['MULTIPLE'], rpt['miRBase'], 100.0 * rpt['miRBase'] / rpt['TOTAL'], rpt['UNMAPPED'])
		print_list(row, fout)

def main():
	args = get_args()
	if args.parallel is None: args.parallel = multiprocessing.cpu_count()
	fqs, samples = test_env(args)
	microRNA = get_microRNA(args.microRNA_fasta_file)
	with multiprocessing.Pool(args.parallel) as pool:
		results = pool.starmap(worker_per_fq, [(args, microRNA, fqs[i], samples[i]) for i in range(len(fqs))])
	summaries = [t['SUMMARY'] for t in results]
	report_all_samples_column(samples, summaries, args, 1)
	report_all_samples_column(samples, summaries, args, 2)
	report_all_samples_MFS(samples, summaries, args, 8)
	reports = [t['REPORT'] for t in results]
	report_all_samples_count(samples, reports, args)

if __name__ == '__main__': main()
