import getopt, sys
try:
   opts, args = getopt.getopt(sys.argv[1:], 'hc:', ['help', 'coverage='])
except getopt.GetoptError as error:
   print str(error)
   sys.exit(2)

cov = 1
for opt, arg in opts:
    if opt in ('-h', '--help'):
       print 'USAGE: python methylCall_from_Bismark.py [options] <sample_name> <input_dir> <output_dir>\n'
       print 'Arguments:\n'
       print '<sample_name> Input fastq names to bismark.\n'
       print '<input_dir>   Input directory where files are stored.\n'
       print '<output_dir>  Write all output files into this directory. \n'
       print 'Options:\n'
       print '-h/--help Display the help file.\n'
       print '-c/--coverage The minimun coverage needed to call methylation for a cytosine. Default 1.\n' 
       sys.exit()
    elif opt in ('-c', '--coverage'): cov = int(arg)

if len(sys.argv)-1-2*len(opts) != 3: raise EOFError('Cannot locate input file or output directors!')
sample = sys.argv[len(sys.argv)-3]
input_dir = sys.argv[len(sys.argv)-2]
output_dir = sys.argv[len(sys.argv)-1]

STATE = 1
CHROM = 2
BASE = 3
MET_label = 4
C_read = 0
T_read = 1

sample_list = sample.split('.')
output_name = output_dir + '/' + 'cpg.' + sample_list[0] + '.mincov' + str(cov) + '.txt'
output_file = open(output_name, 'w')
line = 'chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n'
output_file.write(line)

for strand in ('OT', 'OB'):
    input_name = input_dir + '/' + 'CpG_' + strand + '_' + sample + '_trimmed.fq_bismark.txt'
    print 'Extracting information from', input_name

    label = 'F' if strand == 'OT' else 'R' 
    dict = {}
    error_reads = 0

    input_file = open(input_name)
    for line in input_file:
        line = line.strip()
        string = line.split('\t')
        if len(string) > 1:
           if not (string[CHROM] in dict): dict[string[CHROM]] = {}
           base = int(string[BASE])
           if not (base in dict[string[CHROM]]): dict[string[CHROM]][base] = [0, 0]
           if (string[STATE] == '+') and (string[MET_label] == 'Z'): dict[string[CHROM]][base][C_read] += 1
           elif (string[STATE] == '-') and (string[MET_label] == 'z'): dict[string[CHROM]][base][T_read] += 1
           error_reads += 0
    input_file.close()
    if error_reads > 0: print error_reads, 'error records found!'         
    key_chr = dict.keys()
    key_chr.sort()
    for chr in key_chr:
        key_base = dict[chr].keys()
        key_base.sort()
        for base in key_base:
            C_total = dict[chr][base][C_read]
            T_total = dict[chr][base][T_read]
            read_total = C_total + T_total
            if read_total >= cov:
               line = chr + '.' + str(base) + '\t' + chr + '\t' + str(base) + '\t' + label + '\t' + str(read_total) + '\t'
               freqC = round(float(C_total)/read_total*100, 3)
               line += (str(freqC) + '\t')
               freqT = round(float(T_total)/read_total*100, 3)
               line += (str(freqT) + '\n')
               output_file.write(line)
    print 'Methylation information have been extracted from', input_name
output_file.close()
