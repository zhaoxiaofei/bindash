import os, sys

lines = [line for line in sys.stdin if not line.startswith('#')][1:]
files = [line.strip().split('\t')[19] for line in lines if
	     (line.strip().split('\t')[4] in [
		 'reference genome', 
		 #'representative genome'
		 ])]

i = 0
for f in files:
	#if 'GCF_000691975.1_ASM69197v1' in f: break
	i += 1
#files = files[i:]
sys.stdout.write('>>> Total {} files\n'.format(len(files)))
sys.stdout.flush()
for file in files:
	sys.stdout.write('Getting the file {}\n'.format(file))
	sys.stdout.flush()
	code = os.system('wget {}/*_genomic.fna.gz -R from_genomic.fna.gz'.format(file))
	sys.stdout.write('	Result = {} for file {}\n'.format(code, file))
	sys.stdout.flush()
	
