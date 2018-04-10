import os
import subprocess
import sys

factor = int(sys.argv[1])
fnames = sys.stdin.readlines()
fnamelist = []
for i, fname in enumerate(fnames):
	fname = fname.strip()
	fnamelist.append(fname);
	if i % factor == (factor - 1) or i == len(fnames) - 1:
		outfname = fnamelist[0].replace('.fna.gz', '_combined.fna.gz');
		assert(outfname != fnamelist[0])
		cmd = 'cat {} > {}'.format(' '.join(fnamelist), outfname)
		#print(cmd)
		os.system(cmd)
		sys.stdout.write(outfname)
		for fname in fnamelist:
			x = subprocess.check_output(['zgrep', '-c', "^>", fname])
			sys.stdout.write('\t{}\t{}'.format(fname, x.strip()))
		sys.stdout.write('\n')
		fnamelist = []

