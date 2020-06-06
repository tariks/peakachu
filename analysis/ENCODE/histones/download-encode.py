import subprocess, os

metadata = [line.rstrip().split('\t') for line in open('metadata.tsv')]
link_mode = 'https://www.encodeproject.org/files/{0}/@@download/{0}.bed.gz'

cell = 'GM12878'

label1 = ['bed broadPeak', 'bed narrowPeak']
label2 = 'peaks'
label43 = 'hg19'

cache = set()

for m in metadata[1:]:

    if not m[1] in label1:
        continue
    
    if m[2]!=label2:
        continue
    
    if m[43]!=label43:
        continue
    
    histone = m[18].replace('-human','')
    outfil = '.'.join([cell, histone, label43, 'bed.gz'])

    if outfil in cache:
        continue
    
    if os.path.exists(outfil):
        continue
    
    cache.add(outfil)
    link = link_mode.format(m[0])
    command = ['curl', '-o', outfil, '-L', link]

    try:
        subprocess.check_call(' '.join(command), shell=True)
    except:
        pass


