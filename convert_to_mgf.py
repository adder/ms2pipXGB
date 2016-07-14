import sys

PTMmap = {}
PTMmap['Oxidation'] = '15.994915'
PTMmap['Acetyl'] = '42.010565'
PTMmap['Carbamidomethyl'] = '57.02146'
PTMmap['Pyro-carbamidomethyl'] = ' 39.994915'
PTMmap['Glu->pyro-Glu'] = '-18.010565'
PTMmap['Gln->pyro-Glu'] = '-17.026549'

fpip = open(sys.argv[1]+'.PEPREC','w')
fpip.write("spec_id modifications peptide charge\n")
fmgf = open(sys.argv[1]+'.PEPREC.mgf','w')
fmeta = open(sys.argv[1]+'.PEPREC.meta','w')

specid = 1
with open(sys.argv[1]) as f:
	peptide = None
	charge = None
	parentmz = None
	mods = None
	purity = None
	HCDenergy = None
	read_spec = False
	mgf = ""
	prev = 'A'
	sys.stderr.write(prev)
	for row in f:
		if read_spec:
			l = row.rstrip().split('\t')
			if len(l) != 3:
				if peptide[0] != prev:
					prev = peptide[0]
					sys.stderr.write(prev)

				tmp = mods.split('/')
				
				if tmp[0] != '0':
					specid += 1
					read_spec = False
					mgf = ""
					continue

				if (charge != '2') & (charge != 3):
					specid += 1
					read_spec = False
					mgf = ""
					continue
													
				if tmp[0] != '0':
					m = ""
					for i in range(1,len(tmp)):
						tmp2=tmp[i].split(',')
						if not tmp2[2] in PTMmap:
							print tmp2[2]
							dd
						m += tmp2[0] + '|' + PTMmap[tmp2[2]] + '|'
					fpip.write('%s%i %s %s %s\n'%(sys.argv[2],specid,m[:-1],peptide.replace('(O)',''),charge))
				else:
					fpip.write('%s%i  %s %s\n'%(sys.argv[2],specid,peptide.replace('(O)',''),charge))

				fmeta.write('%s%i %s %s %s %s %s\n'%(sys.argv[2],specid,charge,peptide,parentmz,purity,HCDenergy))
					
				buf = "BEGIN IONS\n"
				buf += "TITLE="+sys.argv[2]+str(specid) + '\n'
				buf += "CHARGE="+str(charge) + '\n'
				buf += "PEPMASS="+parentmz + '\n'
				fmgf.write(buf+mgf+"END IONS"+'\n')
				
				specid += 1
				read_spec = False
				mgf = ""
				continue
			else:
				tt = float(l[1])
				mgf += ' '.join([l[0],l[1]]) + '\n'
				continue
		if row.startswith("Name:"):
			l = row.rstrip().split(' ')
			tmp = l[1].split('/')
			peptide = tmp[0]
			charge = tmp[1].split('_')[0]
			continue
		if row.startswith("Comment:"):
			l = row.rstrip().split(' ')
			for i in range(1,len(l)):
				if l[i].startswith("Mods="):
					tmp = l[i].split('=')
					mods = tmp[1]
				if l[i].startswith("Parent="):
					tmp = l[i].split('=')
					parentmz = tmp[1]
				if l[i].startswith("Purity="):
					tmp = l[i].split('=')
					purity = tmp[1]
				if l[i].startswith("HCD="):
					tmp = l[i].split('=')
					HCDenergy = tmp[1].replace('eV','')
			continue
		if row.startswith("Num peaks:"):
			read_spec = True
			continue

fmgf.close()
fpip.close()
fmeta.close()
