# -*- coding: utf-8 -*-

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Restriction import *
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import numpy as np
import os
import pandas as pd
import re
import requests
import json

# Global variables

# Benchling API key (if not available send a request)
SECRET_KEY = 'SECRET KEY HERE'
SUBMITTOR = 'Christian Choe'
FOLDER = 'Automated Plasmid Design'
COLOR = [255, 255, 0]

# Enzyme list to ignore
ignore_enzyme_list = ['AsuNHI']


def load_m_codon(filename = 'Human_Codon_Frequency.xlsx'):
	""" 	Loads mammalian DNA codon chart 
	
	Args:
		filename: excel file containing the codon information

	Returns:
		aa2codon: dictionary with amino acid (key) to [codon, frequency] (value)
		codon2aa: dictionary with codon (key) to amino acid (value)
	"""
	# Load data and organize
	cwd = os.getcwd()
	cwd = np.array(cwd.split('/'))
	idx = np.where(cwd == 'OneDrive')[0][0]
	f_dir = '/'.join(cwd[:idx+1]) + '/Documents/'
	df = pd.read_excel(f_dir + filename)
	df_1 = df.iloc[:,0:5]
	df_2 = df.iloc[:,5:]
	df_2.columns = df_1.columns
	df = pd.concat([df_1, df_2])
	df.iloc[:,0:5]
	df = df.sort_values(['Amino acid', 'Fraction'], ascending=[True, False])
	df = df.reset_index(drop=True)
	
	# Create 2 dictionaries (aa -> codon & codon -> aa)
	aa2codon = {}
	codon2aa = {}
	for row in df.values:
		temp_codon = row[0]
		temp_aa = row[1]
		frac = row[2]
		codon2aa[temp_codon] = temp_aa
		if temp_aa in aa2codon:
			aa2codon[temp_aa][0].append(temp_codon)
			aa2codon[temp_aa][1].append(frac)
		else:
			aa2codon[temp_aa] = [[temp_codon], [frac]]
	
	# Make sure all probabilities sum to 1
	for key in aa2codon:
		aa2codon[key][1] = np.array(aa2codon[key][1])/sum(aa2codon[key][1])
	return aa2codon, codon2aa

def get_Tm(seq):
	"""	Calculate melting temperature for DNA oligomer

	Melting tempearture is based on Primestar PCR Premix condition
	
	Args:
		seq: string of DNA sequence

	Returns:
		Tm: melting temperature (Celcius)
	"""
	myseq = Seq(seq)

	# note for line below: all mM except dnacl which uses nM
	Tm = mt.Tm_NN(myseq, dnac1=250, Na=0, K=10, Tris=100, Mg=2, dNTPs=0.8) 
	return Tm

def get_rc(seq):
	"""	Returns the reverse complement sequence """
	return str(Seq(seq).reverse_complement())

def get_translation(seq):
	""" Returns the translation of a coding dna or rna sequence """
	assert(len(seq)%3 == 0), 'Sequence length is not divisible by 3!'
	seq = seq.upper()
	seq = seq.replace('U','T')
	aa2codon, codon2aa = load_m_codon()
	aa_seq = ''
	for i in range(len(seq)/3):
		aa_seq += codon2aa[seq[i*3:i*3+3]] 
	return aa_seq

def clean_seq(seq):
	"""	cleans up a sequence by removing all non alphabetic characters """
	pattern = re.compile('[\W_]+', re.UNICODE)
	return pattern.sub('', seq)

def check_GC(seq):
	"""	Measures the GC content of a sequence """
	return (seq.upper().count('G') + seq.upper().count('C'))/float(len(seq))

def check_anneal(seq, target_Tm=60, delta_Tm = 1.5, GC_bound=[0.4, 0.6], GC_3prime=True, verbose=False):
	""" Checks a primer sequence for proper annealing

	Checks all the following
	1. Target melting temperature (usually >60 C and <80 C)
	2. GC content is within bound (usualy >40% and <60%)
	3. 3 prime end contains a G or C

	Args:
		seq: string of DNA sequence
		target_Tm: target melting temperature in Celcius
		delta_Tm: margin of error for target_Tm in Celcius
		GC_bound: boundary for acceptable GC content [low, high]
		GC_3prime: boolean indicating 3 prime end contains a G or C
		verbose: boolean

	Returns:
		Boolean for whether or not the sequence is a good primer
	"""
	fail = []
	if seq[-1] not in ['G', 'C'] and GC_3prime:
		fail.append('Missing 3\' G/C')
	if seq[-3:] in ['GGG', 'CCC']:
		fail.append('Too many 3\' G/C')
	if check_GC(seq) < GC_bound[0] or check_GC(seq) > GC_bound[1]:
		fail.append('GC content out of bound')
	if np.abs(get_Tm(seq) - target_Tm) > delta_Tm:
		fail.append('Tm out of bound')
	if len(fail) == 0:
		return True		
	if verbose:
		print seq
		print 'Tm: {}'.format(get_Tm(seq))
		print 'GC content: {}'.format(check_GC(seq))
		print fail
	return False
	
def check_codon(seq, GC_bound=[0.4, 0.6]):
	""" Checks GC content of a coding sequence

	Args:
		seq: string of DNA sequence
		GC_bound: boundary for acceptable GC content [low, high]

	Returns:
		Boolean for GC content check
	"""

	GC_content = check_GC(seq)
	if GC_content < GC_bound[0] or GC_content > GC_bound[1]:
		return False
	return True

def choose_aa_codon_bound(aa_list, mask, take_first=False):
	""" Picks the codon for an amino acid that matches the mask

	By using the mask, the codon is forced to contain the the same sequence.
	E.g. mask = 'CC*'; Only codons containing 'CC' as the first two bases 
	will be returned

	Args:
		aa_list: string or list of acceptable amino acid codons
			e.g. 'SPG' or ['S', 'P', 'G'] (serine, proline, glycine)
		mask: 3 character string where blanks are represented as *
			e.g. 'C*A' allows for CAA, CCA, CGA, CTA if these codons
			code for an amino acid in aa_list
		take_first: boolean to take the first match possible
			e.g. If you listed the amino acid from most to least desirable

	Returns:
		Random codon for an amino acid in aa_list that fits the mask
	"""

	assert(len(mask)==3), 'Wrong length for codon mask'
	mask = mask.replace('*', '\w')
	p = re.compile(mask, re.IGNORECASE)
	codon_option = []
	aa2codon, codon2aa = load_m_codon()
	for ch in aa_list:
		codon_list = list(aa2codon[ch][0])
		codon_p = list(aa2codon[ch][1])
		for i in range(len(codon_list))[::-1]:
			if not p.match(codon_list[i]):
				codon_list.pop(i)
				codon_p.pop(i)
		if len(codon_list) == 0: continue	  
		random_codon = np.random.choice(codon_list, 1, p=np.array(codon_p)/sum(codon_p))[0]
		if take_first:
			return rando_codon
		codon_option.append(random_codon)
	assert(len(codon_option)>0), 'Expand amino acid choices. None fit the required mask'
	return str(np.random.choice(codon_option))
		

def generate_aa_codon(seq, N=1000, GC_bound=[0.4, 0.6], start='', end=''):
	""" Returns a coding DNA sequence for the amino acid sequence

	Args:
		seq: string of amino acid sequence
		N: number of iteration to try creating a DNA sequence
		GC_bound: boundary for acceptable GC content [low, high]
		start: DNA sequece to be appended to the start
		end: DNA sequence to be appended to the end

	Returns:
		DNA sequence for the corresponding amino acid sequence input
	"""

	seq = [x.upper() for x in seq]
	assert(len(start)<=3 and len(end)<=3), 'Please remove determined codons from requirement'
	count = 0
	aa2codon, codon2aa = load_m_codon()
	seq_list = []
	while True:
		temp_seq = ''
		for i, ch in enumerate(seq):
			if start and i==0:
				temp_seq += choose_aa_codon_bound(ch, start)
				continue
			elif end and i==len(seq):
				temp_seq += choose_aa_codon_bound(ch, end)
				continue
			if len(ch) > 1:
				ch = np.random.choice(list(ch))
			temp_seq += np.random.choice(aa2codon[ch][0], 1,p = aa2codon[ch][1])[0]
		temp_seq = str(temp_seq)
		NN_repeat, repeat = check_repeat(temp_seq)
		GC_ratio = check_GC(temp_seq)
		GC_content = (GC_ratio < GC_bound[0] or GC_ratio > GC_bound[1])
		if (NN_repeat and not repeat and GC_content) or count >= N:
			break
		GC_dif = min(np.abs(GC_ratio-np.array(GC_bound)))*(not GC_content)
		seq_list.append([temp_seq, NN_repeat, not repeat, GC_content, GC_dif])
		count += 1
	if count >= N:
		# choose best from seq_list
		seq_list = np.array(seq_list, dtype=object)
		seq_list = seq_list[seq_list[:,-1].argsort()]
		check_list = np.prod(seq_list[:,1:3], axis=1)
		temp_seq = seq_list[np.where(check_list == 1)[0][0],0]
	return temp_seq

def check_repeat(seq):
	""" Check for repeats in the DNA sequence

	Args:
		seq: string of DNA sequence

	Returns:
		two boolean corresponding to the following
			dinucletide repeat: True if more than 4 repeat
				e.g. TATATATA
			single nucleotide repeat: True if more than 4 repeat
				e.g. AAAAA
	"""

	seq = seq.upper()
	max_count = 1
	
	# Check dinucleotide repeat
	for i in range(len(seq)-1):
		substring = seq[i:i+2]
		count = 1
		tempstring = seq[i+2:]
		while len(tempstring) >= 2:
			if tempstring[:2] == substring:
				count += 1
				max_count = max(count, max_count)
				tempstring = tempstring[2:]
			else:
				break
	
	# Check single nucleotide repeat
	if 'GGGG' not in seq and 'CCCC' not in seq and 'AAAAA' not in seq and 'TTTT' not in seq:
		repeat = False
	else:
		repeat = True
	
	return max_count<4, repeat


# Checks connection to API
def get_connect():
	""" Test API connection
	
	Returns:
		Request. If the status code is fine then you're good to go	
	"""
	url = 'https://benchling.com/api/v2/plates'
	r = requests.get(url, auth=(SECRET_KEY, ''))
	assert(r.status_code == 200), 'Error status code {}'.format(r.status_code)
	return r

def get_folders(pageSize=100, projectId=None, parentFolderId=None,
				archiveReason='NOT_ARCHIVED', onlyRoot=False):
	""" Get list of all folders
	
	Args:
		pageSize: number of results per request (max is 100)
		parentFolderId: will only return children folder of a given parent
		archiveReason: if empty string returns all results
			'NOT_ARCHIVED' ignores all archived folders
		onlyRoot: boolean to return only the root folders

	Retuns:
		folder list
	"""
	if onlyRoot:
		parentFolderId = 'NO_PARENT'
	params = {'pageSize':pageSize, 'projectId':projectId, 'parentFolderId':parentFolderId, 
			  'archiveReason':archiveReason}
	url = 'https://benchling.com/api/v2/folders' 
	r_total = []
	nextToken = 'temp'
	while nextToken: # Keep requesting until no more results
		r = requests.get(url, auth=(SECRET_KEY, ''), params=params)
		r_total += r.json()['folders']
		nextToken = r.json()['nextToken']
		params['nextToken'] = nextToken
	return r_total

def get_folder_id(foldername):
	""" Returns the id of a folder """
	folder = [x for x in get_folders() if x['name'] == foldername][0]
	return folder['id']

def get_project_id(projectFolder):
	""" Returns the projectId for a given folder

	Args:
		projectFolder: folder in the project to be identified (usually root)

	Returns:
		projectId
	"""
	folder2id = {}
	id2folder = {}
	for temp in get_folders():
		folder2id[temp['name']] = {'id':temp['id'], 'projectId':temp['projectId']}
		id2folder[temp['id']] = temp['name']
	return folder2id[projectFolder]['projectId']

def get_folder2dna2id(projectId=None):
	""" For a project returns a dict for folder, dna, and id

	Args:
		projectId: project to search

	Returns
		folder2dna: dict[folder] = [dna_filenames]
		dna2id: dict[dna_filename] = id
	"""
	dna2id = {}
	folder2dna = {}

	# Load all the folder and their respective id
	id2folder = {}
	for temp in get_folders():
		id2folder[temp['id']] = temp['name']

	for temp in get_sequences(projectId=projectId):
		dna2id[temp['name']] = temp['id']
		folder_name = id2folder[temp['folderId']]
		if folder_name not in folder2dna:
			folder2dna[folder_name] = [temp['name']]
		else:
			folder2dna[folder_name].append(temp['name'])
	return folder2dna, dna2id

def get_sequences(pageSize=100, sort=None, modifiedAt=None, name=None, projectId=None, 
				  registryId=None, schemaId=None, archiveReason='NOT_ARCHIVED'):
	""" Get list of all folders
	
	Args:
		pageSize: number of results per request (max is 100)
		sort: how to sort results
			e.g. 'name'
			'' will sort by modified time
		modifiedAt: Datetime, in RFC 3339 format.
			Supports the > operator. Time zone defaults to UTC. Restricts results 
			to those modified in the specified range. e.g. > 2017-04-30
		name: name of DNA sequence
		projectId: will only return sequences of a given project
		registryId: will only return sequences of a given registry
		schemaId: will only return sequences of a given schema
		archiveReason: if empty string returns all results
			'NOT_ARCHIVED' ignores all archived folders
		
	Retuns:
		List of DNA sequences
	"""

	params = {'pageSize':pageSize, 'sort':sort, 'modifiedAt':modifiedAt, 'name':name,
			  'projectId':projectId, 'registryId':registryId, 'schemaId':schemaId,  
			  'archiveReason':archiveReason}
	url = 'https://benchling.com/api/v2/dna-sequences'
	r_total = []
	nextToken = 'temp'
	while nextToken:
		r = requests.get(url, auth=(SECRET_KEY, ''), params=params)
		r_total += r.json()['dnaSequences']
		nextToken = r.json()['nextToken']
		params['nextToken'] = nextToken
	return r_total	

def get_sequence(dna_sequence_id, dna_id=True):
	""" Reads the dna sequence of a given id

	Args:
		dna_sequence_id: the id corresponding to a dna sequence file

	Returns:
		json containing the annotated DNA sequence
	"""
	if not dna_id:
		folder2dna, dna2id = get_folder2dna2id(projectId=None)
		dna_sequence_id = dna2id[dna_sequence_id]
	url = 'https://benchling.com/api/v2/dna-sequences/{}'.format(dna_sequence_id)
	r = requests.get(url, auth=(SECRET_KEY, ''))
	return r.json()

def get_cutsite(dnaseq, label_remove, N=10, delta=0, verbose=False):
	""" Finds all cut sites near the region to remove

	Finds cut sites on both edges of label to remove. It can also search a little
	into the region to remove with delta, otherwise it will only look outside.
	If you want to insert your gene into a single existing cut site then you
	should just specify the restriction enzyme

	Args:
		dnaseq: json containing annotated DNA sequence
		label_remove: the annotated label to remove with restriction enzymes
		N: Number of results to show
			Results are ordered from nearest to farthest from the labeled region
		delta: number of bases to look inside the region
		verbose: boolean

	returns:
		enzyme1, enzyme2, direction of label
	"""
	region = [x for x in dnaseq['annotations'] if x['name'] == label_remove][0]
	direction = region['strand']
	seq = Seq(dnaseq['bases'])
	Ana = Analysis(CommOnly, seq, linear= not dnaseq['isCircular'])
	
	# Get all single cutters
	enzyme_list = []
	for key, value in Ana.with_sites().iteritems():
		if len(value) == 1:
			enzyme_list.append([key, value[0]-1]) # -1 because index starts at 1 for restriction enzymes

	# Remove enzymes from ignore_list
	for i in range(len(enzyme_list))[::-1]:
		if str(enzyme_list[i][0]) in ignore_enzyme_list:
			removed = enzyme_list.pop(i)
			if verbose: print 'Removed:\n', removed
	
	enzyme_list = np.array(enzyme_list)
	enzyme_list = enzyme_list[enzyme_list[:,1].argsort()]
	list1 = enzyme_list[np.where(enzyme_list[:,1] <= region['start']+delta)][::-1]
	list2 = enzyme_list[np.where(enzyme_list[:,1] >= region['end']-delta)][::-1]
	enzyme_list = np.concatenate((list1,list2))

	# Display top N results
	if verbose:
		for pos, d in zip(['start', 'end'],[1,-1]):
			print '\nCut {}: {}'.format(pos, region[pos])
			for x in enzyme_list[::d][:N]:
				dif = int(x[1]) - region[pos]
				loc = 'outside'
				if x[1] > region[pos] and np.abs(dif) <= delta and d == 1: 
					loc = 'inside'
				elif x[1] < region[pos] and np.abs(dif) <= delta and d == -1:
					loc = 'inside'
				loc += ':{0:+}'.format(dif)
				print '{} : {} : {}'.format(x[0], x[1], loc)
	return enzyme_list[0], enzyme_list[-1], direction

def check_cutsite(dna_insert, enzyme_list, label=''):
	""" Checks if any enzyme in the list will cut in the label region
	
	Args:
		dna_insert: dna sequence to check
		enzyme_list: list of enzymes to check
		label: region to check (if blank checks the whole dna_insert)

	Returns:
		Boolean on whether a cut was found
	"""
	if label:
		region = [x for x in dna_insert['annotations'] if x['name'] == label][0]
		seq = dna_insert['bases'][region['start']:region['end']]
	else:
		seq = dna_insert['bases']
	seq = Seq(seq)
	Ana = Analysis(CommOnly, seq, linear=not dna_insert['isCircular']) # Search over all commercially available restricton enzymes
	for key, value in Ana.with_name(enzyme_list).iteritems():
		if len(value)>0: return True
	return False

# Enzyme_list order matters. [start, end]
def get_Gibson_overlap(dnaseq, enzyme_list, direction, preserve_cut=True):
	""" Returns the overlapping dna sequence for Gibson assembly

	From the enzyme list it will grab 60 nts from the start and end for the 
	Gibson overlap
	
	Args:
		dnaseq: dna sequence to check
		enzyme_list: restriction enzymes to cut [start, end] (order matters)
		direction: direction of the inserted region
		preserve_cut: whether or not to preserve the restriction enzyme site

	Returns:
		overhang_list: list of the left and right gibson overhang
		linker_list: extra bases to preserve the restriction site
	"""

	seq = Seq(dnaseq['bases'])
	Ana = Analysis(CommOnly, seq, linear=False)
	
	overhang_list = []
	linker_list = []
	for enzyme, d in zip(enzyme_list, [1,0]):
		exec('cut_template = {}.elucidate()'.format(enzyme))
		# Location of cut
		value = Ana.with_name([enzyme]).values()[0]
		
		# Length of cut to preserve cut_site
		cut_site = cut_template.replace('_', '').find('^')
		cut_site = len(clean_seq(cut_template)) - cut_site
		
		# Correct which side to cut
		assert(len(value)==1), 'Enzyme is not a single cutter'
		if d==0:
			value[0] += cut_template.find('_') - cut_template.find('^') - 1
		overhang = seq[value[0]-1-60*d:value[0]-1+60*(not d)] # -1 due to offset index
		linker = seq[value[0]-1-60*(not d):value[0]-1+60*d] # -1 due to offset index
		linker = str(linker)
		linker = [linker[:cut_site*(2*d-1)],linker[cut_site*(2*d-1):]][::-direction]
		if d == 0:
			overhang = overhang.reverse_complement()
			linker = [get_rc(x) for x in linker]
		overhang_list.append(str(overhang))
		linker_list.append(linker)
	return overhang_list, linker_list

def get_anneal_insert(dna_insert, label_insert, direction,
					insert5p='', insert3p='', target_Tm=60, delta_Tm=1.5, 
					GC_bound=[0.4, 0.6], GC_3prime=True, linker_protein=True, 
					verbose=False):
	""" Returns the annealing sequence (PCR region)
	
	Args:
		dnaseq: dna sequence to check
		label_insert: region to anneal (to PCR)
		direction: direction of label_insert
		insert5p: bases or amino acids to add to the 5' end
		insert3p: bases or amino acids to add to the 3' end
		target_Tm: target melting tempearture for PCR in Celcius
		delta_Tm: margin of error from the target melting temperature
		GC_bound: boundary for acceptable GC content [low, high]
		GC_3prime: boolean indicating 3 prime end contains a G or C
		linker_protein: boolean on whether the insert5p and 3p are amino acids or bases
		verbose: boolean
		
	Returns:
		anneal_list: sequence for the PCR annealing
		insert_list: sequence for for the 5' and 3' insert
		boolean on whether or not it worked
	"""

	aa2codon, codon2aa = load_m_codon()
	region = [x for x in dna_insert['annotations'] if x['name'] == label_insert][0]
	idx_list = [region['start'], region['end']]
	gene_seq = dna_insert['bases'][idx_list[0]:idx_list[1]].upper()
	insert_list = ['', ''] # 5' insert, 3' insert (both are written 5' -> 3')
	anneal_list = ['', ''] # annealing region for 5' end and 3' end (both are written 5' -> 3')

	# Generate linker region
	if linker_protein:
		if insert5p:
			insert_list[0] += generate_aa_codon(insert5p)
		if insert3p:
			insert_list[1] += get_rc(generate_aa_codon(insert3p))
	else:
		insert_list[0] += insert5p
		insert_list[1] += insert3p # Make sure it's reverse complement of what you want


	# Generate annealing region
	if target_Tm < 58:
		print 'Warning: Target melting tempearture is low. Ideally 60-80 °C'
	elif target_Tm > 82:
		print 'Warning: Target melting temperature is high. Ideally 60-80 °C'

	for pos in range(2):
		temp_seq = ''
		for i in range(1,61):
			if pos:
				temp_seq = get_rc(gene_seq[-i:])
			else:
				temp_seq = gene_seq[:i]
			if check_anneal(temp_seq, target_Tm=target_Tm, delta_Tm=delta_Tm, GC_bound=GC_bound, GC_3prime=GC_3prime, verbose=verbose):
				if verbose:
					print 'Annealing GC Content: ', check_GC(temp_seq)
					print 'Annealing Tm: ', get_Tm(temp_seq)
					print 'Length: ', len(''.join(temp_seq))
					print ''.join(temp_seq) + '\n'
				break
			if i == 60:
				print 'Manual design of annealing region is required'
				return anneal_list, insert_list, False
		anneal_list[pos] += temp_seq
	return anneal_list, insert_list, True

def get_primer(overhang_list, linker_list, insert_list, anneal_list, direction, gibson_overhang=20):
	""" Combines all the lists to construct the full primer

	Args:
		overhang_list: list of the left and right gibson overhang
		linker_list: extra bases to preserve the restriction site
		insert_list: sequence for for the 5' and 3' insert (for adding more amino acids, etc)
		anneal_list: sequence for the PCR annealing
		direction: direction to insert
		gibson_overhang: length of desired Gibson overhang (20 is ideal)
			if the primer is too long it will shorten it but it will notify you if < 15 nts

	Returns
		primer_list containing the forward and reverse primer both written as 5' to 3'

	"""
	
	primer_list = []
	for i in range(2):
		primer = [overhang_list[::direction][i], linker_list[::direction][i][(i+1)%2], insert_list[i], anneal_list[i]]
		primer = [primer[0].upper(), primer[1].lower(), primer[2].lower(), primer[3].upper()]# insert and linker is lowercase
		temp_primer = ''.join(primer[1:])
		temp_length = len(temp_primer)

		# Correct for gibson overhang
		if (temp_length + gibson_overhang) > 60:
			gibson_overhang = 60 - temp_length
			if gibson_overhang < 15:
				print 'Warning. Since primer is too long Gibson overlap is now {}.'.format(gibson_overhang)
		
		overlap = primer[0][-gibson_overhang:]
		temp_primer = ''.join([overlap, temp_primer])
		primer_list.append(temp_primer)
	return primer_list

def primer_design(plasmid_base, label_remove, insert, label_insert, delta=3,
				  insert5p='', insert3p='', target_Tm=60, delta_Tm=0.75, 
				  GC_bound=[0.4,0.6], GC_3prime=True,linker_protein=True, 
				  verbose=False, new_plasmid_name='New', create_plasmid=False):

	# Step 0. Load the dna files and their id
	folder2dna, dna2id = get_folder2dna2id(projectId=None)

	# Step 1. Find restriction enzymes

	# Remove DNA for insert
	dnaseq = get_sequence(dna2id[plasmid_base])
	start_cut, end_cut, direction = get_cutsite(dnaseq, label_remove, delta=delta, verbose=verbose)

	# Step 2. Pick out 2 enzymes for cutting the 5' and 3' end
	final_enzyme = [start_cut[0], end_cut[0]]
	for enzyme in final_enzyme:
		print 'Cut: {}'.format(enzyme)

	# Double check restriction site does not exist in the DNA to be inserted
	# dna_insert = get_sequence(dna2id[insert])
	# if check_cutsite(dna_insert, enzyme_list):
	#	 print 'Pick different enzymes'

	# Step 3. Get overlap sequence for Gibson
	overhang_list, linker_list = get_Gibson_overlap(dnaseq, final_enzyme, direction)


	# Step 4. Design annealing region for primer on the insert
	dna_insert = get_sequence(dna2id[insert])
	status = False # Try to search more Tms to find anneal
	anneal_list, insert_list, status = get_anneal_insert(dna_insert, label_insert, direction, 
						insert5p=insert5p, insert3p=insert3p, target_Tm=target_Tm, 
						delta_Tm=delta_Tm, GC_bound=GC_bound, GC_3prime=True, 
						linker_protein=linker_protein, verbose=verbose)

	# Step 5. Put together the primer
	primer_list = get_primer(overhang_list, linker_list, insert_list, anneal_list, direction)
	print 'Forward primer: {} ({} nts)'.format(primer_list[0], len(primer_list[0]))
	print 'Reverse primer: {} ({} nts)'.format(primer_list[1], len(primer_list[1]))

	if create_plasmid:
		return create_plasmid_file(primer_list, direction, plasmid_base, label_remove, insert, label_insert, final_enzyme, new_plasmid_name)
	return primer_list

def create_plasmid_file(primer_list, direction, plasmid_base, label_remove, insert, label_insert, enzyme_list, new_plasmid_name):
	dnaseq = get_sequence(plasmid_base, dna_id=False)
	dna_insert = get_sequence(insert, dna_id=False)
	color = '#%02X%02X%02X' % tuple(COLOR)
	seq = dnaseq['bases']
	Ana = Analysis(CommOnly, Seq(seq), linear=False)

	# Extract data from plasmid_base
	cut_idx = []
	base_annotations = dnaseq['annotations']
	base_translations = dnaseq['translations']
	for i, enzyme in enumerate(enzyme_list):
		exec('cut_template = {}.elucidate()'.format(enzyme))
		# Location of cut
		value = Ana.with_name([enzyme]).values()[0]
		
		assert(len(value)==1), 'enzyme is not cutting once'
		
		# Length of cut to preserve cut_site
		cut_site = cut_template.replace('_', '').find('^')
		cut_site = len(clean_seq(cut_template)) - cut_site
		
		# Correct for 2nd cut site (since you want to preserve the 3' side of the cut)
		if i == 1:
			value[0] += cut_template.find('_') - cut_template.find('^') - 1
		cut_idx.append(value[0] - 1)
	remove_len = cut_idx[1] - cut_idx[0]

	# Extract data from insert
	insert_annotations = dna_insert['annotations']
	insert_translations = dna_insert['translations']
	region = [x for x in dna_insert['annotations'] if x['name'] == label_insert][0]
	idx_list = [region['start'], region['end']] # start and end of insert
	gene_seq = dna_insert['bases'][idx_list[0]:idx_list[1]].upper()
	gene_len = len(gene_seq)
	if direction == -1: gene_seq = get_rc(gene_seq)

	# Combine sequence
	primer_list = primer_list[::direction]
	primer_list = [primer_list[0], get_rc(primer_list[1])] # Flip the 2nd primer
	insert_left = re.sub('[^a-z]', '', primer_list[0])
	insert_right = re.sub('[^a-z]', '', primer_list[1])
	insert_len = len(insert_left) + len(insert_right)
	insert_left_len = len(insert_left)
	new_seq = seq[:cut_idx[0]] + insert_left + gene_seq + insert_right + seq[cut_idx[1]:]

	# Correct annotation in base
	new_annotations = []
	for annot in base_annotations:
		for x in ['start', 'end']:
			skip = False
			if cut_idx[0] <= annot[x] < cut_idx[1]: # start/end is between the cut site
				seq_annot = seq[annot['start']:annot['end']]
				regex_idx = re.search(seq_annot, new_seq, re.IGNORECASE)
				if regex_idx is not None:
					idx_start = regex_idx.start()
				else:
					skip = True
					break
				idx_end = idx_start + len(seq_annot)
				annot['start'] = idx_start
				annot['end'] = idx_end
				break
			if annot[x] > cut_idx[1]: #start/end is after the cut site
				annot[x] += gene_len + insert_len - remove_len
		if (annot['start'] != annot['end']) and (annot['name'] not in [label_remove, label_insert]) and not skip:
			new_annotations.append(annot)
			
	# Correct annotation in insert
	for annot in insert_annotations:
		if idx_list[0] <= annot['start'] <= idx_list[1]:
			if idx_list[0] <= annot['end'] <= idx_list[1]:
				if direction == -1:
					start_new = gene_len - annot['end'] + idx_list[0] 
					end_new = gene_len + idx_list[0] - annot['start']
				else:
					start_new = annot['start'] - idx_list[0]
					end_new = annot['end'] - idx_list[0]
				annot['start'] = start_new + cut_idx[0] + insert_left_len
				annot['end'] = end_new + cut_idx[0] + insert_left_len
				annot['strand'] = direction
				new_annotations.append(annot)
				
	# Annotate primer locations
	for i, primer in enumerate(primer_list):
		idx_start = re.search(primer, new_seq, re.IGNORECASE).start()
		idx_end = idx_start + len(primer)
		annot = {'start':idx_start, 'end':idx_end}
		annot['strand'] = -2*i + 1
		annot['name'] = 'Primer_{}'.format(i+1)
		annot['color'] = color
		annot['type'] = 'misc_feature'
		new_annotations.append(annot)

	new_translations = []
	for trans in base_translations:
		trans_direction = trans['strand']
		trans = trans['regions']
		assert(len(trans) == 1), 'Something is strange in the translation. Check it out.'
		trans = trans[0]
		trans_len = trans['end'] - trans['start']
		trans['strand'] = trans_direction
		
		if trans['start'] <= cut_idx[0] and trans['start'] + trans_len > cut_idx[0]:
			temp_trans = dict(trans)
			temp_trans['end'] = cut_idx[0] - (cut_idx[0] - trans['start'])%3
			new_translations.append(temp_trans)
		if trans['end'] >= cut_idx[1] and trans['end'] - trans_len < cut_idx[1]:
			temp_trans = dict(trans)
			temp_trans['start'] = cut_idx[0] + gene_len
			temp_trans['end'] += gene_len - remove_len
			new_translations.append(temp_trans)
		if trans['end'] <= cut_idx[0]:
			new_translations.append(dict(trans))
		elif trans['start'] >= cut_idx[1]:
			temp_trans = dict(trans)
			temp_trans['start'] += gene_len + insert_len - remove_len
			temp_trans['end'] += gene_len + insert_len - remove_len
			new_translations.append(temp_trans)

	# Get translation from insert
	for trans in insert_translations:
		print trans
		trans_direction = trans['strand']
		trans = trans['regions']
		
		assert(len(trans) == 1), 'Something is strange in the translation. Check it out.'
		trans = trans[0]
		trans_len = trans['end'] - trans['start']
		trans['strand'] = trans_direction * direction
		
		# Check if regions overlap
		print(trans['start'], trans['end'])
		print(idx_list)
		if trans['start'] <= idx_list[0] and idx_list[1] <= trans['end']:
			
			# Cut out the region of translation in the insert
			if trans['start'] <= idx_list[0]:
				trans['start'] = idx_list[0]
			trans['end'] = min(trans['end'], idx_list[1])
			trans['end'] -= (trans['end'] - trans['start'])%3
			
			if direction == -1:
				start_new = gene_len - trans['end'] + idx_list[0] 
				end_new = gene_len + idx_list[0] - trans['start']
			else:
				start_new = trans['start'] - idx_list[0]
				end_new = trans['end'] - idx_list[0]
			
			# Now add the offset from insertion
			trans['start'] = start_new + cut_idx[0] + insert_left_len
			trans['end'] = end_new + cut_idx[0] + insert_left_len
			new_translations.append(trans)
	print new_translations
	#return
	folderId = get_folder_id(FOLDER)
	# Save on benchling
	headers = {'Content-Type' : 'application/json'}
	data = {"aliases": [], "annotations": new_annotations, "bases": new_seq, 
		"customFields": {"Submitter": {"value": SUBMITTOR} }, "fields": {},
		"folderId": folderId, "isCircular": True, "name": new_plasmid_name,
		"schemaId": None, "translations": new_translations,}
	url = 'https://benchling.com/api/v2/dna-sequences'
	response = requests.post(url, headers=headers, data=json.dumps(data), auth=(SECRET_KEY, ''))