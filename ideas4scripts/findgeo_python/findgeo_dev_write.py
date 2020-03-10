#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2010-2011
# findgeo is distributed under the terms of the GNU General Public License
#
"""
    findgeo is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    findgeo is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os
import shutil
import subprocess
#import subprocess
from p3d.protein import Protein 
from optparse import OptionParser
import urllib


# This script analyzes the PDB coordinate files and extracts sites

def main(code, input_file, output_dir, distance, metal):
	
	ref_metalList = ['LI', 'BE', 'NA', 'MG', 'AL', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 
	                'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'RB', 'SR', 'Y', 'ZR', 'NB', 
	                'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'CS', 'BA', 
	                'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 
	                'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 
	                'TL', 'PB', 'BI', 'PO', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP', 'PU', 
	                'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG',
	                'AS']
	
	"""Main function. Prepare the find geometry steps."""

	if input_file:
		code = os.path.join(output_dir, '%s' % code)
		if not os.path.isfile(code):
			sys.exit('%s: file not found' % code)
	else:
		pdb_code = code.lower()
		code = os.path.join(output_dir, 'pdb%s.ent.gz' % pdb_code)
		try:
			urllib.urlretrieve('ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz' % pdb_code, code)
		except:
			sys.exit('Unable to retrieve %s' % pdb_code)
	try:
		pdb = Protein(code)
	except:		
		sys.exit('PDB in bad format.')
		
	if metal:
		metalList = []
		metalList.append(metal)
	else:
		metalList = ref_metalList

	metals = find_metals(pdb, metalList)

	if metals:
		out = open(os.path.join(output_dir, 'FindGeo.summary'), 'w')
		ligands = find_ligands(metals, pdb, metalList, distance, not_donors)
		for key in ligands:
			print("-"*40)
			#print("%s_%s%s_%s_%s --> %d ligands found.",%key.atype.strip(), key.resid, key.altConf, key.idx, key.chain, len(ligands[key]))
			if len(ligands[key]) >= 2 and len(ligands[key]) <= 9 :
				print('Determining coordination geometry...')
				pdb_name = os.path.join(output_dir, '%s_%s_%s_%s.pdb' % (key.atype.strip(), key.resid, key.idx, key.chain))
				outfile = os.path.join(output_dir, pdb_name)
				#temp = open(outfile, 'a')
				temp = open(outfile, 'w')
				for line in pdb.query(key):
				 	temp.writelines('%s\n' % line.output())
					#previous line writes findgeo_input
				for element in ligands[key]:
					temp.writelines('%s\n' % element.output())
					#previous line writes findgeo_input
				temp.close()
				find_geo_input(pdb_name, output_dir, out)
			elif len(ligands[key]) < 2:
				print ("Too few donor atoms identified. At least 2 donor atoms must be present.")
			else:
				print ("Too many donor atoms identified. The maximum allowed coordination number is 9.")
		out.close()
		print ("-"*40)
	else:
		if metal:
			sys.exit('%s not found in the PDB file' % metal)
		sys.exit('No metals found in the PDB file')


def find_metals(pdb, metalList):
	"""Find metals in the pdb"""
	metals = [atom for atom in pdb.query('non-protein and model 1') if atom.elementType.upper().strip() in metalList]
	return metals

def find_ligands(metals, pdb, metalList, maxDist, not_donors):
	"""Find ligands in the pdb"""
	ligands = {}
	
	for metal in metals:
		ligands[metal] = []
		temp = pdb.query('model 1 & within {0} of'.format(maxDist), metal)
		
		for atom in temp:
			if atom.elementType.upper().strip() not in not_donors and atom.elementType.upper().strip() not in metalList:
				ligands[metal].append(atom)					
		
		if not ligands[metal]:
			print("No ligands found for metal %s_%s_%s_%s (try to change distance threshold (-t) or excluded donors (-e))" % (metal.atype.strip(), metal.resid, metal.idx, metal.chain))
	
	return ligands

def find_geo_input(pdb_name, output_dir, out):
	"""Find geometries and print the best"""
	
	geometries = {'lin': 'linear',
              'irr': 'irregular',
              'trv': 'trigonal plane with a vacancy',
              'tri': 'trigonal plane',
              'tev': 'tetrahedron with a vacancy',
              'spv': 'square plane with a vacancy',
              'tet': 'tetrahedron',
              'spl': 'square plane',
              'bva': 'trigonal bipyramid with a vacancy (axial)',
              'bvp': 'trigonal bipyramid with a vacancy (equatorial)',
              'pyv': 'square pyramid with a vacancy (equatorial)',
              'spy': 'square pyramid',
              'tbp': 'trigonal bipyramid',
              'tpv': 'trigonal prism with a vacancy',
              'oct': 'octahedron',
              'tpr': 'trigonal prism',
              'pva': 'pentagonal bipyramid with a vacancy (axial)',
              'pvp': 'pentagonal bipyramid with a vacancy (equatorial)',
              'cof': 'octahedron, face monocapped with a vacancy (capped face)',
              'con': 'octahedron, face monocapped with a vacancy (non-capped face)',
              'ctf': 'trigonal prism, square-face monocapped with a vacancy (capped face)',
              'ctn': 'trigonal prism, square-face monocapped with a vacancy (non-capped edge)',
              'pbp': 'pentagonal bipyramid',
              'coc': 'octahedron, face monocapped',
              'ctp': 'trigonal prism, square-face monocapped',
              'hva': 'hexagonal bipyramid with a vacancy (axial)',
              'hvp': 'hexagonal bipyramid with a vacancy (equatorial)',
              'cuv': 'cube with a vacancy',
              'sav': 'square antiprism with a vacancy',
              'hbp': 'hexagonal bipyramid',
              'cub': 'cube',
              'sqa': 'square antiprism',
              'boc': 'octahedron, trans-bicapped',
              'bts': 'trigonal prism, square-face bicapped',
              'btt': 'trigonal prism, triangular-face bicapped',
              'ttp': 'trigonal prism, square-face tricapped',
              'csa': 'square antiprism, square-face monocapped'
              }

	#next 3 lines are mine
	#shutil.rmtree('findgeo_results', ignore_errors=True)
	#dir_name = 'findgeo_results'
	#os.mkdir(dir_name)

	dir_name = pdb_name.split('.pdb')[0]
	temp = dir_name.split('/')
	name = temp[len(temp)-1]

	#overrode safety check, WILL OVERWRITE PREVIOUS RESULTS
	'''if os.path.exists('%s' % dir_name):
		if overwrite:
			try:
				shutil.rmtree('%s' % dir_name)
			except:
				sys.exit('Cannot remove the existing directory %s: check permissions' % dir_name)
		else:
			sys.exit('%s: is an existing directory... Remove/rename it or use the -o option' % dir_name)
	try:
		os.mkdir(dir_name)
	except:
		sys.exit('Cannot create the directory %s: check permissions' % dir_name)'''

	shutil.move(pdb_name, os.path.join(dir_name, 'findgeo.input'))	
	cmd = './findgeo %s' % dir_name
	#next line is mine
	#cmd = './findgeo findgeo_results'

	geom = subprocess.getoutput(cmd).split()


	if geom[0] != 'irr':
		print('placeholder')
		out.write('%s placeholder' % name)
		
		#print('Best fit geometry: %s %s' % (geometries[geom[0]], geom[1]))
		#out.write('%s: %s - %s %s\n' % (name, geom[0], geometries[geom[0]], geom[1]))
	else:
		print ('Irregular geometry' )
		#out.write('%s: %s - Irregular geometries\n' % (name, geom[0]))
		out.write('%s: - Irregular geometry \n' % (name))


if __name__ == '__main__':

	usage = "usage: %prog -p pdbfile [-c pdbcode] [-t threshold] [-m metal] [-e excluded_donors] [-w workdir] [-o overwrite]"
	
	parser = OptionParser(usage)
	parser.add_option("-w", "--wdir", dest="wd",
					help="Directory where to find or download the input PDB file and to write outputs. Default is ./",
					metavar="WORKDIR")
	
	parser.add_option("-p", "--pdb_file", dest="pdb_file",
					help="Local input PDB file.", 
                                        metavar="FILE")
	
	parser.add_option("-c", "--pdb_code", dest="pdb_code",
					help="PDB code of input PDB file to be downloaded from the web.", 
                                        metavar="PDBCODE")
	
	parser.add_option("-t", "--threshold", dest="threshold",
					help="Coordination distance threshold. Default is 2.8 A.",
					metavar="THRESHOLD")
	
	parser.add_option("-m", "--metal", dest="metal",
					help="Chemical symbol of the metal of interest. Default is all metals.",
					metavar="METAL")
	
	parser.add_option("-e", "--excluded_donors", dest="notdonors",
					help="Chemical symbols of the atoms (separated by commas) excluded from metal ligands. Default is C and H.",
					metavar="NOTDONORS")
	
	parser.add_option("-o", "--overwrite", dest="overwrite",
					action="store_true", default=False,
	                                help="Overwrite existing files and directories.",
                                        metavar="OVERWRITE")
	
	(options, args) = parser.parse_args()
	
	if not options.pdb_file and not options.pdb_code:
		parser.error("A PDB file or PDB code is required.")
	
	if  options.pdb_file and options.pdb_code:
		parser.error("A PDB file or PDB code is required, not both.")
	
	if options.wd:
		wd = os.path.abspath(options.wd) + '/'
	else:
		wd = './'
	
	if options.metal:
		metal = options.metal.upper()
		if metal not in ref_metalList:
			parser.error('%s: invalid metal.' % metal)
	else:
		metal = None
		
	if options.threshold:
		try:
			threshold = float(options.threshold)
		except ValueError:
			sys.exit('Invalid threshold. This must be a number.')
	else:
		threshold = 2.8
		
	if options.pdb_file:
		pdb = options.pdb_file
		input_file = True
	elif options.pdb_code:
		pdb = options.pdb_code
		input_file = False
		if len(pdb) != 4:
			parser.error('Invalid PDB code. Must be four characters.')

	if options.notdonors:
		not_donors = options.notdonors.replace(' ','').split(',')
	else:
	        not_donors = ['C', 'H']
		
	#download = options.download
	overwrite = options.overwrite
	
	main(pdb, input_file, wd, threshold, metal)

