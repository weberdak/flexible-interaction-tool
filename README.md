# flexible-interaction-tool
FLexible InteRaction Tool (FLIRT) VMD script for summarizing interactions in MD simulations.

## Overview

FLIRT (FLexible InteRaction Tool) is an analysis script for implementation with VMD software to rapidly (seconds to minutes) analyse general contacts, hydrogen bonds, cation-π interactions and hydrophobic (non-polar) contacts for a single PDB file, or over entire MD trajectories (thousands of snapshots). Interactions between atoms of two selections, within user-specified geometric constraints (i.e., cutoff distances and angles), are tallied on a per-atom basis and accumulated into any desired attribute (i.e., VMD attributes are index, resid, serial, residue, resname, type, name, structure) for the first selection, or decomposed into pairwise interactions with any attribute of a second selection. The results for either accumulated or pairwise measurements may be summarized as a mean and standard deviation over all frames analysed (*.xy,dat and *.pairxy.dat files, respectively) or time-resolved over an entire trajectory (*.xyz.dat and *.pairxyz.dat files, respectively).

## Usage

The script (flirt.tcl) must be copied into the working directory and loaded into VMD using the Tk Console and command ‘source flirt.tcl’. A measurement may be executed by typing ‘flirt selection1 selection2 measurement [options]’ into the console, selecting only ONE measurement and any number of options, described below:

#### Contacts

* -contact cutoff [options]

Extension of the ‘measure contacts’ command part of VMD. Atom pairs between selection1 and selection2, within cutoff distance (Å), are detected and summarised according to the options provided.

Example:

	source flirt.tcl
	
	set proa [atomselect top "segname PROA"]
	set prob [atomselect top "segname PROB"]
	
	flirt $proa $prob -contact 3.5 -xy -prefix proteins.contacts -zeros
	

#### H-bonds

* -hbond cutoff angle [options]

Extension of the ‘measure hbonds’ command shipped with VMD. A hydrogen bond is detected if the distance between donor and acceptor (non-hydrogen) atoms is within the cutoff distance (Å), and the donor-hydrogen-acceptor angle 2 is less than angle from 180 degrees. Values of 3 Å and 20° are typical. By default, a hydrogen bond is only detected if the donor and acceptor atoms are part of selection1 and selection2, respectively. The ‘-combine’ option may be specified if both donor and acceptor interactions are to be tallied regardless of selection order. It is also a known issue in VMD that ‘hbond’ measurements may produce artefacts if selection1 and selection2 share a common set of atoms. In this case, the ‘-exclude attribute’ option may be used. For example, ‘-exclude resid’ will ignore all hydrogen bonds where donors and acceptor atoms share the same residue ID number.

Example:

	source flirt.tcl
	
	set proa [atomselect top "segname PROA"]
	set prob [atomselect top "segname PROB"]
	
	flirt $proa $prob -hbond 3 20 -xy -pairxy -prefix proteins.hbonds


#### Cation-Pi interactions

* -cationpi cutoff1 cutoff2 [options]

A cation-π interaction is detected between selection1 and selection2 if all carbons of a six-membered ring system are within cutoff1 (Å) to a cationic centre, and if distances measured for all carbons differ no greater than cutoff2 (Å). Values of 7 Å and 1.5 Å for cutoff1 and cutoff2, respectively, have been used by (Grauffel et al., 2013) to approximate the maximum interaction distance and a 120° selection cone for the cation position with respect to the plane of the ring system (Minoux and Chipot, 1999, Petersen et al., 2005). Ring systems that are not part of TRP, TYR or PHE residues, and cations that are not part of ARG, LYS, DPC or ASM residues, must defined in the ‘proc flirt::check_cationpi’ section of the code for the interaction to be detected (note that the code in its current state doesn’t allow for residues which have multiple ring systems). By default, selection1 is scanned for π-systems and selection2 for cations, but, again, this condition may be relaxed using the ‘-combine’ option. Also, a ‘measure_cationpi cutoff1 cutoff2 selection1 selection2’ procedure is made available, which returns a list of π-atom indices (six unique entries) from selection1 and paired cations for selection2 (atom index is entered six times to match size of first list). 

Example:

	source flirt.tcl
	
	set proa [atomselect top "segname PROA"]
	set prob [atomselect top "segname PROB"]
	
	flirt $proa $prob -cationpi 7 1.5 -xy -pairxy -prefix proteins.cationpi
	
#### Hydrophobic interactions

* -nonpolar cutoff

Is a modification of the ‘-contacts’ measurement, but ignores any interaction from atoms that are not assigned as a nonpolar carbon or hydrogen. This function requires that the ‘flirt::check_hydrophobic’ procedure is updated with any non-standard residues (i.e., not a common amino acid). Also, a ‘measure_nonpolar cutoff sel1 sel2’ procedure is introduced, which outputs non-polar contacts in a selected frame. Atom indices of selection1 are in the first list alongside paired atoms of selection2 in a second list.

Example:

	source flirt.tcl
	
	set proa [atomselect top "segname PROA"]
	set prob [atomselect top "segname PROB"]
	
	flirt $proa $prob -nonpolar 3.5 -xy -pairxy -prefix proteins.nonpolar
	

#### Additional options

Options ([options]) may be specified in any order and include:

* -attr1 attribute1

Interactions involving atoms of selection1 are accumulated into a shared attribute1 every frame. Supported attributes include: index, resid, serial, residue, resname, type, name and structure. Default: resid.

* -attr2 attribute2

Accumulate interactions involving selection2 atoms by attribute2. Nothing is done unless -pairxy and/or -pairxyz are also stated. Default resid.

* -prefix prefix

Files outputted by the ‘-xy’, ‘-xyz’, ‘-pairxy’, ‘-pairxyz’ and ‘-total’ flags are named with a prefix prefix. Defaults: ‘total’, ‘contacts’, ‘hbonds’, ‘cationpi’ or ‘nonpolar’.

* -xy

Outputs a file ‘prefix.xy.dat’ with the attribute1(s) of selection1 (first column), the mean number of interactions per frame (second column) and standard deviations (third column).

* -xyz

Outputs a file ‘prefix.xyz.dat’ with frame numbers (first column), attribute1(s) of selection1 (second column) and the number of interactions per frame (third column). A blank line is inserted to create a new block of data for every frame for compatibility with the pm3d 3D-plotting function of Gnuplot.

* -pairxy

Outputs a file ‘prefix.pairxy.dat’ of summarised pairwise interactions between attribute1(s) of selection1 (first column) and attribute2(s) selection2 (second column) according to the mean number of interactions per frame (third column) and standard deviation (fourth column).

* -pairxyz

Outputs a file ‘prefix.pairxyz.dat’ with frame numbers (first column), pairwise interactions between selection1 attribute1(s) and selection2 attribute2(s) (second column, formatted as attribute1.attribute2) and the number of interactions per frame (third column).

* -total

Outputs a file ‘prefix.total.dat’ of the total interactions detected between atoms selection1 and selection2 for every frame.

* -start frame

Starting frame of analysis. Default: 0.

* -stop frame

End frame of analysis. Default: last frame.

* -skip frame

Step every frame. Default: 1.

* -betas

A beta value, useful for colour scaling, is assigned to each atom in selection1 according to the average number of interactions to its corresponding attribute1 (requires option ‘-attr1 attribute1’).

* -zeros

By default, no value is written to ‘.xy.dat’ and ‘.xyz.dat’ files if no interactions are detected. If zero values are necessary, as is the case for some plotting programs, they can be written by specifying this flag. This flag is ignored for pairwise outputs.

* -combine

Is used to consider both donor and acceptor contributions to hydrogen bonding by selection1 and selection2. May also be used to detect cation-π interactions where both selections contain π-systems and cations. This flag is ignored for ‘-contact’ and ‘-nonpolar’ measurements. 

* -exclude attribute

Interactions are ignored when atoms of selection1 and selection2 share the same attribute.

* -mol mol

Specify this flag only if the molecule isn’t ‘top’. Default: top

## References

* GRAUFFEL, C., YANG, B., HE, T., ROBERTS, M. F., GERSHENSON, A. & REUTER, N. 2013. Cation-π interactions as lipid-specific anchors for phosphatidylinositol-specific phospholipase C. Journal of the American Chemical Society, 135, 5740-5750.

* MINOUX, H. & CHIPOT, C. 1999. Cation-π interactions in proteins: Can simple models provide an accurate description? Journal of the American Chemical Society, 121, 10366-10372.

* PETERSEN, F. N. R., JENSEN, M. Ø. & NIELSEN, C. H. 2005. Interfacial tryptophan residues: A role for the cation-π effect? Biophysical Journal, 89, 3985-3996.
