import warnings
warnings.filterwarnings("ignore")

import argparse

ap=argparse.ArgumentParser(description=f"""This program is written for the identification of AA-N-AA stacking interactions on a given mmcif file.
Please provide the required positional arguments correctly.""")
ap.add_argument("-i", dest="input_struct", required=True, type=str, help= "Please enter the name of the mmcif file. Standard format: \"-i example.cif\".")
ap.add_argument("-d", dest="provided_distance", required=True, type=float, help="Please enter the distance. For example: \"-d 5.0\".")
ap.add_argument("-a", dest="tilt_angle", required=True, type=float, help="Please enter the maximum tilt angle. For example: \"-a 40.0\".")

all_args=ap.parse_args()
input_name= all_args.input_struct
input_distance= all_args.provided_distance
input_t_angle= all_args.tilt_angle


from Bio.PDB.MMCIFParser import MMCIFParser
parser= MMCIFParser()
cif_input= parser.get_structure("input", input_name)

def ring_atoms_gen(residue):
  if residue.get_resname()== "TRP":
    return [residue["CG"], residue["CD1"], residue["NE1"] ,residue["CD2"], residue["CE2"], residue["CE3"], residue["CZ2"], residue["CZ3"], residue["CH2"]]
  elif residue.get_resname()=="TYR":
    return [residue["CG"], residue["CD1"], residue["CD2"], residue["CE1"], residue["CE2"], residue["CZ"]]
  elif residue.get_resname()=="PRO":
    return [residue["CA"], residue["CB"], residue["CG"], residue["CD"], residue["N"]]
  elif residue.get_resname()=="PHE":
    return [residue["CG"], residue["CD1"], residue["CD2"], residue["CE1"], residue["CE2"], residue["CZ"]]
  elif residue.get_resname()=="HIS":
    return [residue["CG"], residue["ND1"], residue["CD2"], residue["CE1"], residue["NE2"]]
  elif residue.get_resname() in ["C", "U"]:
    return [residue["N1"], residue["C2"], residue["N3"], residue["C4"], residue["C5"], residue["C6"]]
  elif residue.get_resname() in ["A", "G"]:
    return [residue["N1"], residue["C2"], residue["N3"], residue["C4"], residue["C5"], residue["C6"], residue["N7"], residue["C8"], residue["N9"]]
  else:
    return []

import numpy as np
def centroid_calc(ring_atom_list_variable):
  coordinates= np.array([atom.get_coord() for atom in ring_atom_list_variable])
  return coordinates.mean(axis=0)

from math import sqrt
def distance_calc(centroid_1, centroid_2):
  return sqrt((centroid_1[0]-centroid_2[0])**2 + (centroid_1[1]-centroid_2[1])**2 + (centroid_1[2]-centroid_2[2])**2)


def angle_calc(plane1, plane2):
  v1_amino= plane1[1]-plane1[0]
  v2_amino= plane1[2]-plane1[0]
  n1_amino = np.cross(v1_amino, v2_amino)

  v3_nucleo = plane2[1]-plane2[0]
  v4_nucleo = plane2[2]-plane2[0]
  n2_nucleo = np.cross(v3_nucleo, v4_nucleo)
  angle_radian= np.arccos(np.dot(n1_amino, n2_nucleo)/(np.linalg.norm(n1_amino)*(np.linalg.norm(n2_nucleo))))
  angle_degree= np.degrees(angle_radian)
  return angle_degree


def centroid_to_plane_dist_calc(centroid, plane):
  a= plane[0]-plane[1]
  b= plane[0]-plane[2]
  N= np.cross(a,b)

  d= -np.dot(N, plane[0])
  numerator = np.abs(np.dot(N, centroid)+d)
  denominator = np.linalg.norm(N)
  final_distance= numerator/denominator
  return final_distance



def all_stacks(input_s, distance, angle_t):
  stacking_aminos= []
  nucleotides=[]
  for model in input_s:
    for chain in model:
      for residue in chain:
        if residue.get_resname() in ["PHE", "TYR", "TRP", "HIS", "PRO"]:
          stacking_aminos.append(residue)
        elif residue.get_resname() in ["A", "G", "C", "U"]:
          nucleotides.append(residue)
        else:
          continue
  stacks=[]
  for n in nucleotides:
    ring_atoms_nucleo= ring_atoms_gen(n)
    centroid_nucleo = centroid_calc(ring_atoms_nucleo)

    nearby_stacking_aminos =[]
    for amino in stacking_aminos:
      ring_atoms_amino= ring_atoms_gen(amino)
      centroid_amino= centroid_calc(ring_atoms_amino)
      distance_between_nucleo_amino= distance_calc(centroid_nucleo, centroid_amino)

      if distance_between_nucleo_amino <= distance:
        angle_between_planes= angle_calc(np.array([atom.get_coord() for atom in ring_atoms_amino])[:3], np.array([atom.get_coord() for atom in ring_atoms_nucleo])[:3])
        distance_centro_nucleo_amino_plane= centroid_to_plane_dist_calc(centroid_nucleo, np.array([atom.get_coord() for atom in ring_atoms_amino])[:3])
        distance_centro_amino_nucleo_plane=centroid_to_plane_dist_calc(centroid_amino, np.array([atom.get_coord() for atom in ring_atoms_nucleo])[:3])
        if angle_between_planes <= angle_t or angle_between_planes >= (180.0-angle_t):
          if distance_centro_nucleo_amino_plane >= (0.7* distance_between_nucleo_amino) and distance_centro_amino_nucleo_plane>= (0.7 * distance_between_nucleo_amino):
            nearby_stacking_aminos.append((amino, distance_between_nucleo_amino, angle_between_planes, distance_centro_nucleo_amino_plane, distance_centro_amino_nucleo_plane))
    if len(nearby_stacking_aminos)>=2:
      nearby_stacking_aminos.sort(key=lambda x: x[1])
      amino1= nearby_stacking_aminos[0][0]
      amino2= nearby_stacking_aminos[1][0]

      stacks.append(f"""amino_acid_1 : {amino1.parent.id}.{amino1.get_resname()}.{amino1.get_id()[1]}
nucleotide   : {n.parent.id}.{n.get_resname()}.{n.get_id()[1]}
amino_acid_2 : {amino2.parent.id}.{amino2.get_resname()}.{amino2.get_id()[1]}

calculated values :
amino_acid_1 :
d     : {nearby_stacking_aminos[0][1]}
angle : {nearby_stacking_aminos[0][2]}
n1    : {nearby_stacking_aminos[0][4]}
n2    : {nearby_stacking_aminos[0][3]}

amino_acid_2 :
d     : {nearby_stacking_aminos[1][1]}
angle : {nearby_stacking_aminos[1][2]}
n1    : {nearby_stacking_aminos[1][4]}
n2    : {nearby_stacking_aminos[1][3]}
""")
  return stacks


print("""The output format :
stack no : i
amino_acid_1 : chain.residue_name.residue_number
nucleotide   : chain.residue_name.residue_number
amino_acid_2 : chain.residue_name.residue_number

calculated values :
amino_acid_1 :
d     : distance between the centroids of amino_acid_1 ring/rings and the nucleotide ring/rings (in angstroms)
angle : angle between the planes of amino_acid_1 ring/rings and the nucleotide ring/rings (in degrees)
n1    : normal distance between the centroid of amino_acid_1 ring/rings and the plane of nucleotide ring/rings (in angstroms)
n2    : normal distance between the centroid of nucleotide ring/rings and the plane of amino_acid_1 ring/rings (in angstroms)

amino_acid_2 :
d     : distance between the centroids of amino_acid_2 ring/rings and the nucleotide ring/rings (in angstroms)
angle : angle between the planes of amino_acid_2 ring/rings and the nucleotide ring/rings (in degrees)
n1    : normal distance between the centroid of amino_acid_2 ring/rings and the plane of nucleotide ring/rings (in angstroms)
n2    : normal distance between the centroid of nucleotide ring/rings and the plane of amino_acid_2 ring/rings (in angstroms)
.......................................................................................................................................

For your given parameter, the identified stacks are :
""")

stacks_out= all_stacks(cif_input, input_distance, input_t_angle)
i=1
for stack in stacks_out:
  print(f"""stack no : {i}
{stack}
.............................................................................................................................................
""")
  i=i+1