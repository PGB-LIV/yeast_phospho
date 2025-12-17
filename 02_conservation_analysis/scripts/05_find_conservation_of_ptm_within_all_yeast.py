################################################################################## 
# This script is for conservation analysis                                       # 
# at PTM and proximal pos:                                                       #
# Read in alignment                                                              #
# For target, read in positions of interest                                      #
# Need to loop through alignment for target sequence, then create dictionary     #
# mapping aligned positions to target positions                                  # 
# Check amino acid is as expected                                                #
##################################################################################
from Bio import AlignIO
import os


### Reads protein MSA and extracts the amino acids in the alignment at the given target positions of the target protein
### Returns a dictionary of results
def read_alignment(alignment_file,target_protein,target_pos_dict):

    print("Reading",alignment_file)
    alignment = AlignIO.read(alignment_file, "fasta")
    amino_acids_found_at_pos = {}
    proteins_found_at_pos = {}
    amino_acids_at_target_positions_plus_one = {}
    all_original_protein_positions_for_alignment_pos = {} #Will fill this with the source positions in proteins

    aligned_target_protein_found = False
    aligned_target_protein = None

    #aligned_record_to_pos_in_alignment = -1
    # find target protein in alignment file
    for i in range(0, len(alignment)):
        aligned_record = alignment[i]
        if aligned_record.id == target_protein:
            aligned_target_protein = aligned_record
            aligned_target_protein_found = True
            aligned_record_to_pos_in_alignment = i

    if aligned_target_protein_found == True:
        # yeast prot seq 
        aligned_prot_seq = str(aligned_target_protein.seq)
        prot_seq = aligned_prot_seq.replace("-","")

        # loop through aa's
        aa_counter = 0
        for pos in range(0,len(aligned_prot_seq)):
            aa = aligned_prot_seq[pos:pos+1]
            current_aligned_prot_seq = aligned_prot_seq[:pos+1] # protein seq up to aa
            aa_counter = len(current_aligned_prot_seq) - current_aligned_prot_seq.count("-") -1 #0 based position

            # if aa is a target phosphosite
            if aa_counter+1 in target_pos_dict and aa != "-": #target_pos_dict starts at 1
                #residue_in_dict = target_pos_dict[aa_counter+1]
                #print("aa counter",aa_counter,residue_in_dict ,"pos->",pos)

                # get aa in seq without "-"
                aa = prot_seq[aa_counter:aa_counter+1]
                # get aa in aligned seq
                aligned_aa = aligned_prot_seq[pos:pos+1]
                #print("found",pos,aa_counter,aligned_aa,aa)
                # check aa is the same in both
                if aa != aligned_aa:
                    print("error at pos",pos,"source pos:",aa_counter, target_protein,alignment_file,target_protein," ",aa," doesn't match",aligned_aa)
                #else:
                #    print("match:",aa,aligned_aa,pos,str(aa_counter),"\n")

                proteins_in_alignment = []
                aas_found_at_pos = {}
                aas_found_at_proximal_pos = {}
                original_protein_positions_for_alignment_pos = {} #put PTM pos (starting from zero) into this list
                original_protein_positions_for_alignment_pos[target_protein] = aa_counter

                # Get position of next amino acid
                found = False
                proximal_pos = -1
                proximal_residue = "X"
                j = pos + 1
                while j < len(aligned_target_protein.seq) and found == False:
                    if aligned_target_protein.seq[j] != "-":
                        found = True
                        proximal_pos = j
                        #print("Found proximal pos",proximal_pos)
                        proximal_residue = aligned_target_protein.seq[j]
                        # add to dict, K= yeast prot ID, V= proximal res
                        aas_found_at_proximal_pos[aligned_target_protein.id] = proximal_residue #Always put first
                        #print("proximal ref residue",aligned_target_protein.seq[j],aligned_target_protein.id)
                    j += 1

                aas_found_at_pos[target_protein] = alignment[aligned_record_to_pos_in_alignment].seq[pos]
                #print("target prot",target_protein)
                proteins_in_alignment.append(target_protein)  #Always put first
                # proteins in aln file
                for i in range(0, len(alignment)):
                    # if it is not the target
                    if i != aligned_record_to_pos_in_alignment:
                        # ortho id
                        aligned_prot = alignment[i].id
                        proteins_in_alignment.append(aligned_prot)
                        aas_found_at_pos[aligned_prot] = alignment[i].seq[pos]

                        aligned_seq_up_to_pos = alignment[i].seq[:pos+1]
                        original_protein_positions_for_alignment_pos[alignment[i].id] = len(aligned_seq_up_to_pos) - aligned_seq_up_to_pos.count("-") #Leave as 1 based position

                        if proximal_pos != -1 and i != aligned_record_to_pos_in_alignment: #Proximal residue is always first
                            aas_found_at_proximal_pos[aligned_prot] = alignment[i].seq[proximal_pos]


                amino_acids_found_at_pos[aa_counter+1] = aas_found_at_pos #Target positions start at 1
                proteins_found_at_pos[aa_counter+1] = proteins_in_alignment
                amino_acids_at_target_positions_plus_one[aa_counter+1] = aas_found_at_proximal_pos
                all_original_protein_positions_for_alignment_pos[aa_counter+1] = original_protein_positions_for_alignment_pos


            #else:
            #    print("not found",pos,aa_counter)

    else:
        print("Error, no target protein found",target_protein, " in alignment file",alignment_file)
    return [amino_acids_found_at_pos,proteins_found_at_pos,amino_acids_at_target_positions_plus_one,all_original_protein_positions_for_alignment_pos]

def calculate_conservation_scores_per_ptm_site(protein_id,pos_res_dict,amino_acids_at_target_positions,proteins_at_target_positions,amino_acids_at_target_positions_plus_one,all_original_protein_positions_for_alignment_pos):

    f_out = open(DATA_FOLDER + protein_id + "_stats.tsv","w")
    for pos in amino_acids_at_target_positions:

        protein_list = proteins_at_target_positions[pos]
        aa_list = amino_acids_at_target_positions[pos]
        proximal_aa_list = amino_acids_at_target_positions_plus_one[pos]
        original_positions_for_pos = all_original_protein_positions_for_alignment_pos[pos]

        target_res = ""
        if pos not in pos_res_dict:
            print("error position not in starting dictionary",pos,pos_res_dict)
        else:
            target_res = pos_res_dict[pos]
        #print(pos," ->",str(aa_list), target_res)

        match_count = 0
        soft_match_count = 0
        proximal_res_match_count = 0

        reference_protein = protein_list[0]
        reference_proximal_res = "X"    #Not found code
        if reference_protein in proximal_aa_list:
            reference_proximal_res = proximal_aa_list[reference_protein ]
        #print("ref aa",reference_proximal_res)
        mapped_protein_sites = []
        all_residues = []
        all_proximal_residues = []


        for protein in protein_list:
            if protein != reference_protein:
                res = aa_list[protein]
                all_residues.append(res)

                proximal_res = "-1"  #Error code
                if protein in proximal_aa_list:
                    proximal_res = proximal_aa_list[protein]
                    all_proximal_residues.append(proximal_res)

                soft_match = False
                proximal_match = False
                if res == target_res:
                    match_count += 1
                if (target_res == "S" and res == "T") or (target_res == "T" and res == "S") or res == target_res:
                    soft_match_count +=1
                    soft_match = True
                if proximal_res == reference_proximal_res:
                    proximal_res_match_count += 1
                    proximal_match = True

                if proximal_match and soft_match and protein !=reference_protein:
                    mapped_protein_sites.append(protein + "_" + str(original_positions_for_pos[protein]))

        
         #230924 added soft match count /uniq species in ortho file ie 272
        f_out.write(protein_id + "\t" +
                    str(pos) + "\t" +
                    target_res + "\t" +
                    reference_proximal_res + "\t" +
                    str(len(aa_list)-1) + "\t" +
                    str(match_count) + "\t"+
                    str(soft_match_count) + "\t" +
                    str(proximal_res_match_count) + "\t" +
                    ";".join(all_residues) +"\t"+
                    ";".join(all_proximal_residues) + "\t" +
                    ";".join(protein_list) +"\t"+
                    ";".join(mapped_protein_sites) + "\t" +
                    str(soft_match_count/272) + "\t" +
                    str(match_count/(len(aa_list)-1)) + "\n")
    f_out.close()


### function for reading all target sites from a single file
def get_all_target_sites(protein_ptm_file,prot_col,pos_col,residue_col):
    f = open(protein_ptm_file, "r")
    protein_to_positions = {}

    counter = 0
    for line in f:
        if counter != 0:
            line = line.rstrip()
            cells = line.split("\t")
            protein = cells[prot_col]
            pos = int(cells[pos_col])
            res = cells[residue_col]
            #print(protein, pos, res)

            pos_to_res_dict = {}
            # Is the protein already in dict
            if protein in protein_to_positions:
            #if yes, get back the pos->res dict values for that protein. If no we will use the empty dict defined above
                pos_to_res_dict = protein_to_positions[protein]
            # add new entry, key is the pos and val is residue
            pos_to_res_dict[pos] = res
            # add back to big dict where key is the protein id and the values are the dictionary with the pos->res phospho mod
            protein_to_positions[protein] = pos_to_res_dict
        #else:
            #print("Ignoring non-canonical model",cells[prot_col])
        counter += 1

    return protein_to_positions



DATA_FOLDER = "outputs/conservation_stats_yeast/"
ALIGNMENT_FOLDER = "outputs/fastas_pre_aligment_syntenic_Sept24/fastas_yeast/"

if os.path.exists(DATA_FOLDER) == False:
    os.mkdir(DATA_FOLDER)


TARGET_FILE = "outputs/GSB_STY_conservation_yeast.tsv"

protein_to_targets = get_all_target_sites(TARGET_FILE,4,1,2)
missing_pos = []
count_missing_alignnment=0
test_counter = 0
# For each protein in the GSB file (ie our file where we have found phosphosites)
for protein_id in protein_to_targets:
    # get back the values for the protein- pos->res of mod
    pos_res_dict = protein_to_targets[protein_id]
    #print("prot:",protein_id)
    # get the alignment for that protein
    input_file = ALIGNMENT_FOLDER + protein_id + "_orthos.afa"

    if os.path.exists(input_file) and os.stat(input_file).st_size != 0:
        [amino_acids_at_target_positions,proteins_at_target_positions,amino_acids_at_target_positions_plus_one,all_original_protein_positions_for_alignment_pos] = read_alignment(input_file,protein_id,pos_res_dict)
        calculate_conservation_scores_per_ptm_site(protein_id,pos_res_dict,amino_acids_at_target_positions,proteins_at_target_positions,amino_acids_at_target_positions_plus_one,all_original_protein_positions_for_alignment_pos)
    else:
        print("No alignment file for protein_ID",protein_id,input_file)
        count_missing_alignnment+=1
        missing_pos.append(list(pos_res_dict.keys()))
    test_counter += 1
print(count_missing_alignnment, "missing alignments")
missing_pos_ls= [pos for position in missing_pos for pos in position] # make list of lists into 1 list
print(len(missing_pos_ls)) 
    #if test_counter == 5:
    #    exit()
