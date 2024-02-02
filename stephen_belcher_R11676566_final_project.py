#! /user/bin/env python3
"""
==================================================================================
Title: stephen_belcher_R11676566_final_project.py
Author: stbelche (R#11676566)
Date: 12/1/2023
Version: 1.0
Usage: python python final.py -i SequenceData.fna -o sequences -s BLOSUM62.mtx
Notes: this program is the final project, a fully functioning bioinformatics pipeline
Python Version: 3.11 (i think)
==================================================================================
"""

#===============================================================
#store_to_dict
def store_dict(blosum_mat, score_file):
#function stores the scoring matrix into a 2-d dictionary to
#   ensure constant time lookup
#function also stores the minimum and maximum alignment scores
#===============================================================
    first_line = score_file.readline().split()  #first line of matrix is the metadata
    reference_array = []    #stores the order of the proteins/nucleotides
    smin = 999 #minimum alignment score
    smax = 0 #max alignment score
    for index in first_line:
     blosum_mat[index] = {}  #add it to the dictionary
     reference_array.append(index)
    for lines in score_file:
        x = 1   #keeps track of the column
        lines = lines.split()
        y_value = lines[0]  #first column in each row is metadata
        for index in range(len(lines)-1):
            score = int(lines[x])
            blosum_mat[reference_array[x-1]][y_value] = score    #store the data
            # if we find a new min or max alignment score, save it:
            if (score > smax):
                smax = score
            if (score < smin) and ( lines[0] != '-' and reference_array[x-1][0] != '-' ):
                if(lines[0] != 'X' and lines[x-1][0] != 'X'):
                    smin = score
            x += 1
    return (blosum_mat, smin, smax)
#   end of function
#========================================================================
#==================================================================
# store_old_sequence
#inputs: the sequence data (string),
#        the count of how many bases there are (int)
#        the list of sequences
# function: the sequence is stored into the list of sequences (LOS)
# note: the function sorts the list as it inserts the sequence
def store_old_sequence(sequence, LOS):
    stored = False
    index = 0
    while not stored:
      if (len(LOS) == 0) or (index > len(LOS)): #if we have an empty list
        LOS.append(sequence)
        stored = True
      elif(len(sequence) > len(LOS[index])): #if this seq is longer than indexed seq
        LOS.insert(index, sequence)
        stored = True
      else:
        index += 1
#end of function
#==============================================================================
#   sort_sequences
#input: a .fna file containing sequences
#output: A sorted list of sequences
#function sorts the sequences from the input and writes the sorted sequence to
#   A list of sorted sequences
def sort_sequences(input_file):
    #create the list of lists and the var holding the sequence and number of chars
    LOS = []  #list of sequences
    sequence = []
    #algorithm for processing sequences
    for line in input_file:
      if line[0] == '>':
          store_old_sequence(sequence, LOS)
          sequence = [] #sequence contains ['header', 'base 1', 'base 2',...]
          sequence.append(line.strip())
      else:
          for each in line:
            if(each != '\n'):
                sequence.append(each)
    LOS.pop() #im getting an empty list at the end of my LOS
    return(LOS)
#end of function
#===============================================================================
# write_to_file
#input: list of sequences and an output file
#function takes in the input file and writes to the output in the proper format
#   in the order of the LOS
#output: a properly written output file
def write_to_file(LOS, output_file):
    for each in LOS:
        output_file.write(each[0])
        output_file.write('\n')
        for index in range(1,len(each)):
            output_file.write(each[index])
        output_file.write('\n')
#end of function

#========================================================================
#needleman-wunsch
def needleman_wunsch(seq_1, seq_2, stage):
#function performs needleman_wunsch
#output changes depending on stage
#=========================================================================
    #we dont want to permanently change the sequences
    seq_1 = copy.deepcopy(seq_1)
    seq_2 = copy.deepcopy(seq_2)
    #first, we must initialize the scoring matrix and traceback matrix:
    score_mat = [['orig', 'blank'], ['blank', 0] ]
    traceback_mat = [['orig', 'blank'], ['blank', 'DONE'] ]
    d = 1
    for each in seq_1:
        score_mat[0].append(each) #adding sequence 1 to the top row
        traceback_mat[0].append(each) #apply same thing to traceback matrix
        score_mat[1].append(d * gap_penalty) #adding the multiplier
        d += 1
        traceback_mat[1].append('L') #top row for traceback matrix is always Left
    d = 1
    for each in seq_2:
        score_mat.append([each, d * gap_penalty]) #adding sequence 2 to the scoring matrix
        traceback_mat.append([each, 'U']) #adding sequence 2 to the traceback matrix
        d += 1
    #now, we must calculate the scores and assign the correct traceback designation
    for row in range(2, len(seq_2) + 2):
        for col in range(2, len(seq_1) + 2):
            #calculate each potential score for diag, left, and up
            diag_score = score_mat[row-1][col-1] + sub_mat[ (score_mat[0][col]) ][ (score_mat[row][0])  ]
            left_score = score_mat[row][col-1] + gap_penalty
            up_score = score_mat[row-1][col] + gap_penalty
            #compare them
            if diag_score >= left_score and diag_score >= up_score: #diag has priority 1
                cell_score = diag_score
                designation = 'D'
            elif left_score >= up_score: #left has priority 2
                cell_score = left_score
                designation = 'L'
            else:   #up has lowest priority
                cell_score = up_score
                designation = 'U'
            #now we record the score and character into their respective matrix
            score_mat[row].append(cell_score)
            traceback_mat[row].append(designation)

    #now we go through the traceback matrix and align the sequences

    #we start at the bottom right, so we first must create the index for it
    row = len(traceback_mat) - 1
    col = len(traceback_mat[0]) -1
    designation = traceback_mat[row][col]
    #now we go up and left until we reach the cell that says 'DONE
    while designation != 'DONE':
        if designation == 'D':
            row -= 1
            col -= 1
        elif designation == 'U': #function to gap up
            row -= 1
            seq_1.insert(col-1, '-')
        else:   #function to gap left
            col -= 1
            seq_2.insert(row-1, '-')
        #print(designation)
        designation = traceback_mat[row][col]

    final_score = score_mat[-1][-1] #i love negative indexing

    if stage == 1: #stage 1 = aligning the sequences
        return(seq_1, seq_2)
    if stage == 2: #stage 2 = for the distance matrix (task 4)
        return(len(seq_2), final_score)
    if stage == 3: #to test MSA's (task 4)
        return final_score
# end of function
#===================================================================

#===============================================================================
# percent_identity
#inputs: two sequences, aligned via the needleman-wunsch algorithm
#output: a floating point number
#function takes the two sequences and calculates the percentage of identical
#   bases
#==============================================================================
#percent_identity
def percent_identity(sequence_a, sequence_b):
    identical_count = 0
    for i in range(len(sequence_a)):
        if sequence_a[i] == sequence_b[i]:
            identical_count += 1
    return (identical_count / len(sequence_a))
#end of function
#===============================================================================
#smith-waterman
def smith_waterman(sequence_1, sequence_2):
    #first, we must initialize the scoring matrix and traceback matrix:
    score_mat = [['orig', 'blank'], ['blank', 0] ]
    traceback_mat = [['orig', 'blank'], ['blank', 'X'] ]
    for each in sequence_1:
        score_mat[0].append(each) #adding sequence 1 to the top row
        traceback_mat[0].append(each) #apply same thing to traceback matrix
        score_mat[1].append(0) #adding the edges of the matrix
        traceback_mat[1].append('X') #top row for traceback matrix is always Le
    for each in sequence_2:
        score_mat.append([each, 0]) #adding sequence 2 to the scoring matrix
        traceback_mat.append([each, 'X']) #adding sequence 2 to the traceback matrix

    #now, we must calculate the scores and assign the correct traceback designation
    highest_score = 0
    x_coord = 2
    y_coord = 2
    for row in range(2, len(sequence_2) + 2):
        for col in range(2, len(sequence_1) + 2):
            #calculate each potential score for diag, left, and up
            diag_score = score_mat[row-1][col-1] + sub_mat[ (score_mat[0][col]) ][ (score_mat[row][0])  ]
            left_score = score_mat[row][col-1] + gap_penalty
            up_score = score_mat[row-1][col] + gap_penalty
            #compare them
            if diag_score >= left_score and diag_score >= up_score: #diag has priority 1
                cell_score = diag_score
                designation = 'D'
            elif left_score >= up_score: #left has priority 2
                cell_score = left_score
                designation = 'L'
            else:   #up has lowest priority
                cell_score = up_score
                designation = 'U'
            #now we compare it to zero
            if cell_score <= 0:
                cell_score = 0 #nothing can be less than zero
                designation = 'X'
            #now we record the score and character into their respective matrix
            score_mat[row].append(cell_score)
            traceback_mat[row].append(designation)
            #now we check if the score is greater than the highest score
            if cell_score > highest_score:
                x_coord = row
                y_coord = col
                highest_score = cell_score
    #beginning traceback...
    #we start at the highest score, with an empty pair of sequences
    designation = traceback_mat[x_coord][y_coord]
    seq_1 = ""
    seq_2 = ""
    #now we go up and left until we reach the cell that says 'X'
    while designation != 'X':
        if designation == 'D':
            seq_1 = traceback_mat[0][y_coord] + seq_1
            seq_2 = traceback_mat[x_coord][0] + seq_2
            x_coord -= 1
            y_coord -= 1
        elif designation == 'U': #function to gap up
            seq_1 = '-' + seq_1
            seq_2 = traceback_mat[x_coord][0] + seq_2
            x_coord -= 1
        else:   #function to gap left
            seq_1 = traceback_mat[0][y_coord] + seq_1
            seq_2 = '-' + seq_2
            y_coord -= 1
        #print(designation)
        designation = traceback_mat[x_coord][y_coord]
    return(highest_score, len(seq_1))
#================================================================================
#parse_string
def parse_string(string):
#function takes in the string and separates the sequence headers
#stores the headers into a list and returns the list
#===================================================================
    msa = []
    #print('parsing the string ', string)
    seq_to_find = '>'
    for char in string[1:]:
        if (char != '>'):
            seq_to_find += char
        else:
            msa.append(seq_to_find)
            seq_to_find = '>'
    msa.append(seq_to_find)
    return msa
#end of function
#========================================================================
#calculate_distance
def calculate_distance(seq, msa):
#function takes a sequence and a msa or a pair of msa's
#   returns the distance
#========================================================================
    #print("now i must calculate ", seq, ' with ', msa)
    dist = 0
    msa_numbers = parse_string(msa)
    seq_numbers = parse_string(seq)
    for each in seq_numbers:
        index = sequence_dict[each]
        for every in msa_numbers:
            jndex = sequence_dict[every]
            dist += distance_matrix[index][jndex]
    dist = dist / (len(msa_numbers) * len(seq_numbers) )
    return dist
#end of function
#=========================================================================
#                                MAIN PROGRAM
#===============================================================================
#import functions
import argparse
import copy
import math
#------------------------------------------------------------------------------
#                       TASK 1
#GOAL: take in the .fna file, sort it, then output the result to another file
#------------------------------------------------------------------------------
print("Final Project :: R11676566")
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', required=True)
parser.add_argument('-o', '--output_file',required=True)
parser.add_argument('-s', '--score_matrix', required=True)
#parse the arguments
args = parser.parse_args()
#assign variables to inputs
infile = open(args.input_file, 'r') #input file contains sequences
outfile_string = args.output_file + '.sorted.fna' #output file to write to
outfile = open(outfile_string, 'w')
score_file = open(args.score_matrix, 'r') #scoring matrix

sub_mat = {}    #substitution matrix
(sub_mat, smin, smax) = store_dict(sub_mat, score_file)
score_file.close() #we don't need it anymore
gap_penalty = int(sub_mat['-']['A']) #defines the gap penalty

LOS = sort_sequences(infile) #generate LOS and sort sequences
write_to_file(LOS, outfile) #writes LOS to output file
outfile.close()
'''
for each in LOS:
    print(each)
'''
#-------------------------------------------------------------------------------
#                           TASK 2
#GOAL: organize the sequences into clusters
#--------------------------------------------------------------------------------
#clusters is a list of lists containing the indexes of the sequences from LOS
clusters = [ [0] ] #initialize the list, first cluster contains LOS[0]
cluster_reps = [LOS[0][1:]] #save the representative sequence for each cluster

for i in range(1, len(LOS)): #go through each sequence
    seq_to_examine = LOS[i][1:]
    stored = False
    for j in range(0, len(cluster_reps)): #go through each cluster for the sequence
        representative = cluster_reps[j]
        seq_x = copy.deepcopy(seq_to_examine) #use deepcopy to avoid changing the orignal seq
        seq_y = copy.deepcopy(representative)
        (seq_x, seq_y) = needleman_wunsch(seq_x, seq_y, stage= 1)
        ident = percent_identity(seq_x, seq_y) #calculate percent identicality
        if ident >= .97:
            clusters[j].append(i)
            stored = True
            break
    #if we haven't stopped already, then this sequence is part of a new cluster
    if not stored:
        clusters.append([i])
        cluster_reps.append(LOS[i][1:])
'''
for each in clusters:
    print(each)
'''
#now we must output the representative of each cluster
outfile_string = args.output_file + '.clustered.fna'
outfile = open(outfile_string,'w')
for each in clusters:
    what_to_write = LOS[each[0]][0] + ';size='+ str(len(each)) #add ';size =' to str
    LOS[each[0]][0] = what_to_write
    outfile.write(what_to_write)
    outfile.write('\n')
    for i in range(1, len(LOS[each[0]])):
        outfile.write(LOS[each[0]][i])
    outfile.write('\n')
outfile.close()
#----------------------------------------------------------------
#                       TASK 3
#GOAL: filter for chimeras
#----------------------------------------------------------------
chimeric_string = args.output_file + '.chimeric.fna'
nonchimeric_string = args.output_file + '.nonchimeric.fna'
chimera_file = open(chimeric_string, 'w')
nonchimera_file = open(nonchimeric_string, 'w')
sequences = [] #sequences is used in task 4, we create the data structure here
for i in range(0, len(cluster_reps)):
    similarity_count = 0 #if similarities is more than 1, its a chimera
    for j in range(0, len(cluster_reps)):
        seq_a = copy.deepcopy(cluster_reps[i])
        if(i != j):
            seq_b = copy.deepcopy(cluster_reps[j])
            (score, length) = smith_waterman(seq_a, seq_b)
            if score >= 50 and length >= 60:
                if len(clusters[i]) < len(clusters[j]):
                    similarity_count += 1
    original_index = clusters[i][0]
    if similarity_count >= 2: #this is likely a chimera
        #print('chimera detected, header = ', LOS[original_index][0])
        chimera_file.write(LOS[original_index][0] + '\n') #write to file
        for each in cluster_reps[i]:
            chimera_file.write(each)
    else: #this is likely not a chimera
        nonchimera_file.write(LOS[original_index][0] + '\n')
        sequences.append([LOS[original_index][0]]) #add the header to sequences
        sequences[-1].append(cluster_reps[i]) #add the sequence to sequences
        for each in cluster_reps[i]:
            nonchimera_file.write(each)
        nonchimera_file.write('\n')
chimera_file.close()
nonchimera_file.close()
'''
#-------------------------------------------------------------------------------------
#                                   TASK 4
#GOAL: implement fitch-margoliash MSA on the nonchimeric cluster representatives
#-------------------------------------------------------------------------------------
#FIRST WE SET STUFF UP

#we are now working with the 'sequences' matrix
'''
outfile_string = args.output_file + '.msa.fna'
msa_outfile = open(outfile_string, 'w')
sequence_dict = {} #dictionary matches the header with that sequences index
for i in range(len(sequences)):
    sequence_dict[sequences[i][0]] = i + 1 #this is important for the distance matrix
num_sequences = len(sequences)
#print(sequence_dict)
#print(num_sequences)
alignment_queue = []
distance_matrix = [['null']]

#NOW WE INITIALIZE THE DISTANCE MATRIX
for each in sequences:
    distance_matrix.append([each[0]])
    distance_matrix[0].append(each[0])
#perform needleman_wunsch on every combination of sequences and calculate their distance
for i in range(1, len(sequences) + 1 ):
    for j in range(1, len(sequences) + 1):
        if(i == j): #if we are comparing a sequence to itself
            distance_matrix[i].append(0) #distance from a sequence to itself is 0
        else:
            #first we get the raw alignment score and the length of the alignment
            (seq_len, score) = needleman_wunsch(sequences[i-1][1], sequences[j-1][1], 2)
            #now we normalize the score
            score = (score - (smin * seq_len) ) / ( (smax * seq_len) - (smin * seq_len) )
            #now we convert the normalized score to distance
            score = math.log(score)
            score = -score
            distance_matrix[i].append(score)
'''
for each in distance_matrix:
    print(each)
print('\n\n')
'''

#NOW WE CONSTRUCT THE GUIDE TREE
#we will not touch the distance matrix, we will only manipulate the current dist matrix
cur_dist_mat = copy.deepcopy(distance_matrix) #current distance matrix
while len(cur_dist_mat) > 3: #until we only have 2 things left to combine
    #now we find the minimum alignment score
    min_dist = 100
    i_index = 0 #record the index of i
    j_index = 0 #record the index of j
    for i in range(1, len(cur_dist_mat)):
        for j in range(1, i + 1):
            if(cur_dist_mat[i][j] > 0 and cur_dist_mat[i][j] < min_dist):
                min_dist = cur_dist_mat[i][j]
                row_index = i
                col_index = j

    #print(distance_matrix[0][row_index], distance_matrix[0][col_index], min_dist)
    alpha = min_dist
    beta = 0
    theta = 0
    count = 0
    #calculate beta and theta
    for i in range(1, len(cur_dist_mat)):
        if(i != row_index and i != col_index):
            beta += cur_dist_mat[i][row_index]
            theta += cur_dist_mat[i][col_index]
            count += 1
    beta = beta / count
    theta = theta / count

    #print('alpha = ', alpha, 'beta = ', beta, 'theta = ', theta)

    b = (alpha + beta - theta)/2
    a = (alpha + theta - beta)/2

    #print('b = ', b, 'a = ', a)
    if (a > b):
        alignment_order = cur_dist_mat[0][row_index] + cur_dist_mat[0][col_index]
        alignment_queue.append([ cur_dist_mat[0][row_index], cur_dist_mat[0][col_index]])
    else:
        alignment_order = cur_dist_mat[0][col_index] + cur_dist_mat[0][row_index]
        alignment_queue.append([ cur_dist_mat[0][col_index],cur_dist_mat[0][row_index] ])
    #print(alignment_order)

    #now we update the matrix
    cur_dist_mat[0].pop(row_index) #we no longer use this
    cur_dist_mat[0].pop(col_index) #we no longer use this
    cur_dist_mat[row_index][row_index] = -1 #indicate that its not used
    cur_dist_mat[col_index][col_index] = -1
    cur_dist_mat[0].append(alignment_order) #this is now in play


    msa_sequences = parse_string(alignment_order)
   

    #create a new distance matrix
    next_dist_mat = copy.deepcopy([cur_dist_mat[0]]) #create the header row
    for i in range(1, len(next_dist_mat[0])):
        next_dist_mat.append([next_dist_mat[0][i]]) #fill in the metadata
        for j in range(1, len(next_dist_mat[0]) ):
            if i == j:
                next_dist_mat[i].append(0) #distance from x to x is 0
            else:
                x = next_dist_mat[0][j]
                y = next_dist_mat[i][0]
                if(x in sequence_dict and y in sequence_dict):
                    next_dist_mat[i].append(distance_matrix[sequence_dict[x]][sequence_dict[y]])
                else:
                    next_dist_mat[i].append(calculate_distance(x,y) )

    cur_dist_mat = copy.deepcopy(next_dist_mat)
alignment_order = cur_dist_mat[0][1] + cur_dist_mat[0][2]
alignment_queue.append([cur_dist_mat[0][1], cur_dist_mat[0][2] ])
'''
for each in alignment_queue:
    print(each)
'''
#NOW WE PERFORM THE MSA
for i in range(len(alignment_queue)):
    #print('iteration ', i)
    s1 = parse_string(alignment_queue[i][0])
    s2 = parse_string(alignment_queue[i][1])
    #print('s1 = ', s1)
    #print('s2 = ', s2)
    if len(s1) == 1: #if 1 = seq
        i1 = sequence_dict[s1[0]] - 1
        if len(s2) == 1: #if 1 = seq and 2 = seq
            i2 = sequence_dict[s2[0]] - 1
            (sequences[i1][1], sequences[i2][1]) = needleman_wunsch(sequences[i1][1], sequences[i2][1], 1)
        else: #if 1 = seq and 2 = msa
            max_score = 0
            winner = 0
            for j in range(len(s2)): #in this loop, we choose the best seq from the msa
                i2 = sequence_dict[s2[j]] - 1
                score = needleman_wunsch(sequences[i1][1], sequences[i2][1], 3)
                if score > max_score:
                    max_score = score
                    winner = j
            #print('winner = ', winner)
            i2 = sequence_dict[s2[winner]] - 1
            (sequences[i1][1], sequences[i2][1]) = needleman_wunsch(sequences[i1][1], sequences[i2][1], 1)
            #now we must copy all gaps from the chosen msa
            guide = sequences[sequence_dict[s2[winner]] - 1][1] #the reference
            for k in range(0,len(s2)):
                if k != winner:
                    to_change = sequences[sequence_dict[s2[k]] - 1 ][1] # what youre changing
                    for base in range(len(guide)):
                        if guide[base] == '-':
                            to_change.insert(base, '-')
    else: #1 = msa
        if len(s2) == 1: #1 = msa and 2 = seq
            i2 = sequence_dict[s2[0]] - 1
            max_score = 0
            winner = 0
            for j in range(len(s1)): #in this loop, we choose the best seq from the msa
                i1 = sequence_dict[s1[j]] - 1
                score = needleman_wunsch(sequences[i1][1], sequences[i2][1], 3)
                if score > max_score:
                    max_score = score
                    winner = j
            #print('winner = ', winner)
            i1 = sequence_dict[s1[winner]] - 1
            (sequences[i1][1], sequences[i2][1]) = needleman_wunsch(sequences[i1][1], sequences[i2][1], 1)
            #now we must copy all gaps from the chosen msa
            guide = sequences[sequence_dict[s1[winner]] - 1][1] #the reference
            for k in range(len(s1)):
                if k != winner:
                    to_change = sequences[sequence_dict[s1[k]] - 1 ][1] # what youre changing
                    for base in range(len(guide)):
                        if guide[base] == '-':
                            to_change.insert(base, '-')
        else: #1 = msa and 2 = msa
            max_score = 0
            winner1 = 0
            winner2 = 0
            for j in range(len(s1)): #in this loop, we choose the best seq from the msa
                i1 = sequence_dict[s1[j]] - 1
                for k in range(len(s2)):
                    i2 = sequence_dict[s2[k]] -1
                    score = needleman_wunsch(sequences[i1][1], sequences[i2][1], 3)
                    if score > max_score:
                        max_score = score
                        winner1 = j
                        winner2 = k
            i1 = sequence_dict[s1[winner1]] - 1
            i2 = sequence_dict[s2[winner2]] - 1
            (sequences[i1][1], sequences[i2][1]) = needleman_wunsch(sequences[i1][1], sequences[i2][1], 1)
            guide = sequences[sequence_dict[s1[winner1]] - 1][1] #the reference
            for k in range(len(s1)):
                if k != winner1:
                    to_change = sequences[sequence_dict[s1[k]] - 1 ][1] # what youre changing
                    for base in range(len(guide)):
                        if guide[base] == '-':
                            to_change.insert(base, '-')
            guide = sequences[sequence_dict[s2[winner2]] - 1][1] #the reference
            for k in range(0,len(s2)):
                if k != winner2:
                    to_change = sequences[sequence_dict[s2[k]] - 1 ][1] # what youre changing
                    for base in range(len(guide)):
                        if guide[base] == '-':
                            to_change.insert(base, '-')
    #now we convert each gap to an 'X' so that we know in the future
    for each in sequences:
        for l in range(len(each[1])):
            if each[1][l] == '-':
                each[1][l] = 'X'
        #print(each)
    #print('\n')
    #we will convert them back at the end
#now we switch it back
for each in sequences:
        for l in range(len(each[1])):
            if each[1][l] == 'X':
                each[1][l] = '-'

#NOW WE CALCULATE SUM OF PAIRS
sum_of_pairs = 0
for c in range(len(sequences[0][1])): #for each column
    column_total = 0 #sum for each column
    for i in range(0, num_sequences -1):
        for j in range(i + 1, num_sequences):
            column_total += sub_mat[sequences[i][1][c]][ sequences[j][1][c]]
    sum_of_pairs += column_total
#print(sum_of_pairs)


for each in sequences:
    i = 0
    each[0] += "; score =" + str(sum_of_pairs)
    msa_outfile.write(each[0])
    msa_outfile.write('\n')
    for every in each[1]:
        msa_outfile.write(every)
    msa_outfile.write('\n')


#I think this is the longest program I have ever written