'''
Creates a set of two toy sequences and aligns them using the Smith-Waterman
global alignment algorithm, similar to the one used in Blast. The concept at 
first relied on the idea that alignments could quickly be created with bitwise
operations. Unfortunately, in practice this implementation did not succeed, 
as it became more complex than a more simple script, as well as slower due to
the added complexity. Analysis of the script compared to a simple script,
including runtime, and the simple script I created are available on request. 

I do not plan on expanding further on this script.
'''


#Create sequences to align here, These can be changed but must only contain:
# A,C,G,T and no whitespace, symbols, etc.
seq1 = 'TGAAT'
seq2 = 'AAAGGAAT'
'''
BinConvert() takes an input sequence, and converts the sequence into a
binary number for fast processing. The binary number can be viewed as a
concatenated series of binary pairs, where '00' represents 'A', '01' represents
'C', '10' represents 'G', and '11' represents 'T'. The sequence is started with
a 1 so that any 'A' sequences at the start of the sequence are not lost.

Example: ATCGGAT = (header 1) 00 11 01 10 10 00 11 = 100110110100011
''' 

def BinConvert(seq):
    output = 0b1
    for i in seq:
        if i == 'A':
            output = output << 2
        elif i == 'C':
            output = output << 2
            output +=1
        elif i == 'G':
            output = output << 2
            output +=2
        elif i == 'T':
            output = output << 2
            output +=3
        else:
            return 'Input error at: ' + str(i)
    return output    


'''
BinKeyGenerator() creates a key for checking if any two sequences match.
The input is the 'pos'ition you want to start checking for matches,
and the 'n' number of matches you want. The position starts from the END
of the sequence. For example, to check for matches for the 'x's in the sequence
oooooxxooo, BinKeyGenerator(4,2) would produce the appropriate input.
'''
def BinKeyGenerator(pos, n): 
    key = 0b0
    for i in range(n):
        key = key << 2
        key += 3
    key = key << ((pos-1)*2)
    return key

'''
BinCompare() takes two sequences to be compared followed by the input for
BinKeyGenerator(), which are simply passed. The first step is creating a
new sequence called 'mismatches' which adds '1's wherever a mismatch occurs,
next posChecker checks if there are any mismatches that overlap the area we
are checking for, as indicated by 'pos' and 'n'. If there are no mismatches,
posChecker will equal 0.
'''
def BinPointCompare(seq1,seq2, pos, n):
    mismatches = seq1 ^ seq2 
    posChecker = mismatches & BinKeyGenerator(pos, n) 
    if posChecker == 0: #no mismatch at given position
        return True # a match
    else:
        return False #mismatch

#convert seqs to binary
bin_seq1 = BinConvert(seq1)
bin_seq2 = BinConvert(seq2)

#check for all of the last 100 seqs
#for i in range(1,100):
#    print BinPointCompare(bin_seq1,bin_seq2,i, 3),


#simply prints an array, useful for NeedlemanWunsch() troubleshooting
def ArrayPrinter(Arr):
    for i in range(len(Arr)):
        print Arr[i]

'''
BinSmith() takes two input sequences, the outputs of BinConvert
'''
def BenSmith(x_seq,y_seq):
    Arr = []
    row = [0] #starts with a single 0 for a header
    x_temp = x_seq
    y_temp = y_seq
    while x_temp != 1: #until sequence is empty, add a unit to 'row' and remove a unit from x_temp
        row.append(0) 
        x_temp = x_temp >> 2
    Arr.append(row) #starts with header row for array
    while y_temp != 1: #adds rows for each instance in y sequence
        Arr.append([0 for x in range(len(row))])
        y_temp = y_temp >> 2
    #
    #fill array
    for x in range(1,len(Arr[1])):
        for y in range(1,len(Arr)):
            x_temp = x_seq >> (len(Arr[1])-1 - x)*2
            y_temp = y_seq >> (len(Arr)-1 - y)*2
            if BinPointCompare(x_temp, y_temp, 1, 1): 
		Arr[y][x] = Arr[y-1][x-1] + 2
            else:
                if (Arr[y-1][x-1] >= Arr[y][x-1]) and (Arr[y-1][x-1] >= Arr[y-1][x]):
                    Arr[y][x] = Arr[y-1][x-1] - 1
                elif Arr[y-1][x] >= Arr[y][x-1]:
                    Arr[y][x] = Arr[y-1][x] - 1
                else:
                    Arr[y][x] = Arr[y][x-1] - 1
            if Arr[y][x] < 0:
                Arr[y][x] = 0
    ArrayPrinter(Arr)
    return BinTrace(Arr, x_seq, y_seq)

def BinTrace(Arr, x_seq, y_seq):
    #find max value in array
    max = 0
    for y in range(len(Arr)):
        for x in range(len(Arr[1])):
            if Arr[y][x] >= max:
                max = Arr[y][x]
                x_pos = x
                y_pos = y
    output_y = ''
    output_x = ''
    suffix_key = 3
    y_jump_len = len(Arr) - y_pos -1
    x_jump_len = len(Arr[1]) - x_pos -1
    y_temp = y_seq >> y_jump_len
    x_temp = x_seq >> x_jump_len
    while y_seq != 1: 
        while x_seq != 1:
            if  Arr[y][x] == 0:
                return output_y,output_x
            elif (Arr[y-1][x-1] >= Arr[y-1][x]) and (Arr[y-1][x-1] >= Arr[y][x-1]):
                y_single = y_seq & suffix_key
                x_single = x_seq & suffix_key
                output_y = BinReturn(y_single) + output_y
                output_x = BinReturn(x_single) + output_x
                y_seq = y_seq >> 2
                x_seq = x_seq >> 2
                x -=1
                y -=1
            elif Arr[y-1][x] >= Arr[y][x-1]:
                y_single = y_seq & suffix_key
                output_y = BinReturn(y_single) + output_y
                output_x = '-' + output_x
                y_seq = y_seq >> 2
                y-=1
            elif Arr[y][x-1] > Arr[y-1][x]:
                x_single = x_seq & suffix_key
                output_x = BinReturn(x_single) + output_x
                output_y = '-' + output_y
                x_seq = x_seq >> 2
                x-=1
            else:
                print 'Error: misdirected traceback'
    return output_x,output_y
def BinReturn(seq):
    if seq == 0:
        return 'A'
    elif seq == 1:
        return 'C'
    elif seq == 2:
        return 'G'
    elif seq == 3:
        return 'T'
    else:
        print 'Error: improper input into BinReturn'



print BenSmith(bin_seq1,bin_seq2)

