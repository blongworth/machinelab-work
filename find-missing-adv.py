# find missing adv reads in sequence from surface and lander data

import re

test = [*range(0,256)]
testgood = test + test

testbad = testgood[:1] + testgood[2:]

# find missing numbers in an array
def findMissingNumbers(numlist):
    misscount = 0
    for i in range( 1, len(numlist)):
        if numlist[i] == 0:
            misscount += 255 - numlist[i-1]
        else:
            misscount += numlist[i] - 1 - numlist[i-1]
    return misscount

def missFrac(numlist):
    return findMissingNumbers(numlist) / len(numlist)

findMissingNumbers(testgood)
findMissingNumbers(testbad)

missFrac(testgood)
missFrac(testbad)


                sequence.append(int(re.search(r'\d+', line).group()))
    return sequence

file = '2023_05_05_13_13_09.txt'
fileseq = readSequence(file)
findMissingNumbers(fileseq)
missFrac(fileseq)

file = '2023_05_05_12_13_22_sfc.txt'
fileseq = readSequence(file)
findMissingNumbers(fileseq)
missFrac(fileseq)
