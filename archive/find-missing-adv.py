# find missing adv reads in sequence from surface and lander data

import re

# find missing numbers in an array
def findMissingNumbers(numlist):
    misscount = 0
    for i in range( 1, len(numlist)):
        if numlist[i] < numlist[i-1]:
            misscount += 255 + numlist[i] - numlist[i-1]
        else:
            misscount += numlist[i] - 1 - numlist[i-1]
    return misscount

def missFrac(numlist):
    return findMissingNumbers(numlist) / len(numlist)

def readSequence(data_file):
    with open(data_file) as f:
        sequence = []
        for line in f:
            if line.startswith("D:"):
                sequence.append(int(re.search(r'\d+', line).group()))
    return sequence

### Tests
# Test data
test = [*range(0,256)]
testgood = test + test
testbad = testgood[:1] + testgood[2:]
testbad2 = testgood[:1] + testgood[3:]
testbad3 = test + test[1:]
testbad4 = test + test[2:]


# tests
findMissingNumbers(testgood) #expect 0
findMissingNumbers(testbad)  #expect 1
findMissingNumbers(testbad2) #expect 2
findMissingNumbers(testbad3) #expect 1
findMissingNumbers(testbad4) #expect 1

missFrac(testgood) # expect 0
missFrac(testbad)  # expect ~0.002

# check on lander data
file = '2023_05_05_13_13_09.txt'
fileseq = readSequence(file)
findMissingNumbers(fileseq)
missFrac(fileseq)

# check on surface data
file = '2023_05_05_12_13_22_sfc.txt'
fileseq = readSequence(file)
findMissingNumbers(fileseq)
missFrac(fileseq)
