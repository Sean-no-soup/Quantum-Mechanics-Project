#Sean Heffley 5/8/23
#generate all Clebsch–Gordan coefficients for a given j1 and j2
from sympy.physics.quantum.cg import CG
from sympy import S,factorial,sqrt
import numpy as np
#S used to create number objects which can be fractions
#   rather than exclusively integers or floats:these are easier to read
#
#either the following function or sympy's included function, CG, may be used: check lines 82-83 with all the ##############s-----------------------------------

def calculateCoefficient(j1, m1, j2, m2, j3, m3):
    delta = 1 if (m3 == (m1 + m2)) else 0

    numerator = (2*j3+1) *      \
        factorial(j3+j1-j2) *   \
        factorial(j3-j1+j2) *   \
        factorial(j1+j2-j3) *   \
        factorial(j3+m3) *      \
        factorial(j3-m3) *      \
        factorial(j1-m1) *      \
        factorial(j1+m1) *      \
        factorial(j2-m2) *      \
        factorial(j2+m2)

    denominator = factorial(j1+j2+j3+1)
   
    summation = 0
    for k in range(0, j3 + m3 + 2):
        try:
            temp = (((-1)**k) /         \
                (factorial(k) *         \
                factorial(j1+j2-j3-k) * \
                factorial(j1-m1-k) *    \
                factorial(j2+m2-k) *    \
                factorial(j3-j2+m1+k)*  \
                factorial(j3-j1-m2+k)))
            summation += temp


        except:pass #factorials of negatives return an error (souldn't be counted so this is an easy way to exclude them)
         
    return delta * sqrt(numerator/denominator) * summation

#choose j1 and j2----------------------------------------------------------------------------------------------------------------------------------------------
def getHalfInt(wantedname):
    while True:
        print("------------\ninput {} as an integer '#' or fraction in format '#/2'".format(wantedname))
        j = input()
        if type(j) == str:
            j = j.split('/')
            if len(j) == 2 and j[0].isnumeric() and j[1] == "2":
                return S(j[0])/2
            if len(j) == 1 and j[0].isnumeric():
                return S(j[0])
        print('\tcannot fit to a valid input')

j1 = getHalfInt('j1')
j2 = getHalfInt('j2')

if j1 < j2: #swap if needed for convention
    print('swapping j1 and j2 for convention j1 > j2')
    temp = j1
    j1 = j2
    j2 = temp

#go through J, m1,m2, and M. create list of coefficients-------------------------------------------------------------------------------------------------------
coefficients = [] #[j1,m1,j2,m2,j3,m3,val] per entry

                     #S fractions don't work with range to create a for loop-> create values first then increment. exclude cgc(coefficients) that should = 0
j3 = j1 + j2         #max J
while j3 >= j1 - j2: #min J: convention j1 > j2 required

    m1 = -j1             #min m1
    while m1 <= j1 :     #max m1

        m2 = -j2             #min m2
        while m2 <= j2:      #max m2
            m3 = m1 + m2         #all other m3 -> cgc==0

            if np.absolute(m3) <= j3: #all others -> cgc==0

                #calculate coefficient#########################################################################################################################
                coefficient = calculateCoefficient(j1, m1, j2, m2, j3, m3) #function at very top of file
                ##coefficient = CG(j1, m1, j2, m2, j3, m3).doit()           #sympy function

                #print('\u27E8{} {} {} {} \u007c {} {}\u27E9'.format(j1, m1, j2, m2, j3, m3),end=' = ')#print what's going on, debug stuff
                #print(coefficient)
               
                if str(coefficient)[0] == "-": #formatting things to look nicer
                    coefficient =   "-sqrt({})".format(str(coefficient**2))
                else: coefficient = " sqrt({})".format(str(coefficient**2))

                coefficients += [[j1, m1, j2, m2, j3, m3,coefficient]] #put in a list for sorting later

            m2 += 1
        m1 += 1
    j3 -= 1

#trying to make nice display of coefficients
def sortCoef(n,sortlist):
    if n in ['j1','m1','j2','m2','J','M']:
        index = ['j1','m1','j2','m2','J','M'].index(n)
        return(sorted(sortlist,key = lambda x:x[index]))

sorted_list = sortCoef('m1',coefficients) #should start on the same m1(and pared m2) value for each chart
sorted_list = sortCoef('J',sorted_list)   #and follow a column of J
sorted_list = sortCoef('M',sorted_list)   #in groups of M. all compared to the wikipedia chart

maxlen = [0,0,0,0,0,0]  #find max character length of each j1,m1,j2,m2,j3,m3 (use a monospce font)
for i in sorted_list:
    for j in [0,1,2,3,4,5]:
        if len(str(i[j])) > maxlen[j]:
            maxlen[j] = len(str(i[j]))
            
running = 0
for j in sorted_list: #go through list of cgcs-----------------------------------------------------------------------------------------------------------------
    if running != j[5]:
        running = j[5]
        print('\n--  M={}  --'.format(running))#separate groups of different M

    #why is formatting so pain, this is disgusting code: print formatted <j1 m1 j2 m2 |j3 m3> = cgc with EVEN SPACINGS between each
    print("\u27E8{0:>{6}s} {1:>{7}s} {2:>{8}s} {3:>{9}s} \u007c {4:>{10}s} {5:>{11}s}\u27E9 = {12}\
          ".format(str(j[0]), str(j[1]), str(j[2]), str(j[3]), str(j[4]), str(j[5]),\
                   maxlen[0], maxlen[1], maxlen[2], maxlen[3], maxlen[4], maxlen[5],\
                    j[6]))










