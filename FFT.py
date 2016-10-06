import time
import math
import cmath
import random

'''  In this program, we represent a polynomial 
     A(x) = a0 + a1x + a2x^2 + ... + a_{k-1}x^{k-1}$ as a Python list 
     A = [a0, a1, a2, ..., a_{n-1}].  Similarly, we represent a sequence
     of points where it is evaluated as as a python list [x0, x1, .., x_{n-1}],
     and the value of the polynomial at those points with a python list
     list [y0, y1, ..., y_{n-1}].  These are essentially vectors,
     as described in the readings, and I will sometimes refer to them
     as vectors.
'''

''' A is a vector of any size >= 1, and inverse is False if you want the
    DFT, and True if you want the inverse DFT.  To get the length of A
    to be a power of 2, we can pad it with zeros.  The inverse is just
    the same as the forward one except that it winds around the circle
    in the opposite way (in other words it uses omega_n^{-1} in place
    of omega_n), and it must be divided by n'''
def FFT(A, inverse):
    l = nextPower2(len(A))         # smallest power of two >= this degree bound
    A = A + [0] * (l - len(A))
    if inverse == False:
        omega_n = cmath.exp(complex(0, 2*cmath.pi/len(A)))
    else:
        omega_n = cmath.exp(complex(0, -2*cmath.pi/len(A)))
    Y = FFTAux(A, omega_n)
    if inverse == False:
        return Y
    else:
        return [Y[i] / len(Y) for i in range(len(Y))]

''' You have to fill out the following method.  A precondition is that the 
    length n of A is a power of two and omega_n is the principal nth root of 
    1  if you want the DFT, and its multiplicative inverse if you want the
    inverse DFT.  

    It would not be interesting to adhere too closely to the book's pseudocode.
    Instead, when you get YE, the DFT of the even subscripted elements
    of A, double it up by letting YE = YE + YE.  Do a similar thing
    for the DFT of the odd coefficients, and then rotate them by
    the appropriate amounts, as we discussed in class.  Return the
    vector sum of AE and AO. 

    Here's a useful python trick:  A[0:len(A):2] gives the elements in
    even positions of A and A[1:len(A):2] gives the elements in 
    odd positions.  Try it out in the interpreter on a list.  Try changing
    the value after the third colon to see what role it has.  '''
def FFTAux(A, omega_n):
    if (len(A) == 1):
        # Basecase, return A
        return A
    AE = A[0:len(A):2]
    AO = A[1:len(A):2]
    # List size / 2, double up omega_n
    Even_result = FFTAux(AE, omega_n * omega_n)
    Odd_result = FFTAux(AO, omega_n * omega_n)
    # Double up even result and odd result
    Even_result = Even_result * 2
    Odd_result = Odd_result * 2
    # Start to rotate each entry of odd_result to make them right
    current_omega = omega_n
    for i in range(len(Odd_result)-1):
        Odd_result[i+1] *= current_omega
        current_omega *= omega_n
    # Add up even result and odd result
    Return_result = [0] * len(Even_result)
    for i in range(len(Even_result)):
        Return_result[i] = Even_result[i] + Odd_result[i]
    return Return_result

'''  This is to work just like FFT, except that it pads A so that its
     length is the next power of three, so that it can use a three-way
     recursion instead of a two-way one. '''
def FFT3(A, inverse):
    l = nextPower3(len(A))  # smallest power of three >= this degree bound
    A = A + [0] * (l - len(A))
    if inverse == False:
        omega_n = cmath.exp(complex(0, 2 * cmath.pi / len(A)))
    else:
        omega_n = cmath.exp(complex(0, -2 * cmath.pi / len(A)))
    Y = FFT3Aux(A, omega_n)
    if inverse == False:
        return Y
    else:
        return [Y[i] / len(Y) for i in range(len(Y))]

'''  This has to work like FFTAux, except that the length of A
     is a power of three, and it makes three recursive calls of size n/3
     instead of two recursive calls of size n/2.  The standard FFT is
     based on the expression A(x) = A^e(x) + xA^o(x), where A^e and
     A^o are the even- and odd-indexed coefficients.  This one
     should be based on the expression A(x) = A^0(x) + xA^1(x) + x^2A^2(x),
     where A^0 is the coefficients whose subscripts are divisible
     by three, A^1 is the coefficients whose subscripts are 1 (mod 3),
     and A^2 is the coefficients whose subscripts are 2 (mod 3).  It
     should triple up the results of the recursive calls using
     concatenation, rather than doubling them up, as you did in
     FFTAux.  Add any additional functions to the file that you need to make
     it work.
'''
def FFT3Aux(A, omega_n):
    if (len(A) == 1):
        # Basecase, return A
        return A
    A0 = A[0:len(A):3]
    A1 = A[1:len(A):3]
    A2 = A[2:len(A):3]
    # List size / 3, triple up omega_n
    A0_result = FFT3Aux(A0, omega_n * omega_n * omega_n)
    A1_result = FFT3Aux(A1, omega_n * omega_n * omega_n)
    A2_result = FFT3Aux(A2, omega_n * omega_n * omega_n)
    # Triple up even result and odd result
    A0_result = A0_result * 3
    A1_result = A1_result * 3
    A2_result = A2_result * 3
    # Start to rotate each entry of A1_result to make them right
    current_omega = omega_n
    for i in range(len(A1_result)-1):
        A1_result[i+1] *= current_omega
        current_omega *= omega_n
    # Start to rotate each entry of A2_result to make them right
    current_omega = omega_n * omega_n
    for i in range(len(A2_result)-1):
        A2_result[i+1] *= current_omega
        current_omega *= (omega_n * omega_n)
    # Add up even result and odd result
    Return_result = [0] * len(A0_result)
    for i in range(len(A0_result)):
        Return_result[i] = A0_result[i] + A1_result[i] + A2_result[i]
    return Return_result

def logCeil(n):
    result = 0
    while (n > 1):
        result += 1
        n = n / 2 + n % 2
    return result

def logCeil3(n):
    result = 0
    while (n > 1):
        result += 1
        n = n / 3 + n % 3
    return result

def twoToThe(i):
    result = 1
    while (i > 0):
        result = result * 2
        i -= 1
    return result

def threeToThe(i):
    result = 1
    while (i > 0):
        result = result * 3
        i -= 1
    return result

''' Naive multiplication algorithm for two polynomials.  Implement the
    O(n^2) algorithm for multiplying polynomials A1 and A2.  '''
def naiveMult(A1, A2):
    Result = [0] * (len(A1) + len(A2))
    for i in range(len(A1)):
        for j in range(len(A2)):
            Result[i+j] += A1[i] * A2[j]
    return Result

def nextPower2(i):
    return twoToThe(logCeil(i))

def nextPower3(i):
    return threeToThe(logCeil3(i))

'''  Fill in the rest of this, so that it computes the product of
     two polynomials that have integer coefficients in O(n log n) time.  
     Use the approach suggested by the book chapter.

     The output list should match the one given by naiveMult; to
     achieve this, you will need to round the real part of each coefficient
     of the inverse DFT, since it will have small roundoff errors.   '''
def polyMult(A1, A2):
    l = len(A1) + len(A2) - 1  # number of terms in product polynomial
    l2 = nextPower2(l)         # smallest power of two >= this
    # Populate A1 and A2 to length of l2
    A1 = A1 + [0] * (l2 - len(A1))
    A2 = A2 + [0] * (l2 - len(A2))
    A1_Point_Value = FFT(A1, False)
    A2_Point_Value = FFT(A2, False)
    Result_Point_Value = [0] * l2
    for i in range(l2):
        Result_Point_Value[i] = A1_Point_Value[i] * A2_Point_Value[i]
    # Now do inverse FFT
    Result = FFT(Result_Point_Value, True)
    for i in range(len(Result)):
        Result[i] = int(Result[i].real + 0.00001)
    return Result

def polyMult3(A1, A2):
    l = len(A1) + len(A2) - 1  # number of terms in product polynomial
    l2 = nextPower3(l)         # smallest power of two >= this
    # Populate A1 and A2 to length of l2
    A1 = A1 + [0] * (l2 - len(A1))
    A2 = A2 + [0] * (l2 - len(A2))
    A1_Point_Value = FFT3(A1, False)
    A2_Point_Value = FFT3(A2, False)
    Result_Point_Value = [0] * l2
    for i in range(l2):
        Result_Point_Value[i] = A1_Point_Value[i] * A2_Point_Value[i]
    # Now do inverse FFT
    Result = FFT3(Result_Point_Value, True)
    for i in range(len(Result)):
        Result[i] = int(Result[i].real + 0.00001)
    return Result

if __name__ == "__main__":

    # A1 = [2,9,4,2]
    # A2 = [7,9,6,0]
    '''  The following commented-out code gives a way to generate long 
         sequences of digits to experiment with.  Notice how the times for the 
         two methods change as you play with 'length'.  If you have done your
         code correctly, your polyMult should be dramatically faster than
         naiveMult when the length gets large. '''
    length =  20470
    A1 = [random.randint(0,9) for i in range(length)]
    print A1
    A2 = [random.randint(0,9) for i in range(length)]
    print "\n", A2
    t1 = time.clock()
    A3 = polyMult(A1, A2)
    t2 = time.clock()
    A4 = naiveMult(A1, A2)
    t3 = time.clock()
    print "result of naiveMult:  ", A4
    print "result of polyMult:   ", A3
    print "Time for polyMult: ", t2 - t1
    print "Time for naiveMult:", t3 - t2

    # print
    # flag = True
    # for i in range(len(A4)):
    #     if A3[i] != A4[i]:
    #         print i
    #         flag = False
    #         break
    #
    # print "Identical: " + str(flag)
