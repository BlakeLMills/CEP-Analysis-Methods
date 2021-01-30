import math

RADIUS_OF_EARTH_METERS = 6378137.0


def DistanceBetweenTwoLLA(lat1, lon1, alt1, lat2, lon2, alt2):
    #Haversine Formula
    InsideAsin = math.sqrt(math.sin((lat1-lat2)/2.0) * math.sin((lat1-lat2)/2.0) + math.cos(lat1) * math.cos(lat2) * math.sin((lon1-lon2)/2.0) * math.sin((lon1-lon2)/2.0))
    distance = 2.0 * RADIUS_OF_EARTH_METERS * math.asin(InsideAsin)
    alt = alt2 - alt1
    return distance, alt



def SecantMethod(f, x0, x1, iterations):
    # Return the root calculated using the secant method.
    i = 0
    while i < iterations:
        x2 = x1 - f(x1) * (x1 - x0) / float(f(x1) - f(x0))
        x0, x1 = x1, x2
        i += 1
    return x2

def f_example(x):
    return x ** 2 - 612

def CalculateProbability(SigX, SigY, MeanX, MeanY, R, iterations):

    MIN_NUMBER_OF_ITERATIONS = 10
    x_k_sum = 0
    y_j_sum = 0
    Summation = 0
    
    k = 0
    kLoopCheck = 1
    while kLoopCheck > 0.0001:
    #for k in range(0, iterations):
        j_sum = 0
        for j in range(0, iterations):

            x_k_1 = math.factorial(2 * k) / (math.factorial(k + 1) * math.factorial(k) ** 2)
            x_k_2 = ((- R ** 2) / (8 * SigX ** 2)) ** k
            x_k_sum = 0
            for l in range(0, k):
                x_k_l_exp = ((-2 * MeanX ** 2)/ (SigX ** 2)) ** l
                x_k_sum += (math.factorial(k) / (math.factorial(k - l) * math.factorial(2 * l))) * x_k_l_exp
            x_k = x_k_1 * x_k_2 * x_k_sum

            y_j_1 = math.factorial(2 * j) / (math.factorial(j + 1) * math.factorial(j) ** 2)
            y_j_2 = ((- R ** 2) / (8 * SigY ** 2)) ** j
            y_j_sum = 0
            for i in range(0, j):
                y_j_i_exp = ((-2 * MeanY ** 2)/ (SigY ** 2)) ** i
                y_j_sum_addition = (math.factorial(j) / (math.factorial(j - i) * math.factorial(2 * i))) * y_j_i_exp
                y_j_sum += y_j_sum_addition
                if abs(y_j_sum_addition) < 0.0000001:
                    break
            y_j = y_j_1 * y_j_2 * y_j_sum
            if j > 500 and abs(y_j) < 0.0000001:
                print(j, y_j)
                break
            #print(x_k, y_j)
            factorials = (math.factorial(k + 1) * math.factorial(j + 1)) / math.factorial(k + j + 1)
            Addition = x_k * y_j * factorials
            j_sum += Addition
        Summation += j_sum

        k += 1
         
        if k > MIN_NUMBER_OF_ITERATIONS:
            kLoopCheck = abs(j_sum)
            #print("Additional summations would be pointless", Addition, k)
            #break
    
    exponent = (MeanX ** 2) / (2 * SigX ** 2) + (MeanY ** 2) / (2 * SigY ** 2)
    D = 1/(SigX * SigY) * math.exp(-exponent)
    Probablity = ((R ** 2) / 2) * D * Summation

    return Probablity



def ExactMethod(SigX, SigY, MeanX, MeanY, InitialCep, iterations, ProbabilityTarget):

    CEP_i_minus_1 = InitialCep
    Prob_i_minus_1 = CalculateProbability(SigX, SigY, MeanX, MeanY, CEP_i_minus_1, iterations)
    if Prob_i_minus_1 < 0.5:
        CEP_i = 1.05 * CEP_i_minus_1
    else:
        CEP_i = 0.95 * CEP_i_minus_1
    Prob_i = CalculateProbability(SigX, SigY, MeanX, MeanY, CEP_i, iterations)
    while abs(Prob_i - ProbabilityTarget) > 0.00000001:

        try: 
            CEP_i_plus_1 = CEP_i + ((CEP_i_minus_1 - CEP_i) * (ProbabilityTarget - Prob_i)/(Prob_i_minus_1 - Prob_i))
            #print(CEP_i_minus_1, CEP_i, CEP_i_plus_1, Prob_i, Prob_i_minus_1)
        except:
            CEP_i_plus_1 = CEP_i
            print("Unable to calculate CEP i+1")
            break
        CEP_i_minus_1 = CEP_i
        Prob_i_minus_1 = Prob_i
        CEP_i = CEP_i_plus_1
        
        Prob_i = CalculateProbability(SigX, SigY, MeanX, MeanY, CEP_i_plus_1, iterations)
        print(Prob_i)
    return CEP_i
        
def GrubbsPatnaikWilsonHilferty(SigX, SigY, MeanX, MeanY, Corr):

    m = SigX ** 2 + SigY ** 2 + MeanX ** 2 + MeanY ** 2
    v = 2 * (SigX ** 4 + 2 * (Corr ** 2) * (SigX ** 2) * (SigY ** 2) + SigY ** 4) + 4 * ((MeanX ** 2) * (SigX ** 2) + 2 * MeanX * MeanY * SigX * SigY + (MeanY ** 2) * (SigY ** 2))

    CEP = math.sqrt(m * (1 - (v/(9 * (m ** 2))) ** 3))
    return CEP


if __name__ == "__main__":
    root = SecantMethod(f_example, 10, 30, 5)
    print("Root: {}".format(root))

    x_bar = 120.4
    Sig_x = 165.5
    y_bar = 134.4
    Sig_y = 255.9
    Corr = 0.78

    x_bar_prime = 37.1
    y_bar_prime = 281.8
    Sig_x_prime = 91.1
    Sig_y_prime = 290.7
    Corr_prime = 0.0
    

    #CEP_GPWH = GrubbsPatnaikWilsonHilferty(Sig_x_prime, Sig_y_prime, x_bar_prime, y_bar_prime, 0)
    CEP_GPWH = GrubbsPatnaikWilsonHilferty(Sig_x, Sig_y, x_bar, y_bar, Corr)
    CEP_Exact = ExactMethod(Sig_x_prime, Sig_y_prime, x_bar_prime, y_bar_prime, CEP_GPWH, 300, 0.5)
    #CEP_Exact = ExactMethod(Sig_x, Sig_y, x_bar, y_bar, 280.34, 2000, 0.5)
    print("CEP Exact:", CEP_Exact)
    print("CEP GPWH:", CEP_GPWH)
    print(.614 * Sig_x + .563 * Sig_y)


