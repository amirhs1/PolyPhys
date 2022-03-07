def sphere_cylinder_intersection(R , r, b):
    """
    computes the volume of intersection of  a sphere and a cylinder. This function is written 
    based on the following article:
    https://www.sciencedirect.com/science/article/pii/0010465590901843
    
    This function can be used to find the local volume fraction of a spherical beads in the radial direction of a cylindrical space where the centers of sphere and the cylinder are not on the same axis.
    
    Here R is smaller than r
    Parameters:
    R: the radius of a sphere
    r: the radius of a cylinder
    b: the smallest distance between the centers of the sphere and the axis of the cylinder
    
    Caution:
    r here is the radius of cylinder but it is the radius of sphere in the article.
    
    Returns:
    V_i: the volume of intersection
    
    Requierment: 
    Numpy and SymPi packages.
    
    """

    xmax = b + R 
    xmin = b - R
    if r > xmax: # For the purpose of the local volume fraction calculation, this case is meaningful.
        #raise Exception("This funtion only works for the systems where the cylinder's radius is smaller than the sphere radius.")
        V_in = 0
    if 0 <= xmin - r: # There  is no intersection if R =< d-r. This works only if d resides on x+ axis.
        #raise Exception('R <= (b-r): There is no intersection.')
        V_in = 0
    else:
        def intersec_cons(R, r, b):
            """
            computes the quantities needed for calculating the volume of intersection of a sphere and a cylinder. This
            function is written based on the following article:
            https://www.sciencedirect.com/science/article/pii/0010465590901843

            Parameters:
            R: the radius of a sphere
            r: the radius of a cylinder
            b: the smallest distance between the centers of the sphere and the axis of the cylinder

            Caution:
            This implementation of the intersection volume is used to computing the volume fraction of several
            spherical monomers in a cylinderical geometry (the bins are cylindrical shells). As a result, r is
            always smaller than R+b: r < R + b

            Returns:
            A: the maximum of R^2 and (r+b)^2
            B: the minimum of R^2 and (r+b)^2
            C: (b-r)^2
            m: the parameter of the elliptic integral (In the main reference, it is shown by k^2 but in Python's 
            scipy package and Wolfram Mathematica the m is used. As a result m=k^2). k is called the elliptic moduleus
            while m is called

            n: the second parameter of the elliptic integral of the third kind. In the main article, it is shown by
            -\alpha^2. Pay attention to the location of the -\alpha^2 in the equation of the elliptic integral of the
            third kind in the article.
            s: (b+r)(b-r)
            """
            A = max(R ** 2, (b + r) ** 2)
            #print ("A:", A)
            B = min(R ** 2, (b + r) ** 2)
            #print ("B:", B)
            C = (b - r) ** 2
            #print ("C:", C)

            if A == C and A != B: # According to my calculation, if A != B and A == C, then k2=+inf.
                k2 = np.inf
                print ("k2:", k2)
            elif A == C and A == B: 
                k2 = np.nan
                #print("k2 cannot be measured.")
            else:
                k2 = (B - C) / (A - C)
                #print ("k2:", k2)

            if C == 0.:
                neg_alpha2 = np.inf
                #print ("neg_alpha2:", neg_alpha2)
            else:
                neg_alpha2 = (B - C) / C # -\alpa^2 in the article
                #print ("neg_alpha2:", neg_alpha2)
            s = (b + r) * (b - r)
            #print ("s:", s)
            return (A, B, C, k2, neg_alpha2, s)

        A, B, C, k2, neg_alpha2, s = intersec_cons(R, r, b)
        V_1 = 4 / 3 * np.pi * (R ** 3) * np.heaviside(r - b, 1.)

        if R > (b + r):
            if C==0.:
                V_2_pi = A ** 2 * s * elliptic_pi(-neg_alpha2, k2).evalf() / C 
                #print('V_2_pi:',V_2_pi)
            else:
                V_2_pi = A ** 2 * s * elliptic_pi(-neg_alpha2, k2).evalf() / C 
            # evalf method converts the result to floating point numerics
            V_2_k = (A * s - (A - B) * (A - C) / 3) * special.ellipk(k2)
            V_2_e = (A - C) * (s + (4 * A - 2 * B - 2 * C) / 3) * special.ellipe(k2)
            V_2 = 4 / (3 * np.sqrt(A - C)) * (V_2_pi - V_2_k - V_2_e)
            V_in = V_1 + V_2
        elif R < (b + r):
            if C==0.:
                V_2_pi = B ** 2 * s * elliptic_pi(-neg_alpha2, k2).evalf() / C
                #print('V_2_pi:',V_2_pi)
            else:
                V_2_pi = B ** 2 * s * elliptic_pi(-neg_alpha2, k2).evalf() / C
            V_2_k = ((A - 2 * B) * s + (A - B) * (3 * B - C - 2 * A) / 3) * special.ellipk(k2)
            V_2_e = (A - C) * (- s + (2 * A + 2 * C - 4 * B) / 3) * special.ellipe(k2)
            V_2 = 4 / (3 * np.sqrt(A - C)) * (V_2_pi + V_2_k + V_2_e)
            V_in = V_1 + V_2
        else :
            # np.arctan or np.arctan2: ?
            # There is warning here if C=0 due to the division by zero in arctan's argument.
            V_2_arc = 4 * R ** 3 / 3 * np.arctan(2 * np.sqrt(b * r) / (b - r)) 
            V_2_sqrt = 4 / 3 * np.sqrt(A - C) * (s + 2 / 3 * (A - C))
            V_2 = V_2_arc - V_2_sqrt
            V_in = V_1 + V_2
            
    return V_in