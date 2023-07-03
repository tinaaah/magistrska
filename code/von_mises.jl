function rand_von_Mises(N, mu, kappa)

    """
        rand_von_Mises(N,mu,kappa)
        ==========================

        Generates theta an Nx1 array of samples of a von Mises distribution
        with mean direction mu and concentration kappa.

        INPUT:

            * N - number of samples to be generated (integer)
            * mu - mean direction (float)
            * kappa - concentration (float)

        OUTPUT:

            * theta - an Nx1 array of samples of a von Mises distribution with mean
            direction mu and concentration kappa.

         References:
         ===========

         Algorithm first given in

         [1] D. J. Best and N. I. Fisher, Efficient Simulation of the von Mises
         Distribution, Applied Statistics, 28, 2, 152--157, (1979).

         Also given in the following textbook/monograph

         [2] N. I. Fisher, Statistical analysis of circular data, Cambridge University Press, (1993).

    """

    # MAIN BODY OF ALGORITHM
    # =======================

    # Used same notation as Ref.~[2], p49

    a = 1.0 + sqrt(1.0 + 4.0 * kappa^2)
    b = (a - sqrt(2.0 * a)) / (2.0 * kappa)
    r = (1.0 + b^2) / (2.0 * b)

    counter = 1
    theta = zeros((N,1))

    while counter <= N

        # Pseudo-random numbers sampled from a uniform distribution [0,1]
        U1 = rand()
        U2 = rand()
        U3 = rand()

        z = cos(pi * U1)
        f = (1.0 + r *z) / (r + z)
        c = kappa * (r - f)

        if ((c * (2.0 - c) - U2) > 0.0) | ((log(c/U2) + 1.0 - c) > 0.0)

            theta[counter] = mod(sign(U3 - 0.5) * acos(f) + mu, 2*pi)     
            counter += 1
        end
    end

    return theta
end