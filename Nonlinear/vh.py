def myNewtonSolver(gfu,a,f,eps=1e-13,Nmax=30, callback=lambda gfu, it: None,
    printrates=True, printError = True, converge_func=lambda *args : False, 
    dampfactor = 1, inverse="umfpack", pre = None, itToClacJacobi=1):
    """Newton solver for nonlinear problems

    :param gfu: solution vector
    :param a: bilinear form
    :param f: linear form
    :param eps: tolerance, defaults to 1e-13
    :param Nmax: maximum number of iteations, defaults to 10
    :param callback: callback function, defaults to lambda gfu, it: None
    :param printrates: print error in each step, defaults to True
    :param printError: print error if not converging, defaults to True
    :param converge_func: convergance function, defaults to lambda*args:False
    :param dampfactor: damping factor, defaults to 1
    :param inverse: string for the inverse if precon is not set, defaults to "pardiso"
    :param pre: optional preconditioner to solve the system in each iteration, defaults to None
    :return: convergence, it
    """




    from ngsolve import sqrt, InnerProduct, solvers, GridFunction
    res = gfu.vec.CreateVector()
    

    fes = gfu.space
    du = GridFunction(fes)
    gfu_o = GridFunction(fes)
    callback(gfu, -1)
    conv = True
    
    for it in range(Nmax):
        # print ("Iteration {:3}  ".format(it),end="")
        gfu_o.vec.data = gfu.vec

        # if True:
        a.Apply(gfu.vec, res)
        # else:
        #     a.Assemble()
        #     res = a.mat * gfu.vec
        if it%itToClacJacobi == 0:
            a.AssembleLinearization(gfu.vec)

        f.Assemble()

        if type(pre) == type(None):
            du.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse=inverse) * (res - f.vec)
        else:
            solvers.BVP(bf = a, lf = f, gf = du, pre = pre)
        du.vec.data *= dampfactor
        gfu.vec.data -= du.vec
        callback(gfu, it)
        #stopping criteria


        stopcritval = sqrt(abs(InnerProduct(du.vec,res - f.vec)))

        if printrates:
            print ("<A u",it,", A u",it,">_{-1}^0.5 = ", stopcritval)

        if stopcritval < eps or converge_func(gfu, gfu_o):
            break
    else:
        if printError:            
            print(f"warning: too many iterations {it} with {stopcritval:f} >= {eps:f}")
        conv = False

    return conv, it
