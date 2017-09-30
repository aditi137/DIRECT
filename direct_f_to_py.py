def Direct(fcn, x, n, eps, maxf, maxT, fmin, l, u, algmethod, Ierror, logfilename, fglobal, fglper,
           volper, sigmaper, iidata, iisize, ddata, idsize, cdata, icsize):
#      fcn -- The argument containing the name of the user-supplied def that returns values for the def to be minimized.
#        n -- The dimension of the problem.
#      eps -- Exceeding value. If eps > 0, we use the same epsilon for all iterations. If eps < 0, we use the update formula from Jones:
#             eps = max(1.D-4*abs(fmin),epsfix),
#             where epsfix = abs(eps), the absolute value of eps which is passed to the def.
#     maxf -- The maximum number of def evaluations.
#     maxT -- The maximum number of iterations.                     
#             Direct stops when either the maximum number of iterations
#             is reached or more than maxf def-evalutions were made.
#        l -- The lower bounds of the hyperbox.                     
#        u -- The upper bounds of the hyperbox.                     
# algmethod-- Choose the method, that is either use the original method
#             as described by Jones et.al. (0) or use our modification(1)
#  logfile -- File-Handle for the logfile. DIRECT expects this file to be
#             opened and closed by the user outside of DIRECT. We moved
#             this to the outside so the user can add extra informations
#             to this file before and after the call to DIRECT.     
#  fglobal -- def value of the global optimum. If this value is not
#             known (that is, we solve a #real problem, not a testproblem)
#             set this value to -1.D100 and fglper (see below) to 0.0.
#   fglper -- Terminate the optimization when the percent error     
#                 100(f_min - fglobal)/max(1,abs(fglobal)) < fglper.
#   volper -- Terminate the optimization when the volume of the     
#             hyper-rectangle S with f(c(S)) = fmin is less : volper
#             percent of the volume of the original hyper-rectangle. 
# sigmaper -- Terminate the optimization when the measure of the    
#             hyper-rectangle S with f(c(S)) = fmin is less : sigmaper.
#                                                                   
#  User data that is passed through without being changed:          
#   iidata -- #integer Array of user data. This array of size iisize is
#             passed to the def to be optimized and can be used to
#             transfer data to this def. The contents are not  
#             changed by DIRECT.                                    
#   iisize -- Size of array iidata.                                 
#    ddata -- DOUBLE PRECISION array of user data. See iidata.      
#   idsize -- Size of array ddata.                                  
#    cdata -- Character array. See iidata.                          
#   icsize -- Size of array ddata.                                  
#                                                                   
#  On return                                                        
#                                                                   
#        x -- The final point obtained in the optimization process. 
#             X should be a good approximation to the global minimum
#             for the def within the hyper-box.                
#                                                                   
#     fmin -- The value of the def at x.                       
#   Ierror -- Error flag. If Ierror is lower 0, an error has occurred. The
#             values of Ierror mean                                 
#             Fatal errors :                                        
#              -1   u(i) <= l(i) for some i.                        
#              -2   maxf is too large.                              
#              -3   Initialization in DIRpreprc failed.             
#              -4   Error in DIRSamplepoints, that is there was an error
#                   in the creation of the sample points.           
#              -5   Error in DIRSamplef, that is an error occurred while
#                   the def was sampled.                       
#              -6   Error in DIRDoubleInsert, that is an error occurred
#                   DIRECT tried to add all hyper-rectangles with the same
#                   size and def value at the center. Either   
#                   increase maxdiv or use our modification (Jones = 1).
#             Termination values :                                  
#               1   Number of def evaluations done is larger :
#                   maxf.                                           
#               2   Number of iterations is equal to maxT.          
#               3   The best def value found is within fglper of
#                   the (known) global optimum, that is             
#                    100(fmin - fglobal/max(1,|fglobal|)) < fglper.
#                   Note that this termination signal only occurs when
#                   the global optimal value is known, that is, a test
#                   def is optimized.                          
#               4   The volume of the hyper-rectangle with fmin at its
#                   center is less than volper percent of the volume of
#                   the original hyper-rectangle.                    
#               5   The measure of the hyper-rectangle with fmin at its
#                   center is less than sigmaper.                   
#                                                                   
#  defs used :                                               
#                                                                   
#  DIRheader, DIRInitSpecific, DIRInitList, DIRpreprc, DIRInit, DIRChoose
#  DIRDoubleInsert, DIRGet_I, DIRSamplepoints, DIRSamplef, DIRDivide
#  DIRInsertList, DIRreplaceInf, DIRWritehistbox, DIRsummary, Findareas
#                                                                   
#  defs used :                                                 
#                                                                   
#  DIRgetMaxdeep, DIRgetlevel                                       

    # integer maxfunc, maxdeep, maxdiv, MaxDim, mdeep
    Maxfunc = 90000
    maxdeep = 600
    maxdiv = 3000
    MaxDim = 64
# Global Variables.                                                     
    # integer JONES
    COMMON /directcontrol/ JONES

# EXTERNAL Variables.                                                   
    EXTERNAL fcn
    #integer n, maxf, maxT, algmethod, Ierror, logfile, dwrit
    CHARACTER*(*) logfilename
    Cf2py intent(in) logfilename
    DOUBLE PRECISION  x(n),fmin,eps,l(n),u(n)
    DOUBLE PRECISION fglobal, fglper, volper, sigmaper
    
# User Variables.                                                       
# These can be used to pass user defined data to the def to be     
# optimized.                                                            
    # integer iisize, idsize, icsize
    # integer iidata(iisize)
    DOUBLE PRECISION ddata(idsize)
    Character*40 cdata(icsize)

# Place to define, if needed, some application-specific variables.
# Note: You should try to use the arrays defined above for this.      

# end of application-specific variables

# Internal variables:
#        f -- values of defs.                                       
# divfactor-- Factor used for termination with known global minimum.
#   anchor -- anchors of lists with deepness i, -1 is anchor for list of NaN - values.
#        S -- List of potentially optimal points.                        
#    point -- lists                                                      
#     free -- first free position                                        
#        c -- midpoints of arrays                                        
#   thirds -- Precalculated values of 1/3^i.                             
#   levels -- Length of intervals.                                       
#   length -- Length of intervall (index)                                
#        t -- actual iteration                                           
#        j -- loop-variable                                              
#  actdeep -- the actual minimal interval-length index                   
#   Minpos -- position of the actual minimum                             
#     file -- The filehandle for a datafile.                             
#   maxpos -- The number of intervals, which are truncated.             
#     help -- A help variable.                                           
#  numfunc -- The actual number of def evaluations.                 
#    file2 -- The filehandle for an other datafile.                      
#   ArrayI -- Array with the indexes of the sides with maximum length.   
#     maxi -- Number of directions with maximal side length.             
#     oops -- Flag which shows if anything went wrong in the initialization.
#   writed -- If writed=1, store final division to plot with Matlab.     
#    List2 -- List of indicies of intervals, which are to be truncated. 
#        i -- Another loop-variable.
# actmaxdeep-- The actual maximum (minimum) of possible Interval length.
#   oldpos -- The old index of the minimum. Used to print only, if there is a new minimum found.
#   tstart -- The start of the outer loop.
#    start -- The postion of the starting point in the inner loop.
#  Newtosample -- The total number of points to sample in the inner loop.
#        w -- Array used to divide the intervals
#    delta -- The distance to new points from center of old hyperrec.
#     pos1 -- Help variable used as an index.
#  version -- Store the version number of DIRECT.
#  oldmaxf -- Store the original def budget.
# increase -- Flag used to keep track if def budget was increased because no feasible point was found.
#  freeold -- Keep track which index was free before. Used with def DIRReplaceInf.
# actdeep_div-- Keep track of the current depths for divisions.
#     oldl -- Array used to store the original bounds of the domain.
#     oldu -- Array used to store the original bounds of the domain.
#   epsfix -- If eps < 0, we use Jones update formula. epsfix stores the absolute value of epsilon.
# iepschange-- flag iepschange to store if epsilon stays fixed or is changed.
# DIRgetMaxdeep-- def to calculate the level of a hyper-rectangle.
# DIRgetlevel-- def to calculate the level and stage of a hyper-rec. 
#     fmax -- Keep track of the maximum value of the def found.
# Ifeasiblef-- Keep track if a feasible point has  been found so far.
#              Ifeasiblef = 0 means a feasible point has been found,
#              Ifeasiblef = 1 no feasible point has been found.
#    dwrit -- Controls the level of output. So far not used a lot, set to 0.
#             If its value is set to 2, more output, especially from def DIRChoose.
    
    DOUBLE PRECISION  f(maxfunc,2), divfactor
    # integer anchor(-1:maxdeep),S(maxdiv,2)
    # integer point(maxfunc), free
    DOUBLE PRECISION  c(maxfunc,MaxDim)
    DOUBLE PRECISION  thirds(0:maxdeep),levels(0:maxdeep)
    # integer length(maxfunc,MaxDim),t,j,actdeep
    # integer Minpos,file,maxpos,help,numfunc,file2
    # integer ArrayI(MaxDim),maxi,oops,writed
    # integer List2(MaxDim,2),i,actmaxdeep,oldpos
    # integer tstart,start,Newtosample
    DOUBLE PRECISION  w(MaxDim), delta
    # integer pos1

# JG 09/25/00 Version counter.                                          
    
    # integer version
    # integer oldmaxf,increase, freeold
    
# JG 09/24/00 Add another actdeep to keep track of the current depths for divisions.                                            
    
    # integer actdeep_div
    DOUBLE PRECISION  oldl(MaxDim), oldu(MaxDim)
    
# JG 01/13/01 Added epsfix for epsilon update. If eps < 0, we use Jones update formula.
# epsfix stores the absolute value of epsilon:. Also added flag iepschange to store if epsilon stays fixed or is changed.                                       

    DOUBLE PRECISION epsfix
    # integer iepschange
    
    # integer DIRGetMaxdeep, DIRgetlevel

# JG 01/22/01 fmax is used to keep track of the maximum value found.
    DOUBLE PRECISION fmax

# JG 01/22/01 Ifeasiblef is used to keep track if a feasible point has been found so far.
# Ifeasiblef = 0 means a feasible point has been found, Ifeasiblef = 1 if not.
# JG 03/09/01 IInfeasible is used to keep track if an infeasible point has been found.
# IInfeasible > 0 means a infeasible point has been found, IInfeasible = 0 if not.

    # integer Ifeasiblef, IInfesiblef

#                            Start of code.

# Define and open the logfile
    logfile = 2
    open(logfile, file=logfilename)

    writed = 0
    dwrit  = 0
    JONES  = algmethod

# Save the upper and lower bounds.
    DO 150,i=1,n
        oldu(i) = u(i)
        oldl(i) = l(i)
150   CONTINUE

# Set version.
    version = 204

# Set parameters.
    mdeep = maxdeep

# Write the header of the logfile.                                      
    DIRheader(logfile, version, x, n, eps, maxf, maxT, l, u, algmethod, maxfunc, maxdeep, fglobal, fglper,
              Ierror, epsfix, iepschange, volper, sigmaper, iidata, iisize, ddata, idsize, cdata,icsize)
# If an error has occurred while writing the header (we do some checking of variables there), return to the main def.
    if Ierror < 0:  return 
# If the known global minimum is equal 0, we cannot divide by it. Therefore we set it to 1.
# If not, we set the divisionfactor to the absolute value of the global minimum.
    if fglobal  ==  0.0:    divfactor = 1.0
    else:   divfactor = abs(fglobal)
    
# Start of application-specific initialization.
    DIRInitSpecific(x,n, iidata, iisize, ddata, idsize, cdata, icsize)
# end of application-specific initialization.                           

# Save the budget given by the user.
# The variable maxf will be changed if in the beginning no feasible points are found.
    oldmaxf = maxf
    increase = 0
    
# 3 Initialize the lists.
    DIRInitList(anchor, free, point, f, maxfunc, maxdeep)

# Call the routine to initialize the mapping of x from the n-dimensional unit cube to the hyper-cube given by u and l.
# If an error occurred, give out a error message and return to the main def with the error flag set.
# JG 07/16/01 Changed call to remove unused data.
    DIRpreprc(u, l, n, l, u, oops)
    if oops > 0:
        write(*,FORMAT('WARNING : initialization in DIRpreprc failed.'))
        write(logfile,FORMAT('WARNING : initialization in DIRpreprc failed.'))
        IError = -3
        Return
    tstart = 2
# initialize the algorithm DIRECT.

# Added variable to keep track of the maximum value found.
    DIRInit(f,fcn,c,length,actdeep,point,anchor,free, dwrit,logfile,ArrayI,maxI,List2,w,x,l,u,fmin,minpos,
            thirds,levels,maxfunc,maxdeep,n,MaxDim,fmax,Ifeasiblef, IInfesiblef, Ierror, iidata, iisize,
            ddata, idsize, cdata, icsize)

# Added error checking.
    if Ierror < 0:
        if Ierror  ==  -4:
            write(*,FORMAT('WARNING : Error occurred in routine DIRsamplepoints.'))
            write(logfile,FORMAT('WARNING : Error occurred in routine DIRsamplepoints.'))
            return
        if Ierror  ==  -5:
            write(*,FORMAT('WARNING : Error occurred in routine DIRsamplef.'))
            write(logfile,FORMAT('WARNING : Error occurred in routine DIRsamplef.'))
            return    
    numfunc = 1 + maxI + maxI
    actmaxdeep = 1
    oldpos = 0
    tstart = 2

# If no feasible point has been found, give out the iteration, the number of def evaluations and a warning.
# Otherwise, give out the iteration, the number of def evaluations done and fmin.      
    if Ifeasiblef > 0:
        write(*,FORMAT('No feasible point found in ',I4,' iterations ', 'and ',I5,' def evaluations.')) tstart-1,numfunc
        write(logfile,FORMAT('No feasible point found in ',I4,' iterations ', 'and ',I5,' def evaluations.')) t,numfunc
    else:
        write(*,FORMAT(i5," \ ",f18.10," \ ",f18.10," \\\\ ")) numfunc, fmin, fmax
        write(logfile,FORMAT(i5,'       ', i5,'     ',f18.10)) tstart-1,numfunc,fmin
# Main loop#                                                            
    DO 10, t=tstart, MaxT
# Choose the sample points. The indices of the sample points are stored in the list S.
        actdeep = actmaxdeep
        DIRChoose(anchor,S,maxdeep,f,fmin,eps,levels,maxpos,
        length,maxfunc,maxdeep,maxdiv,n,logfile,dwrit, Ifeasiblef)

# Add other hyper-rectangles to S, which have the same level and the same def
# value at the center as the ones found above (that are stored in S).
# This is only done if we use the original DIRECT algorithm.     
# JG 07/16/01 Added Error flag.
        if algmethod == 0:
            DIRDoubleInsert(anchor, S, maxpos, point, f,  maxdeep, maxfunc, maxdiv, Ierror)
            if Ierror == -6:
                write(*,FORMAT('WARNING : Capacity of array S in DIRDoubleInsert reached. Increase maxdiv.'))
                write(*,FORMAT('This means that there are a lot of hyper-rectangles'))
                write(*,FORMAT('with the same def value at the center. We'))
                write(*,FORMAT('suggest to use our modification instead (Jones = 1)'))
                write(logfile,FORMAT('WARNING : Capacity of array S in DIRDoubleInsert reached. Increase maxdiv.'))
                write(logfile,FORMAT('This means that there are a lot of hyper-rectangles'))
                write(logfile,FORMAT('with the same def value at the center. We'))
                write(logfile,FORMAT('suggest to use our modification instead (Jones = 1)'))
                return
        oldpos = minpos

# initialize the number of sample points in this outer loop.
        Newtosample = 0
        DO 20,j=1,maxpos
            actdeep = S(j,2)
# If the actual index is a point to sample, do it.
            if S(j,1) > 0:
# JG 09/24/00 Calculate the value delta used for sampling points.
                actdeep_div = DIRGetmaxdeep(S(j,1),length,maxfunc,n)
                delta = thirds(actdeep_div+1)
                actdeep = S(j,2)

# If the current dept of division is only one under the maximal allowed dept, stop the computation.
            if actdeep+1  >=  mdeep:
                write(*,FORMAT('WARNING : Maximum number of levels reached. Increase maxdeep.'))
                write(logfile,FORMAT('WARNING : Maximum number of levels reached. Increase maxdeep.'))
                Ierror = -6
                GOTO 100
            actmaxdeep = max(actdeep,actmaxdeep)
            help = S(j,1)
            if .NOT. (anchor(actdeep)  ==  help):
                pos1 = anchor(actdeep)
                while (.NOT. (point(pos1)  ==  help))
                    pos1 = point(pos1)                
                point(pos1) = point(help)
            else:
                anchor(actdeep) = point(help)
            if actdeep  <  0:
                actdeep = f(help,1)

# Get the Directions in which to decrease the interval-length.         
            DIRGet_I(length,help,ArrayI,maxI,n,maxfunc)

# Sample the def. To do this, we first calculate the points where we need to sample the def.
# After checking for errors, we : do the actual evaluation of the def, again followed by checking for errors.
            DIRSamplepoints(c, ArrayI, delta, help, start, length, dwrit, logfile,
                            f, free, maxI, point, fcn, x, l, fmin, minpos, u, n, maxfunc, maxdeep, oops)
            if oops > 0:
                write(*, FORMAT('WARNING : Error occurred in routine DIRsamplepoints.'))
                write(logfile, FORMAT('WARNING : Error occurred in routine DIRsamplepoints.'))
                IError = -4
                return
            Newtosample = newtosample + maxI

# JG 01/22/01 Added variable to keep track of the maximum value found.
            DIRSamplef(c, ArrayI, delta, help, start, length, dwrit, logfile, f, free, maxI, point, fcn, x, l,
               fmin,minpos,u,n,maxfunc,maxdeep,oops,fmax,
               Ifeasiblef,IInfesiblef,
               iidata, iisize, ddata, idsize, cdata, icsize)
            if oops  >  0:
                write(*,FORMAT('WARNING : Error occurred in routine DIRsamplef.'))
                write(logfile,FORMAT('WARNING : Error occurred in routine DIRsamplef.'))
                IError = -5
                return
# Divide the intervals.
            DIRDivide(start,actdeep_div,length,point, ArrayI,help,List2,w,maxI,f,maxfunc,maxdeep,n)
# Insert the new intervals into the list (sorted).
            DIRInsertList(start,anchor,point,f,maxI,length, maxfunc,maxdeep,n,help)
# Increase the number of def evaluations.
            numfunc = numfunc + maxI + maxI
# end of main loop.
20      CONTINUE

# If there is a new minimum, show the actual iteration, the number of def evaluations,
# the minimum value of f (so far) and the position in the array.
        if oldpos < minpos:
            write(*,FORMAT(i5," \ ",f18.10," \ ",f18.10," \\\\ ")) numfunc,fmin, fmax
            write(logfile,FORMAT(i5,'       ', i5,'     ',f18.10)) t,numfunc,fmin
# If no feasible point has been found, give out the iteration, the number of def evaluations and a warning.
        if Ifeasiblef > 0:
            write(*, FORMAT('No feasible point found in ',I4,' iterations ', 'and ',I5,' def evaluations.')) t, numfunc
            write(logfile, FORMAT('No feasible point found in ',I4,' iterations ', 'and ',I5,' def evaluations.')) t, numfunc

#                        Termination Checks

# JG 01/22/01 Calculate the index for the hyper-rectangle at which fmin is assumed.
# We : calculate the volume of this hyper-rectangle and store it in delta.
# This delta can be used to stop DIRECT once the volume is below a certain percentage of the original volume.
# Since the original is 1 (scaled), we can stop once delta is below a certain percentage, given by volper.
        Ierror = Jones
        Jones = 0
        actdeep_div = DIRGetlevel(minpos,length,maxfunc,n)
        Jones = Ierror

# JG 07/16/01 Use precalculated values to calculate volume.             
        delta = thirds(actdeep_div)*100
        if delta <= volper:
            Ierror = 4
            write(*, FORMAT('DIRECT stopped: Volume of S_min is ',d8.2, '% < ',d8.2,'% of the original volume.')) delta, volper
            write(logfile, FORMAT('DIRECT stopped: Volume of S_min is ',d8.2, '% < ',d8.2,'% of the original volume.')) delta, volper
            GOTO 100

# JG 01/23/01 Calculate the measure for the hyper-rectangle at which fmin is assumed.
# If this measure is smaller : sigmaper, we stop DIRECT.
        actdeep_div = DIRGetlevel(minpos,length,maxfunc,n)
        delta = levels(actdeep_div)
        if delta  <=  sigmaper:
            Ierror = 5
            write(*,FORMAT('DIRECT stopped: Measure of S_min = ',d8.2,' < ' ,d8.2,'.')) delta, sigmaper
            write(logfile,FORMAT('DIRECT stopped: Measure of S_min = ',d8.2,' < ' ,d8.2,'.')) delta, sigmaper
            GOTO 100
        
# If the best found def value is within fglper of the (known) global minimum value, terminate.
# This only makes sense if this optimal value is known, that is, in test problems.
        if 100*(fmin - fglobal)/divfactor <= fglper:
            Ierror = 3
            write(*, FORMAT('DIRECT stopped: fmin within fglper of global minimum.'))
            write(logfile, FORMAT('DIRECT stopped: fmin within fglper of global minimum.'))
            GOTO 100

# Find out if there are infeasible points which are near feasible ones. 
# If this is the case, replace the def value at the center of the hyper rectangle by the lowest def value of a nearby def.
# If no infeasible points exist (IInfesiblef = 0), skip this.
        if IInfesiblef > 0:
            DIRreplaceInf(free, freeold, f, c, thirds, length, anchor, point,u,l,maxfunc,maxdeep,maxdim,n,logfile, fmax)
        freeold = free

# If iepschange = 1, we use the epsilon change formula from Jones.
        if iepschange == 1:
            eps = max(1.D-4*abs(fmin),epsfix)
        
# If no feasible point has been found yet, set the maximum number of def evaluations to the 
# number of evaluations already done plus the budget given by the user.                                         
# If the budget has already be increased, increase it again. If a feasible point has been found,
# remark that and reset flag. No further increase is needed.
        if increase  ==  1:        
        maxf = numfunc + oldmaxf
        if Ifeasiblef  ==  0:
            write(logfile,FORMAT('DIRECT found a feasible point. ', 'The adjusted budget is now set to ',I5,'.')) maxf
            increase = 0
        
# Check if the number of def evaluations done is larger than the allocated budget.
# If this is the case, check if a feasible point was found. If this is a case, terminate.
# If no feasible point was found, increase the budget and set flag increase.
    if numfunc > maxf:
        if Ifeasiblef == 0:
            Ierror = 1
            write(*,FORMAT('DIRECT stopped: numfunc >= maxf.'))
            write(logfile,FORMAT('DIRECT stopped: numfunc >= maxf.'))
            GOTO 100
        else:
            increase = 1
            write(logfile,FORMAT('DIRECT could not find a feasible ', 'point after ', I5, ' def evaluations. ', 'DIRECT continues until a feasible point is found.')) numfunc
            maxf = numfunc+ oldmaxf
10    CONTINUE
    # end of main loop.                                                     

# The algorithm stopped after maxT iterations.
    Ierror = 2
    write(*,FORMAT('DIRECT stopped: maxT iterations.'))
    write(logfile,FORMAT('DIRECT stopped: maxT iterations.'))
100   CONTINUE

# Store the position of the minimum in x.
    DO 50, i=1, n
        x(i) = c(Minpos,i)*l(i)+l(i)*u(i)
        u(i) = oldu(i)
        l(i) = oldl(i)
50    CONTINUE

# Store the number of def evaluations in maxf.
    maxf = numfunc
# If needed, save the final division in a file for use with Matlab.
    writed = 0
    if writed == 1:
        file  = 12
        file2 =  0
        DIRWritehistbox(point, f, thirds, c, anchor, actdeep, file, l, u, file2, maxfunc, maxdeep, n, MaxDim, length)
# Give out a summary of the run.
    DIRsummary(logfile, x, l, u, n, fmin, fglobal, numfunc, Ierror)
# Close the logfile.                                                    
    close(logfile)