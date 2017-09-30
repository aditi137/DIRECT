import numpy

def Direct(Problem, bounds, opts, *args):
# Function   : Direct Version 1.0
# Written by : Aditi Dutta (aditi011@e.ntu.edu.sg)
# Created on : 09/12/2017
# Last Update: 09/12/2017
# Purpose    : Direct optimization algorithm.
#
#          x, fmin, history = Direct(Problem, bounds, opts, *args)
#
# Parameters:
#         IN: Problem   - Structure containing problem
#                    Problem.f              = Objective function handle
#
#                    NOTE: If you problem has no constraints (other than
#                          those on the bounds, this is the only field you
#                          need to add to Problem
#
#                    Problem.numconstraints = number of constraints
#                    Problem.constraint(i).func    = i-th constraint handle
#                    Problem.constraint(i).penalty = penalty parameter for
#                                                    i-th constraint
#                    Note: If f returns objective function AND constraints
#                          then set Problem.constraint(1).func = Problem.f
#             bounds    - an n x 2 vector of the lower and upper bounds.
#                         The first column is the lower bounds, and the second
#                         column contains the upper bounds
#             opts      - (optional) MATLAB structure.
#                    opts.ep        = Jones factor                      (default is 1e-4)
#                    opts.maxevals  = max. number of function evals     (default is 20)
#                    opts.maxits    = max. number of iterations         (default is 10)
#                    opts.maxdeep   = max. number of rect. divisions    (default is 100)
#                    opts.testflag  = 1 if globalmin known, 0 otherwise (default is 0)
#                    opts.showits   = 1 if disp. stats shown, 0 oth.
#                                      (default is 1)
#                    opts.globalmin = globalmin (if known)
#                                      (default is 0)
#                    opts.tol       = tolerance for term. if tflag=1
#                                      (default is 0.01)
#                    opts.impcons   = turns on implicit constraint capability
#                                      (default is 0)
#                                     If set to one, objective function
#                                     is expected to return a flag which represents
#                                     the feasibility of the point sampled
#             *args  - (optional) additional arguments to be passed to
#                         objective function
#
#                NOTE: If opts.tflag == 0, maxevals, maxevals and maxdeep are ignored.
#                      DIRECT will stop when the absolute error is less
#                      than tol. Also, preallocation will not occur, and the algorithm
#                      can run slower than if opts.tflag == 1
#                NOTE: opts.maxevals is an approximate stopping condition. DIRECT will
#                      exceed this budget by a slight amount
#
#        OUT: minval    -  minimum value found
#             xatmin    - (optional) location of minimal value
#             history      - (optional) array of iteration historyory, useful for tables and plots
#                The three columns are iteration, fcn evals, and min value found.
#
#                Direct may be called by
#
#                minval = Direct(Problem,bounds)
#                          or
#                with any variation of the optional arguments

#------------------------------------------------------------------
#
# Implementation taken from:
# D.R. Jones, C.D. Perttunen, and B.E. Stuckman. "Lipschitzian
# Optimization Without the Lipschitz Constant". Journal of
# Optimization Theory and Application, 79(1):157-181, October 1993
#
#------------------------------------------------------------------

#-- Initialize the variables --------------------------------------
    lengths      = numpy.array([])
    c            = numpy.array([])
    fc           = numpy.array([])

    con          = numpy.array([])
    szes         = numpy.array([])
    feas_flags   = numpy.array([])
    
    om_lower     = bounds[ : , 0]
    om_upper     = bounds[ : , 1]
    fcncounter   = 0
    perror       = 0
    itctr        = 1
    done         = 0
    g_nargout    = nargout
    n            = bounds.shape[0]
    
    # Determine option values
    if len(args) < 3:
        opts = numpy.array([])
    elif opts.shape[1] == 0:
        opts = numpy.array([])
    getopts(opts, 
            maxits   =  20,        # maximum of iterations
            maxevals =  10,         # maximum # of function evaluations
            maxdeep  =  100,        # maximum number of side divisions
            testflag =  0,          # terminate if within a relative tolerance of f_opt
            globalmin=  0,          # minimum value of function
            ep       =  1e-4,       # global/local weight parameter.
            tol      =  0.01,       # allowable relative error if f_reach is set
            showits  =  1,          # print iteration stats
            impcons  =  0,          # flag for using implicit constraint handling
            pert     =  1e-6)       # perturbation for implicit constraint handling
    theglobalmin = globalmin
    tflag        = testflag
    
    #-- Pre-allocate memory for storage vectors
    if tflag == 0:
        lengths    = numpy.zeros((n, maxevals + floor(.10*maxevals)), float)
        c          = lengths
        fc         = numpy.zeros((1, maxevals + floor(.10*maxevals)), float)
        szes       = fc
        con        = fc
        feas_flags = fc
    
    #-- Call DIRini ---------------------------------------------------
    thirds, lengths, c ,fc, con, feas_flags, minval, xatmin, perror, history, szes, fcncounter, calltype = \
    DIRini(Problem, n, bounds[ : , 0], bounds[ : , 1], lengths, c, fc, con, feas_flags, szes,
           theglobalmin, maxdeep, tflag, g_nargout, impcons, *args)
    
    ret_minval = minval
    ret_xatmin = xatmin
    #-- MAIN LOOP -----------------------------------------------------
    minval = fc(1) + con(1)
    while perror > tol:
        #-- Create list S of potentially optimal hyper-rectangles
        S = find_po(fc[0:fcncounter] + con[0:fcncounter], lengths[:, 0:fcncounter],
                    minval, ep, szes[0:fcncounter])
    
        #-- Loop through the potentially optimal hrectangles -----------
        #-- and divide -------------------------------------------------
        for i in range(S.shape[1]):
            lengths, fc, c, con, feas_flags, szes, fcncounter, success = \
            DIRdivide(bounds[ : , 0], bounds[ : , 1], Problem, S[0,i], thirds, lengths,
                      fc, c, con, feas_flags, fcncounter, szes, impcons, calltype, *args)
       
        #-- update minval, xatmin --------------------------------------
        minval, fminindex =  min(fc[0:fcncounter] + con[0:fcncounter])
        penminval = minval + con(fminindex)
        xatmin = (om_upper - om_lower) * c[:, fminindex-1] + om_lower
        if con(fminindex) > 0 or feas_flags(fminindex) != 0:
            pass #--- new minval is infeasible, don't do anything
        else:
            #--- update return values
            ret_minval = minval
            ret_xatmin = xatmin

       #--see if we are done ------------------------------------------
        if tflag == 1:
          #-- Calculate error if globalmin known
          if theglobalmin != 0:
              perror = 100*(minval - theglobalmin)/abs(theglobalmin)
          else:
              perror = 100*minval
          
        else:
          #-- Have we exceeded the maxits?
            if itctr >= maxits:
                print('Exceeded max iterations. Increase maxits')
                done = 1
          #-- Have we exceeded the maxevals?
            if fcncounter > maxevals:
                print('Exceeded max fcn evals. Increase maxevals')
                done = 1
            if done == 1:
                perror = -1
       
        if max(max(lengths)) >= maxdeep:
          #-- We've exceeded the max depth
            print('Exceeded Max depth. Increase maxdeep')
            perror = -1
       
        if g_nargout == 3:
          #-- Store History
            maxhist = history.shape[0]
            history[maxhist,0] = itctr
            history[maxhist,1] = fcncounter
            history[maxhist,2] = minval
      
        #-- Call replaceinf if impcons flag is set to 1
        if impcons == 1:
          fc = replaceinf(lengths[:, 0:fcncounter], c[:, 0:fcncounter], fc[0:fcncounter],
                          con[0:fcncounter], feas_flags[0:fcncounter], pert)
      
        #-- show iteration stats
        if showits == 1:
            if con(fminindex) > 0 or feas_flags(fminindex) == 1:
                print('Iter: %4i   f_min: %15.10f*    fn evals: %8i\n', itctr, minval, fcncounter)
        else:
            print('Iter: %4i   f_min: %15.10f    fn evals: %8i\n', itctr, minval, fcncounter)
        itctr  = itctr + 1
    
    #-- Return values
    if g_nargout == 2:
        #-- return x*
        final_xatmin = ret_xatmin
    elif g_nargout == 3:
        #-- return x*
        final_xatmin = ret_xatmin
    
        #-- chop off 1st row of history
        history[0:history.shape[0]-1,:] = history[0:history.shape[0],:]
        history = history[0:history.shape[0]-1,:]
    return x, fmin, history

#------------------------------------------------------------------
# Function:   DIRini                                               
# Purpose   : Initialization of Direct                             
#             to eliminate storing floating points                 
#------------------------------------------------------------------
def DIRini(Problem, n, a, b, p_lengths, p_c, p_fc, p_con, p_feas_flags, p_szes,
           theglobalmin, maxdeep, tflag, g_nargout, impcons, *args):
    l_lengths    = p_lengths
    l_c          = p_c
    l_fc         = p_fc
    l_con        = p_con
    l_feas_flags = p_feas_flags
    szes         = p_szes
    
    #-- start by calculating the thirds array
    #-- here we precalculate (1/3)^i which we will use frequently
    l_thirds[0] = 1/3
    for i in range(1, maxdeep):
        l_thirds[i-1] = (1/3)*l_thirds[i-3]
    
    #-- length array will store # of slices in each dimension for
    #-- each rectangle. dimension will be rows; each rectangle
    #-- will be a column
    
    #-- first rectangle is the whole unit hyperrectangle
    l_lengths[:, -1] = numpy.zeros((n,1))
    
    #-- store size of hyperrectangle in vector szes
    szes[0, 0] = 1
    
    #-- first element of c is the center of the unit hyperrectangle
    l_c[:, 0] = numpy.ones((n,1))/2
    
    #-- Determine if there are constraints
    calltype = DetermineFcnType(Problem, impcons)
    
    #-- first element of f is going to be the function evaluated
    #-- at the center of the unit hyper-rectangle.
    #om_point   = abs(b - a).*l_c(:,1)+ a
    #l_fc(1)    = feval(f,om_point,*args{:}
    [l_fc[0], l_con[0], l_feas_flags[0]] = \
        CallObjFcn(Problem, l_c[:, 0], a, b, impcons, calltype, *args)
    fcncounter = 1
    
    #-- initialize minval and xatmin to be center of hyper-rectangle
    xatmin = l_c[:, 0]
    minval   = l_fc[0]
    if tflag == 1:
        if theglobalmin != 0:
            perror = 100*(minval - theglobalmin)/abs(theglobalmin)
        else:
            perror = 100*minval
    else:
        perror = 2

    #-- initialize history
    #if g_nargout == 3
        history[0, 0] = 0
        history[0, 1] = 0
        history[0, 2] = 0
    
    return l_thirds, l_lengths, l_c, l_fc, l_con, l_feas_flags, minval, xatmin, perror,\
            history, szes, fcncounter, calltype

#------------------------------------------------------------------
# Function   :  find_po                                            
# Purpose    :  Return list of PO hyperrectangles                  
#------------------------------------------------------------------
def find_po(fc, lengths, minval, ep, szes):
    
    #-- 1. Find all rects on hub
    diff_szes = sum(lengths, 1)
    tmp_max = max(diff_szes)
    j=1
    sum_lengths = sum(lengths, 1)
    for i in range(1, tmp_max+2):
        tmp_idx = find(sum_lengths==i-1)
        tmp_n, hullidx = min(fc(tmp_idx))
        if hullidx.shape[1] > 0:
            hull(j) = tmp_idx(hullidx)
            j=j+1
            #-- 1.5 Check for ties
            ties = find(abs(fc(tmp_idx)-tmp_n) <= 1e-13)
            if ties.shape[1] > 1:
                mod_ties = find(tmp_idx(ties) != hull(j-1))
                hull = [hull, tmp_idx(ties(mod_ties))]
                j = hull.shape[1] + 1     
    
    #-- 2. Compute lb and ub for rects on hub
    lbound = calc_lbound(lengths, fc, hull,szes)
    ubound = calc_ubound(lengths, fc, hull,szes)
    
    #-- 3. Find indices of hull who satisfy
    #--    1st condition
    maybe_po = find(lbound-ubound <= 0)
    
    #-- 4. Find indices of hull who satisfy
    #--    2nd condition
    t_len  = hull(maybe_po).shape[1]
    if minval != 0:
        po = find((minval-fc(hull(maybe_po)))/abs(minval) + \
                  szes(hull(maybe_po))*ubound(maybe_po)/abs(minval) >= ep)
    else:
        po = find(fc(hull(maybe_po)) - szes(hull(maybe_po))*ubound(maybe_po) <= 0)
    final_pos = hull(maybe_po(po))
    rects = [final_posszes(final_pos)]
    return rects

#------------------------------------------------------------------
# Function   :  calc_ubound                           
# Purpose    :  calculate the ubound used in determining  
#               potentially optimal hrectangles
#------------------------------------------------------------------
def calc_ubound(lengths, fc, hull, szes):   
    hull_length  = hull.shape[1]
    hull_lengths = lengths(:,hull)
    for i in range(1, hull_length+1):
        tmp_rects = find(sum(hull_lengths,1)<sum(lengths(:,hull(i))))
        if tmp_rects.shape[1] > 0
            tmp_f     = fc(hull(tmp_rects))
            tmp_szes  = szes(hull(tmp_rects))
            tmp_ubs   = (tmp_f-fc(hull(i)))./(tmp_szes-szes(hull(i)))
            ub(i)        = min(tmp_ubs)
        else:
            ub(i)=1.976e14   
    return ub

#------------------------------------------------------------------
# Function   :  calc_lbound                                        
# Purpose    :  calculate the lbound used in determining  
#               potentially optimal hrectangles                                
#------------------------------------------------------------------
def calc_lbound(lengths, fc, hull, szes):    
    hull_length  = hull.shape[1]
    hull_lengths = lengths[:, hull-1]
    for i in range(1, hull_length+1):
        tmp_rects = find(sum(hull_lengths,1)>sum(lengths(:,hull(i))))
        if tmp_rects.shape[1] > 0:
            tmp_f     = fc(hull(tmp_rects))
            tmp_szes  = szes(hull(tmp_rects))
            tmp_lbs   = (fc(hull(i))-tmp_f)./(szes(hull(i))-tmp_szes)
            lb(i)     = max(tmp_lbs)
        else:
            lb(i)     = -1.976e14    
    return lb

#------------------------------------------------------------------
# Function   :  DIRdivide                                          
# Purpose    :  Divides rectangle i that is passed in              
#------------------------------------------------------------------
def DIRdivide(a, b, Problem, index, thirds, p_lengths, p_fc, p_c, p_con,
              p_feas_flags, p_fcncounter, p_szes, impcons, calltype, *args):
    lengths    = p_lengths
    fc         = p_fc
    c          = p_c
    szes       = p_szes
    fcncounter = p_fcncounter
    con        = p_con
    feas_flags = p_feas_flags
    
    #-- 1. Determine which sides are the largest
    li     = lengths[:, index-1]
    biggy  = min(li)
    ls     = find(li==biggy)
    lssize = ls.shape[1]
    j = 0
    
    #-- 2. Evaluate function in directions of biggest size
    #--    to determine which direction to make divisions
    oldc       = c(:,index)
    delta      = thirds(biggy+1)
    newc_left  = oldc(:,ones(1,lssize))
    newc_right = oldc(:,ones(1,lssize))
    f_left     = zeros((1,lssize))
    f_right    = zeros((1,lssize))
    for i in range(1, lssize+1):
        lsi               = ls(i)
        newc_left(lsi,i)  = newc_left(lsi,i) - delta
        newc_right(lsi,i) = newc_right(lsi,i) + delta
        [f_left(i), con_left(i), fflag_left(i)]    = CallObjFcn(Problem,newc_left(:,i),a,b,impcons,calltype,*args)
        [f_right(i), con_right(i), fflag_right(i)] = CallObjFcn(Problem,newc_right(:,i),a,b,impcons,calltype,*args)
        fcncounter = fcncounter + 2
    
    w = [min(f_left, f_right)' ls]
    
    #-- 3. Sort w for division order
    [V,order] = sort(w,1)
    
    #-- 4. Make divisions in order specified by order
    for i in range(1, order.shape[0]+1):
       newleftindex  = p_fcncounter+2*(i-1)+1
       newrightindex = p_fcncounter+2*(i-1)+2
       #-- 4.1 create new rectangles identical to the old one
       oldrect = lengths[:, index-1]
       lengths[:, newleftindex-1]   = oldrect
       lengths[:, newrightindex-1]  = oldrect

       #-- old, and new rectangles have been sliced in order(i) direction
       lengths(ls(order(i,1)), newleftindex)  = lengths(ls(order(i,1)), index) + 1
       lengths(ls(order(i,1)), newrightindex) = lengths(ls(order(i,1)), index) + 1
       lengths(ls(order(i,1)), index)         = lengths(ls(order(i,1)), index) + 1
    
       #-- add new columns to c
       c[:, newleftindex-1]  = newc_left[:, order(i)-1]
       c[:, newrightindex-1] = newc_right[:, order(i)-1]
    
       #-- add new values to fc
       fc(newleftindex)  = f_left(order(i))
       fc(newrightindex) = f_right(order(i))
    
       #-- add new values to con
       con(newleftindex)  = con_left(order(i))
       con(newrightindex) = con_right(order(i))
    
       #-- add new flag values to feas_flags
       feas_flags(newleftindex)  = fflag_left(order(i))
       feas_flags(newrightindex) = fflag_right(order(i))
    
       #-- 01/21/04 Dan Hack
       #-- store sizes of each rectangle
       szes(1,newleftindex)  = 1/2*norm((1/3*ones(lengths.shape[0],1)).^(lengths(:,newleftindex)))
       szes(1,newrightindex) = 1/2*norm((1/3*ones(lengths.shape[0],1)).^(lengths(:,newrightindex)))
    
    szes(index) = 1/2*norm((1/3*ones(lengths.shape[0],1)).^(lengths(:,index)))
    pass = 1
    
    return lengths,fc,c,con,feas_flags,szes,fcncounter,pass

#------------------------------------------------------------------
# Function   :  CallConstraints                                    
# Purpose    :  Evaluate Constraints at pointed specified          
#------------------------------------------------------------------
def CallConstraints(Problem, x, a, b, *args):
    #-- Scale variable back to original space
    point = abs(b - a)*x + a
    ret_value = 0
    if isfield(Problem,'constraint'):
        if not isempty(Problem.constraint):
            for i in range(1, Problem.numconstraints+1):
                if (Problem.constraint(i).func).shape[1] == (Problem.f).shape[1]:
                    if double(Problem.constraint(i).func) == double(Problem.f):                   
                        #-- Dont call constraint; value was returned in obj fcn
                        con_value = 0
                    else:
                        con_value = feval(Problem.constraint(i).func, point, *args)
                else:
                    con_value = feval(Problem.constraint(i).func, point, *args)
                if con_value > 0:
                    #-- Infeasible, punish with associated pen. param
                    ret_value = ret_value + con_value*Problem.constraint(i).penalty
    return ret_value

#------------------------------------------------------------------
# Function   :  CallObjFcn                                         
# Purpose    :  Evaluate ObjFcn at pointed specified               
#------------------------------------------------------------------
def CallObjFcn(Problem,x,a,b,impcon,calltype,*args):
    con_value = 0
    feas_flag = 0  
    #-- Scale variable back to original space
    point = abs(b - a).*x+ a
    if calltype == 1:
        #-- No constraints at all
        fcn_value = feval(Problem.f,point,*args{:})
    elif calltype == 2:
        #-- f returns all constraints
        [fcn_value, cons] = feval(Problem.f,point,*args{:})
        for i in range(1, cons.shape[1]+1):
            if cons > 0:
                con_value = con_value + Problem.constraint(i).penalty*cons(i)
    elif calltype == 3:
        #-- f returns no constraint values
        fcn_value = feval(Problem.f, point, *args)
        con_value = CallConstraints(Problem,x,a,b,*args)
    elif calltype == 4:
        #-- f returns feas flag
        [fcn_value,feas_flag] = feval(Problem.f,point,*args)
    elif calltype == 5:
        #-- f returns feas flags, and there are constraints
        [fcn_value,feas_flag] = feval(Problem.f,point,*args)
        con_value = CallConstraints(Problem,x,a,b,*args)
    
    if feas_flag == 1:
        fcn_value = 10^9
        con_value = 0
    
    return fcn_value, con_value, feas_flag        

#------------------------------------------------------------------
# Function   :  replaceinf                                         
# Purpose    :  Assign R. Carter value to given point              
#------------------------------------------------------------------
def replaceinf(lengths,c,fc,con,flags,pert):
    
    #-- Initialize fcn_values to original values
    fcn_values = fc
    
    #-- Find the infeasible points
    infeas_points = find(flags == 1)
    
    #-- Find the feasible points
    feas_points   = find(flags == 0)
    
    #-- Calculate the max. value found so far
    if not isempty(feas_points):
        maxfc = max(fc(feas_points) + con(feas_points))
    else:
        maxfc = max(fc + con) 
    for i in range(1, infeas_points.shape[1]+1):
        if isempty(feas_points):
            #-- no feasible points found yet
            found_points = numpy.array([])
            found_pointsf = numpy.array([])
            index = infeas_points(i)
        else:
            index = infeas_points(i)
    
            #-- Initialize found points to be entire set
            found_points  = c(:,feas_points)
            found_pointsf = fc(feas_points) + con(feas_points)
    
            #-- Loop through each dimension, and find points who are close enough
            for j in range(1, lengths.shape[0]+1):
                neighbors = find(abs(found_points(j, :) - c(j,index)) <= \
                    3^(-lengths(j,index)))
                if not isempty(neighbors):
                    found_points  = found_points(:, neighbors)
                    found_pointsf = found_pointsf(neighbors)
                else:
                    found_points = numpy.array([])
                    found_pointsf = numpy.array([])
                    break
        #-- Assign Carter value to the point
        if not isempty(found_pointsf):
            #-- assign to index the min. value found + a little bit more
            fstar = min(found_pointsf)
            if fstar != 0:
                fcn_values(index) = fstar + pert*abs(fstar)
            else:
                fcn_values(index) = fstar + pert*1
        else:
            fcn_values(index) = maxfc+1
            maxfc             = maxfc+1
    return fcn_values

#------------------------------------------------------------------
# Function   :  DetermineFcnType                                   
# Purpose    :  Determine how constraints are handled              
#------------------------------------------------------------------
def DetermineFcnType(Problem, impcons):    
    retval = 0
    if not isfield(Problem,'constraint') and not impcons:
        #-- No constraints at all
        retval = 1
    if isfield(Problem,'constraint'):
        #-- There are explicit constraints. Next determine where
        #-- they are called
        if not isempty(Problem.constraint):
            if (Problem.constraint(1).func).shape[1] == (Problem.f).shape[1]:
                #-- Constraint values may be returned from objective
                #-- function. Investigate further
                if double(Problem.constraint(1).func) == double(Problem.f):
                    #-- f returns constraint values
                    retval = 2
                else:
                    #-- f does not return constraint values
                    retval = 3
            else:
                #-- f does not return constraint values
                retval = 3
        else:
            if impcons:
                retval = 0
            else:
                retval = 1
    if impcons:
        if not retval:
            #-- only implicit constraints
            retval = 4
        else:
            #-- both types of constraints
            retval = 5
    return retval

#------------------------------------------------------------------
# GETOPTS Returns options values in an options structure
# USAGE
#   [value1,value2,...]=getopts(options,field1,default1,field2,default2,...)
# INPUTS
#   options  : a structure variable
#   field    : a field name
#   default  : a default value
# OUTPUTS
#   value    : value in the options field (if it exists) or the default value
#
# Variables with the field names will be created in the caller's workspace
# and set to the value in the option variables field (if it exists) or to the
# default value.
#
# Example called from a function:
#   getopts(options,'tol',1e-8,'maxits',100)
# where options contains the single field 'tol' with value equal to 1
# The function have two variable defined in the local workspace, tol with a
# value of 1 and maxits with a value of 100.
#
# If options contains a field name not in the list passed to getopts, a
# warning is issued.
#------------------------------------------------------------------
def getopts(options, *args):
    K=fix(len(args)/2)
    if len(args)/2==K:
        error('fields and default values must come in pairs')
    if isa(options,'struct'):
        optstruct=1
    else:
        optstruct=0
    varargout=cell(K,1)
    k=0
    ii=1
    for i in range(1,K):
        if optstruct and isfield(options,*args{ii}):
            assignin('caller',*args{ii},getfield(options,*args{ii}))
            k=k+1
        else:
            assignin('caller',*args{ii},*args{ii+1})
        ii=ii+2
    if optstruct and k != fieldnames(options).shape[0]:
        warning('options variable contains improper fields')
    return varargout