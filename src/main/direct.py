import numpy as np
    
class GlobalMin():
    def __init__(self, minimize = True, known = False, val = None):
        self.minimize = minimize
        self.known = known
        self.value = val


class Rectangle():
    def __init__(self, center, f_val, sides):
        self.center  = center
        self.f_val   = f_val
        self.sides   = sides
        self.__str__ = f_val

    @property
    def d2(self):
        return np.sum((self.sides / 2.) ** 2)


class Direct():
    def __init__(self, f, bounds, epsilon=1e-4, max_feval=200, max_iter=10,
                 max_rectdiv=200, globalmin=GlobalMin(), tol = 1e-2):
        self.f           = f        # should take (D,) nparray as input
        self.epsilon     = epsilon  # global/local weight parameter
        self.max_feval   = max_feval
        self.max_iter    = max_iter
        self.max_rectdiv = max_rectdiv
        self.globalmin   = globalmin
        self.tolerance   = tol      # allowable relative error if globalmin is known

        if not self.globalmin.minimize:     # means maximization problem
            self.f_wrap = lambda x: -self.f(x)
        else:
            self.f_wrap = self.f

        self.scale = bounds[:,1] - bounds[:,0]
        self.shift = bounds[:,0]
        self.D = bounds.shape[0]
        
        self.curr_opt      = np.inf
        self.x_at_opt      = None
        self.x_at_opt_unit = None
        self.n_feval       = 0
        self.n_rectdiv     = 0
        self.d_rect        = {}
        self.l_hist        = []

        assert isinstance(bounds, np.ndarray)
        assert len(bounds.shape) == 2
        assert bounds.shape[1]   == 2
        assert np.all(self.scale > 0.)
        #TODO: other assertions

    def u2r(self, unit_coord):
        """unit to real: map a coordinate in unit hyper-cube to one in the actual rectangle"""
        return unit_coord * self.scale + self.shift
    
    def true_sign(self, val):
        return val if self.globalmin.minimize else -val
    
    def divide_rectangle(self, po_rect):
        maxlen      = np.max(po_rect.sides)
        gap         = maxlen / 3.
        d_new_rects = {}

        self.d_rect[po_rect.d2].remove(po_rect)
        maxlen_sides = list(np.nonzero(po_rect.sides == maxlen)[0]) # only the longest sides are divided

        # evaluate points near center
        for side_idx in maxlen_sides:
            d_new_rects[side_idx]   = []
            new_center_u            = po_rect.center.copy()
            new_center_u[side_idx]  += gap
            new_fval_u              = self.f_wrap(self.u2r(new_center_u))

            d_new_rects[side_idx].append(Rectangle(new_center_u, new_fval_u, po_rect.sides.copy()))
            self.l_hist.append((self.u2r(new_center_u), self.true_sign(new_fval_u)))

            self.n_feval   += 1
            self.n_rectdiv += 1
            if new_fval_u < self.curr_opt:
                self.curr_opt      = new_fval_u
                self.x_at_opt_unit = new_center_u.copy()
                self.x_at_opt      = self.u2r(self.x_at_opt_unit)
            if self.globalmin.known and self.error < self.tolerance:
                self.TERMINATE     = True
                break
            if not self.globalmin.known and (self.n_feval >= self.max_feval or self.n_rectdiv >= self.max_rectdiv):
                self.TERMINATE     = True
                break
            
            new_center_l           = po_rect.center.copy()
            new_center_l[side_idx] -= gap
            new_fval_l             = self.f_wrap(self.u2r(new_center_l))
            d_new_rects[side_idx].append(Rectangle(new_center_l, new_fval_l, po_rect.sides.copy()))
            self.l_hist.append((self.u2r(new_center_l), self.true_sign(new_fval_l)))
            self.n_feval   += 1
            self.n_rectdiv += 1
            if new_fval_l < self.curr_opt:
                self.curr_opt      = new_fval_l
                self.x_at_opt_unit = new_center_l.copy()
                self.x_at_opt      = self.u2r(self.x_at_opt_unit)
            if self.globalmin.known and self.error < self.tolerance:
                self.TERMINATE     = True
                break
            if not self.globalmin.known and (self.n_feval >= self.max_feval or self.n_rectdiv >= self.max_rectdiv):
                self.TERMINATE     = True
                break
        if self.TERMINATE:  return

        # axis with better function value get divided first
        maxlen_sides = sorted(maxlen_sides, key=lambda x: min([t.f_val for t in d_new_rects[x]]))
        for side_idx in maxlen_sides:   # po_rect gets divided in every (longest) dimension
            po_rect.sides[side_idx] /= 3.
        for i in maxlen_sides:
            for each_rect in d_new_rects[i]:
                for side_idx in maxlen_sides[i:]:
                    each_rect.sides[i] /= 3.

        for l_rect in d_new_rects.values():
            for each_rect in l_rect:
                d2 = each_rect.d2
                if d2 not in self.d_rect:
                    self.d_rect[d2] = [each_rect]
                else:
                    if each_rect.f_val < self.d_rect[d2][0].f_val:
                        self.d_rect[d2].insert(0, each_rect)
                    else:
                        self.d_rect[d2].append(each_rect)

        # insert po_rect
        if po_rect.d2 not in self.d_rect:
            self.d_rect[po_rect.d2] = [po_rect]
        else:
            if po_rect.f_val < self.d_rect[po_rect.d2][0].f_val:
                self.d_rect[po_rect.d2].insert(0, po_rect)
            else:
                self.d_rect[po_rect.d2].append(po_rect)
        
        # remove empty list from the dictionary
        for dd in [key for key in self.d_rect if len(self.d_rect[key]) == 0]:
            self.d_rect.pop(dd)


    def calc_lbound(self, border):
        lb     = np.zeros(len(border))
        border = np.asarray(border)
        for i in range(len(border)):
            tmp_rects    = [j for j, val in enumerate(border[:,0]) if val < border[i,0]]
            if len(tmp_rects):
                tmp_f    = border[tmp_rects,1]
                tmp_szes = border[tmp_rects,0]
                tmp_lbs  = (border[i,1]-tmp_f)/(border[i,0]-tmp_szes)
                lb[i]    = max(tmp_lbs)
            else:
                lb[i] = -1.976e14
        return lb

    def calc_ubound(self, border):
        ub     = np.zeros(len(border))
        border = np.asarray(border)
        for i in range(len(border)):
            tmp_rects    = [j for j, val in enumerate(border[:,0]) if val > border[i,0]]
            if len(tmp_rects):
                tmp_f    = border[tmp_rects,1]
                tmp_szes = border[tmp_rects,0]
                tmp_ubs  = (tmp_f-border[i,1])/(tmp_szes-border[i,0])
                ub[i]    = min(tmp_ubs)
            else:
                ub[i] = 1.976e14
        return ub

    def get_potentially_optimal_rects(self):
        # among rects with the same size, choose the one with the smallest function value
        border = [(key, l[0].f_val) for key, l in self.d_rect.items()]
        border = sorted(border, key=lambda t:t[0])

        l_po_key = []
        for i in range(len(border)-1):
            if border[i][1] <= border[i+1][1]:
                l_po_key.append(border[i][0])
        l_po_key.append(border[-1][0])

        lbound = self.calc_lbound(border)
        ubound = self.calc_ubound(border)
 
        # find indices of hull that satisfy first condition
        maybe_po = [i for i in range(len(border)) if lbound[i] <= ubound[i]]

        # find indices of hull that satisfy second condition
        po = []
        if self.curr_opt:
            for j in range(len(maybe_po)):
                cond = (self.curr_opt - border[maybe_po[j]][1])/abs(self.curr_opt) + \
                        border[maybe_po[j]][0]*ubound[maybe_po[j]]/abs(self.curr_opt)
                if cond >= self.epsilon:
                    po.append(j)
        else:
            for j in range(len(maybe_po)):
                cond = border[maybe_po[j]][-1] - border[maybe_po[j]][0]*ubound[maybe_po[j]]
                if cond <= 0:
                    po.append(j)

        for i in range(len(po)):
            l_po_key.append(border[maybe_po[po[i]]][0])

        return [self.d_rect[key][0] for key in l_po_key]

            
    def run(self, file):
        D     = self.D
        c     = np.array([0.5] * D)
        f_val = self.f_wrap(self.u2r(c))
        s     = np.array([1.]*D)
        rect  = Rectangle(c, f_val, s)
        
        self.d_rect[rect.d2] = [rect]
        self.l_hist.append((self.u2r(c), self.true_sign(f_val)))
        self.n_feval       += 1
        self.n_rectdiv     += 1
        self.curr_opt      = f_val
        self.error         = self.tolerance
        self.x_at_opt_unit = c
        self.x_at_opt      = self.u2r(c)
        self.TERMINATE     = False
        
        for i in range(self.max_iter):
            # select potentially optimal rectangles
            l_potentially_optimal = self.get_potentially_optimal_rects()
            for po_rect in l_potentially_optimal:
                self.divide_rectangle(po_rect)
                if self.globalmin.value:
                    self.error = (self.curr_opt - self.globalmin.value)/abs(self.globalmin.value)
                else:
                    self.error = self.curr_opt
                if self.TERMINATE:
                    break
            if self.TERMINATE:
                break

        print("number of function evaluations =", self.n_feval)
        file.write("number of function evaluations = "+str(self.n_feval)+"\n")
        opt, x_at_opt, hist = self.true_sign(self.curr_opt), self.x_at_opt, self.l_hist
        print("optimum =", opt, ",  x_at_opt =", x_at_opt, "\n")
        file.write("optimum = "+str(opt)+" ,  x_at_opt = "+str(x_at_opt)+"\n\n")