(in-package :cl-user)

(defpackage :numerical.interpolation
  (:use :common-lisp :numerical)
  (:import-from :numerical
                :/_
                :+tolerance+
                :+limit+)
  (:export :lagrange-polynomial
           :lagrange/order
           :hermite/cubic
           :cubic-spline
           :cubic-spline-prime
           :first-derivative/2
           :first-derivative/4
           :second-derivative/3
           :second-derivative/5
           :first-derivative
           :first-derivative+
           :second-derivative
           :second-derivative+))

(in-package :numerical.interpolation)

(defconstant +step+ .04d0) ; initial step size to use for derivative determination

(defun bracket/3 (xx lx lfx ldfx &optional (order 2))
  (let* ((lx2 (mapcar (lambda (n) (abs (- xx n))) lx))
         (z (loop for x in lx
               for x2 in lx2
               for dfx in ldfx
               for fx in lfx
               as i = 0 then (1+ i) collect (list x2 fx x dfx i)))
         (sz (sort (subseq (sort z (lambda (x y) (< (car x) (car y))))
                           0 order)
                   (lambda (x y)
                     (< (third x) (third y))))))
    (list (loop for x in sz collect (fifth x))
          (loop for x in sz collect (third x))
          (loop for x in sz collect (second x))
          (loop for x in sz collect (fourth x)))))

(defun lagrange-polynomial (lx lfx)
  "Create a lambda to approximate f(xx) given a list of values and their function evaluations."
  (flet ((ljn (xx xj j)
           (let ((v (loop for x in lx as i = 0 then (1+ i)
                       when (not (equal i j)) collect x)))
             (/_ (apply #'* (mapcar (lambda (x) (- xx x)) v))
                 (apply #'* (mapcar (lambda (x) (- xj x)) v))))))
    (lambda (xx)
      (loop
         for x in lx
         for fx in lfx
         as j = 0 then (1+ j)
         summing (* (ljn xx x j) fx)))))

(defun lagrange/order (lx lfx order)
  "Create a lambda to approximate f(xx) given a list of values and their function evaluations.
Only the evaluations associated with the nearest `ORDER` values will be used to construct
the approximating polynomial."
  (flet ((get-order (xx)
           (let* ((lx2 (mapcar (lambda (n) (abs (- xx n))) lx))
                  (z (loop for x in lx
                        for x2 in lx2
                        for fx in lfx collect (list x2 fx x)))
                  (sz (sort (subseq (sort z (lambda (x y) (< (car x) (car y))))
                                    0 order)
                            (lambda (x y)
                              (< (third x) (third y))))))
             (list (loop for x in sz collect (third x))
                   (loop for x in sz collect (second x))))))
    (lambda (xx)
      (funcall (apply #'lagrange-polynomial (get-order xx)) xx))))

(defun hermite/cubic (lx lfx ldfx)
  "Create a lambda to approximate f(xx) using a Hermite cubic
polynomial given a list of values, their function evaluations, and
their function derivate evaluations."
  (flet ((get-order (xx order)
           (let* ((lx2 (mapcar (lambda (n) (abs (- xx n))) lx))
                  (z (loop for x in lx
                        for x2 in lx2
                        for dfx in ldfx
                        for fx in lfx collect (list x2 fx x dfx)))
                  (sz (sort (subseq (sort z (lambda (x y) (< (car x) (car y))))
                                    0 order)
                            (lambda (x y)
                              (< (third x) (third y))))))
             (list (loop for x in sz collect (third x))
                   (loop for x in sz collect (second x))
                   (loop for x in sz collect (fourth x))))))
    (lambda (xx)
      (destructuring-bind ((x1 x2) (fx1 fx2) (dfx1 dfx2)) (get-order xx 2)
        (+ (* (/_ (* (+ 1 (* 2 (- xx x1))) (- xx x2) (- xx x2))
                  (* (- x1 x2) (- x1 x2))) fx1)
           (* (/_ (* (- 1 (* 2 (- xx x2))) (- xx x1) (- xx x1))
                  (* (- x1 x2) (- x1 x2))) fx2)
           (/_ (* (- xx x1) (- xx x2) (- xx x2) dfx1)
               (* (- x1 x2) (- x1 x2)))
           (/_ (* (- xx x2) (- xx x1) (- xx x1) dfx2)
               (* (- x1 x2) (- x1 x2))))))))

(defun tridiagonal-solver (a b c r)
  (labels ((beta (a b c)
             (labels ((rec (a b pc pbeta &optional res)
                        (if (null a)
                            res
                            (if (zerop pbeta)
                                (make-list (length a) :initial-element 0)
                                (let ((v (- (car b) (* (car a) (car pc) (/ 1 pbeta)))))
                                  (rec (cdr a) (cdr b) (cdr pc) v (cons v res)))))))
               (append (rec (cdr a) (cdr b) c (car b)) (list (car b)))))
           (rho (r a beta)
             (labels ((rec (r a pbeta prho &optional res)
                        (if (null a)
                            res
                            (let ((v (- (car r) (* (car a) prho (/_ 1 (car pbeta))))))
                              (rec (cdr r) (cdr a) (cdr pbeta) v (cons v res))))))
               (append (rec (cdr r) (cdr a) beta (car r)) (list (car r)))))
           (backsub (beta rho c px)
             (when beta
               (let ((v (/_ (- (car rho) (* (car c) px)) (car beta))))
                 (cons v (backsub (cdr beta) (cdr rho) (cdr c) v))))))
    (if (find 0 b)
        (values b nil)
        (let ((beta (beta a b c)))
          (if (find 0 beta)
              (values beta nil)
              (let ((rho (rho r a (reverse beta))))
                (values
                 (reverse (backsub beta
                                   rho
                                   (reverse c)
                                   (/_ (car rho) (car beta))))
                 t)))))))

(defun clamped-cubic-spline (x f fp1 fpn)
  (flet ((h (j) (- (nth (1+ j) x) (nth j x)))
         (p (j) (- (nth (1+ j) f) (nth j f))))
    (let* ((n (1- (length x)))
           (ac (loop for i from 0 to (1- n) collect (h i))))
      (tridiagonal-solver
       (cons 0 ac)
       (append (list (* 2 (h 0)))
               (loop for i from 1 to (- n 1) collect
                    (* 2 (- (nth (1+ i) x) (nth (1- i) x))))
               (list (* 2 (h (1- n)))))
       (append ac (cons 0 nil))
       (append (list (* 6 (- (/_ (p 0) (h 0)) fp1)))
               (loop for i from 1 to (- n 1) collect
                    (- (* 6 (/_ (p i) (h i))) (* 6 (/ (p (1- i)) (h (1- i))))))
               (list (* -6 (- (/_ (p (1- n)) (h (1- n))) fpn))))))))

(defun natural-cubic-spline (x f)
  (flet ((h (j) (- (nth (1+ j) x) (nth j x)))
         (p (j) (- (nth (1+ j) f) (nth j f))))
    (let* ((n (1- (length x)))
           (ac (cons 0 (append (loop for i from 1 to (- n 1) collect (h i)) (list 0)))))
      (tridiagonal-solver
       ac
       (append (list 1)
               (loop for i from 1 to (- n 1) collect
                    (* 2 (- (nth (1+ i) x) (nth (1- i) x))))
               (list 1))
       ac
       (append (list 0)
               (loop for i from 1 to (- n 1) collect
                    (- (* 6 (/_ (p i) (h i))) (* 6 (/_ (p (1- i)) (h (1- i))))))
               (list 0))))))

(defun cubic-spline (x f &optional fp1 fpn)
  "Create a lambda to approximate f(xx) using a cubic spline through a
list of values and their function evaluations. The derivates at the
endpoints may optionally be provided, in which case a clamped cubic
spline will be used. Otherwise a natural cubic spline will be used."
  (flet ((h (j) (- (nth (1+ j) x) (nth j x)))
         (p (j) (- (nth (1+ j) f) (nth j f))))
    (let ((second (if fpn (clamped-cubic-spline x f fp1 fpn) (natural-cubic-spline x f))))
      (lambda (xv)
        (destructuring-bind ((j j+1) (xj xj+1) (pj pj+1) (qj qj+1)) (bracket/3 xv x f second)
          (declare (ignore j+1 xj+1 pj+1))
          (+ pj
             (* (- xv xj)
                (- (/_ (p j) (h j))
                   (* (h j) qj+1 1/6)
                   (* (h j) qj 1/3)))
             (* (- xv xj) (- xv xj) 1/2 qj)
             (* (- xv xj) (- xv xj) (- xv xj) (- qj+1 qj) (/_ 1/6 (h j)))))))))

(defun cubic-spline-prime (x f &optional fp1 fpn)
  "Create a lambda to approximate the derivative of f(xx) using a
cubic spline through a list of values and their function
evaluations. The derivates at the endpoints may optionally be
provided, in which case a clamped cubic spline will be used. Otherwise
a natural cubic spline will be used."
  (flet ((h (j) (- (nth (1+ j) x) (nth j x)))
         (p (j) (- (nth (1+ j) f) (nth j f))))
    (let ((second (if fp1 (clamped-cubic-spline x f fp1 fpn) (natural-cubic-spline x f))))
      (lambda (xv)
        (destructuring-bind ((j j+1) (xj xj+1) (pj pj+1) (qj qj+1)) (bracket/3 xv x f second)
          (declare (ignore j+1 xj+1 pj+1 pj))
          (+ (/_ (p j) (h j))
             (* 1d0 -1/6 (h j) qj+1)
             (* -1/3 (h j) qj)
             (* qj (- xv xj))
             (* (- xv xj) (- xv xj) (/_ (- qj+1 qj) (* 2 (h j))))))))))

(defun first-derivative/2 (f &key (step +step+))
  "Create a lambda to approximate the first derivative of `F` at any point. `F` will be called
2 times per approximation."
  (lambda (x)
    (/ (- (funcall f (+ x step)) (funcall f (- x step)))
       (* 2d0 step))))

(defun first-derivative/4 (f &key (step +step+))
  "Create a lambda to approximate the first derivative of `F` at any point. `F` will be called
4 times per approximation."
  (lambda (x)
    (/ (+ (funcall f (- x step step))
          (* -8d0 (funcall f (- x step)))
          (* 8d0 (funcall f (+ x step)))
          (* -1d0 (funcall f (+ x step step))))
       (* 12d0 step))))

(defun second-derivative/3 (f &key (step +step+))
  "Create a lambda to approximate the second derivative of `F` at any point. `F` will be called
3 times per approximation."
  (lambda (x)
    (/ (+ (funcall f x) (* -2d0 (funcall f (+ x step))) (funcall f (+ x step step)))
       (* step step))))

(defun second-derivative/5 (f &key (step +step+))
  "Create a lambda to approximate the second derivative of `F` at any point. `F` will be called
5 times per approximation."
  (lambda (x)
    (/ (+ (* -1d0 (funcall f (- x step step)))
          (* 16d0 (funcall f (- x step)))
          (* -30d0 (funcall f x))
          (* 16d0 (funcall f (+ x step)))
          (* -1d0 (funcall f (+ x step step))))
       (* 12d0 step step))))

(defun richardson-extrapolation (f df h tol limit)
  (labels ((di+1 (i dih di2h)
             (/_ (- (* (expt 2d0 (+ i i)) dih) di2h)
                 (1- (expt 2d0 (+ i i)))))
           (rec (x &optional (i 1) (h h) (dh (funcall (funcall df f :step h) x)))
             (let ((x (coerce x 'double-float)))
               (if (>= i limit)
                   (values dh nil i)
                   (let ((nv (di+1 (1+ i) (funcall (funcall df f :step (/ h 2d0)) x) dh)))
                     (if (< (abs (/_ (- nv dh) nv)) tol)
                         (values nv t i)
                         (rec x (1+ i) (/ h 2d0) nv)))))))
    #'rec))

(defun first-derivative (f &key (step +step+) (tolerance +tolerance+) (limit +limit+))
  "Create a lambda to approximate the first derivative of a function
at any point. Richardson extrapolation will be used to converge on the
solution down to a specified tolerance."
  (richardson-extrapolation f #'first-derivative/2 step tolerance limit))

(defmacro first-derivative+ ((x &key (step +step+) (tolerance +tolerance+) (limit +limit+)) f)
  "Create a lambda to approximate the first derivative of a function
at any point. Richardson extrapolation will be used to converge on the
solution down to a specified tolerance."
  `(richardson-extrapolation #'(lambda (,x) ,f) #'first-derivative/2 ,step ,tolerance ,limit))

(defun second-derivative (f &key (step +step+) (tolerance +tolerance+) (limit +limit+))
  "Create a lambda to approximate the second derivative of a function
at any point. Richardson extrapolation will be used to converge on the
solution down to a specified tolerance."
  (richardson-extrapolation f #'second-derivative/3 step tolerance limit))

(defmacro second-derivative+ ((x &key (step +step+) (tolerance +tolerance+) (limit +limit+)) f)
  "Create a lambda to approximate the second derivative of a function
at any point. Richardson extrapolation will be used to converge on the
solution down to a specified tolerance."
  `(richardson-extrapolation #'(lambda (,x) ,f) #'second-derivative/3 ,step ,tolerance ,limit))
